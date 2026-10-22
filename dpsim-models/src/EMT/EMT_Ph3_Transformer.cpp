/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>

using namespace CPS;

EMT::Ph3::Transformer::Transformer(String uid, String name,
                                   Logger::Level logLevel,
                                   Bool withResistiveLosses)
    : Base::Ph3::Transformer(mAttributes),
      CompositePowerComp<Real>(uid, name, true, true, logLevel) {
  mPhaseType = PhaseType::ABC;
  if (withResistiveLosses)
    setVirtualNodeNumber(3);
  else
    setVirtualNodeNumber(2);

  setTerminalNumber(2);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(1, 1);
}

/// DEPRECATED: Delete method
SimPowerComp<Real>::Ptr EMT::Ph3::Transformer::clone(String name) {
  auto copy = Transformer::make(name, mLogLevel);
  copy->setParameters(mNominalVoltageEnd1, mNominalVoltageEnd2, mRatedPower,
                      std::abs(**mRatio), std::arg(**mRatio), mResistance,
                      mInductance);
  copy->setZeroSequencePassthrough(mZeroSequencePassthrough); 
  return copy;
}

void EMT::Ph3::Transformer::setParameters(Real nomVoltageEnd1,
                                          Real nomVoltageEnd2, Real ratedPower,
                                          Real ratioAbs, Real ratioPhase,
                                          Matrix resistance,
                                          Matrix inductance) {

  Base::Ph3::Transformer::setParameters(nomVoltageEnd1, nomVoltageEnd2,
                                        ratedPower, ratioAbs, ratioPhase,
                                        resistance, inductance);

  SPDLOG_LOGGER_INFO(
      mSLog, "Nominal Voltage End 1 = {} [V] Nominal Voltage End 2 = {} [V]",
      mNominalVoltageEnd1, mNominalVoltageEnd2);
  SPDLOG_LOGGER_INFO(mSLog, "Rated Apparent Power  = {} [VA]", mRatedPower);
  SPDLOG_LOGGER_INFO(mSLog, "Tap Ratio = {} [ ] Phase Shift = {} [deg]",
                     std::abs(**mRatio), std::arg(**mRatio));

  mParametersSet = true;
}

void EMT::Ph3::Transformer::createSubComponents() {
  if (mSubCompCreated)
    return;
  mSubCompCreated = true;

  // Switch terminals so terminal 0 is always the higher-voltage side.
  if (Math::abs(**mRatio) < 1.) {
    // refer the series impedance to the new terminal-0 side: Z' = Z / |a|²
    Real ratioNormSq = std::norm(**mRatio); 
    mResistance /= ratioNormSq;
    mInductance /= ratioNormSq;
    **mRatio = 1. / **mRatio;
    std::shared_ptr<SimTerminal<Real>> tmp = mTerminals[0];
    mTerminals[0] = mTerminals[1];
    mTerminals[1] = tmp;
    Real tmpVolt = mNominalVoltageEnd1;
    mNominalVoltageEnd1 = mNominalVoltageEnd2;
    mNominalVoltageEnd2 = tmpVolt;
    SPDLOG_LOGGER_INFO(mSLog, "Switching terminals to have first terminal at "
                              "higher voltage side. Updated parameters: ");
    SPDLOG_LOGGER_INFO(
        mSLog, "Nominal Voltage End 1 = {} [V] Nominal Voltage End 2 = {} [V]",
        mNominalVoltageEnd1, mNominalVoltageEnd2);
    SPDLOG_LOGGER_INFO(mSLog, "Tap Ratio = {} [ ] Phase Shift = {} [deg]",
                       std::abs(**mRatio), std::arg(**mRatio));
    SPDLOG_LOGGER_INFO(mSLog,
                       "Resistance (re-referred to new terminal 0) = {} [Ohm]",
                       Logger::matrixToString(mResistance));
    SPDLOG_LOGGER_INFO(mSLog,
                       "Inductance (re-referred to new terminal 0) = {} [H]",
                       Logger::matrixToString(mInductance));
  }

  // Create T matrix: implements T_ij = N · ( k0/3 + (2/3)·cos(θ + (j − i)·2π/3) ), i, j = 0..2
  Real N = std::abs(**mRatio); 
  Real theta = std::arg(**mRatio); 
  Real k0 = mZeroSequencePassthrough ? 1. : 0.;

  if (theta == 0. && k0 == 1.){
    mTurnsRatio = N * Matrix::Identity(3,3);
  } else{
    mTurnsRatio = Matrix(3,3); 
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        mTurnsRatio(i,j) = N * (k0 / 3 + (2./3.) * std::cos(theta + (j - i) * (2*PI/3) ) );
      }
    }
  }
  SPDLOG_LOGGER_INFO(mSLog,
                     "Turns ratio matrix T (N = {}, phase shift = {} [deg], "
                     "k0 = {}) = {}",
                     N, theta * 180. / PI, k0,
                     Logger::matrixToString(mTurnsRatio));

  // Create series sub components
  mSubInductor =
      std::make_shared<EMT::Ph3::Inductor>(**mName + "_ind", mLogLevel);
  mSubInductor->setParameters(mInductance);
  addMNASubComponent(mSubInductor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  if (mNumVirtualNodes == 3) {
    mSubResistor =
        std::make_shared<EMT::Ph3::Resistor>(**mName + "_res", mLogLevel);
    addMNASubComponent(mSubResistor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
    mSubResistor->setParameters(mResistance);
    mSubResistor->connect({node(0), mVirtualNodes[2]});
    mSubInductor->connect({mVirtualNodes[2], mVirtualNodes[0]});
  } else {
    mSubInductor->connect({node(0), mVirtualNodes[0]});
  }

  // Create parallel sub components (three-phase power)
  Real pSnub = P_SNUB_TRANSFORMER * mRatedPower;

  // A snubber conductance is added on the higher voltage side
  Real snubberResistance1 = std::pow(std::abs(mNominalVoltageEnd1), 2) / pSnub;
  mSnubberResistance1 =
      Math::singlePhaseParameterToThreePhase(snubberResistance1);
  mSubSnubResistor1 =
      std::make_shared<EMT::Ph3::Resistor>(**mName + "_snub_res1", mLogLevel);
  mSubSnubResistor1->setParameters(mSnubberResistance1);
  mSubSnubResistor1->connect({node(0), EMT::SimNode::GND});
  SPDLOG_LOGGER_INFO(
      mSLog,
      "Snubber Resistance 1 (connected to higher voltage side {}) = {} [Ohm]",
      node(0)->name(), Logger::matrixToString(mSnubberResistance1));
  addMNASubComponent(mSubSnubResistor1,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  // A snubber conductance is added on the lower voltage side
  Real snubberResistance2 = std::pow(std::abs(mNominalVoltageEnd2), 2) / pSnub;
  mSnubberResistance2 =
      Math::singlePhaseParameterToThreePhase(snubberResistance2);
  mSubSnubResistor2 =
      std::make_shared<EMT::Ph3::Resistor>(**mName + "_snub_res2", mLogLevel);
  mSubSnubResistor2->setParameters(mSnubberResistance2);
  mSubSnubResistor2->connect({node(1), EMT::SimNode::GND});
  SPDLOG_LOGGER_INFO(
      mSLog,
      "Snubber Resistance 2 (connected to lower voltage side {}) = {} [Ohm]",
      node(1)->name(), Logger::matrixToString(mSnubberResistance2));
  addMNASubComponent(mSubSnubResistor2,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  // LV-side snubber capacitor created here; its omega-dependent value is set in initializeParentFromNodesAndTerminals().
  mSubSnubCapacitor2 =
      std::make_shared<EMT::Ph3::Capacitor>(**mName + "_snub_cap2", mLogLevel);
  mSubSnubCapacitor2->connect({node(1), EMT::SimNode::GND});
  addMNASubComponent(mSubSnubCapacitor2,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
}

void EMT::Ph3::Transformer::initializeParentFromNodesAndTerminals(
    Real frequency) {
  // Set initial voltage of virtual node in between
  mVirtualNodes[0]->setInitialVoltage(initialSingleVoltage(1) * **mRatio);

  // Static calculations from load flow data
  Real omega = 2. * PI * frequency;

  Real qSnub = Q_SNUB_TRANSFORMER * mRatedPower;

  // // A snubber capacitance is added to higher voltage side (not used as capacitor at high voltage side made it worse)
  // Real snubberCapacitance1 = qSnub / std::pow(std::abs(mNominalVoltageEnd1),2) / omega;
  // mSnubberCapacitance1 = Math::singlePhaseParameterToThreePhase*(snubberCapacitance1);
  // mSubSnubCapacitor1 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_snub_cap1", mLogLevel);
  // mSubSnubCapacitor1->setParameters(mSnubberCapacitance1);
  // mSubSnubCapacitor1->connect({ node(0), EMT::SimNode::GND });
  // SPDLOG_LOGGER_INFO(mSLog, "Snubber Capacitance 1 (connected to higher voltage side {}) = \n{} [F] \n ", node(0)->name(), Logger::matrixToString(mSnubberCapacitance1));
  // mSubComponents.push_back(mSubSnubCapacitor1);

  Real snubberCapacitance2 =
      qSnub / std::pow(std::abs(mNominalVoltageEnd2), 2) / omega;
  mSnubberCapacitance2 =
      Math::singlePhaseParameterToThreePhase(snubberCapacitance2);
  mSubSnubCapacitor2->setParameters(mSnubberCapacitance2);
  SPDLOG_LOGGER_INFO(
      mSLog,
      "Snubber Capacitance 2 (connected to lower voltage side {}) = {} [F]",
      node(1)->name(), Logger::matrixToString(mSnubberCapacitance2));
  MatrixComp impedance = MatrixComp::Zero(3, 3);
  impedance << Complex(mResistance(0, 0), omega * mInductance(0, 0)),
      Complex(mResistance(0, 1), omega * mInductance(0, 1)),
      Complex(mResistance(0, 2), omega * mInductance(0, 2)),
      Complex(mResistance(1, 0), omega * mInductance(1, 0)),
      Complex(mResistance(1, 1), omega * mInductance(1, 1)),
      Complex(mResistance(1, 2), omega * mInductance(1, 2)),
      Complex(mResistance(2, 0), omega * mInductance(2, 0)),
      Complex(mResistance(2, 1), omega * mInductance(2, 1)),
      Complex(mResistance(2, 2), omega * mInductance(2, 2));

  SPDLOG_LOGGER_INFO(mSLog,
                     "Resistance (referred to higher voltage side) = {} [Ohm]",
                     Logger::matrixToString(mResistance));
  SPDLOG_LOGGER_INFO(mSLog,
                     "Inductance (referred to higher voltage side) = {} [H]",
                     Logger::matrixToString(mInductance));
  SPDLOG_LOGGER_INFO(mSLog,
                     "Reactance (referred to higher voltage side) = {} [Ohm]",
                     Logger::matrixToString(omega * mInductance));

  MatrixComp vInitABC = MatrixComp::Zero(3, 1);
  vInitABC(0, 0) =
      RMS3PH_TO_PEAK1PH *
      (mVirtualNodes[0]->initialSingleVoltage() - initialSingleVoltage(0));
  vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
  vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;

  MatrixComp iInit = impedance.inverse() * vInitABC;
  **mIntfCurrent = iInit.real();
  **mIntfVoltage = vInitABC.real();

  if (mNumVirtualNodes == 3)
    mVirtualNodes[2]->setInitialVoltage(initialSingleVoltage(0));

  SPDLOG_LOGGER_INFO(
      mSLog,
      "\n--- Initialization from powerflow ---"
      "\nVoltage across: {:s}"
      "\nCurrent: {:s}"
      "\nTerminal 0 voltage: {:s}"
      "\nTerminal 1 voltage: {:s}"
      "\nVirtual Node 1 voltage: {:s}"
      "\n--- Initialization from powerflow finished ---",
      Logger::matrixToString(**mIntfVoltage),
      Logger::matrixToString(**mIntfCurrent),
      Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
      Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)),
      Logger::phasorToString(RMS3PH_TO_PEAK1PH *
                             mVirtualNodes[0]->initialSingleVoltage()));
}

void EMT::Ph3::Transformer::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  SPDLOG_LOGGER_INFO(
      mSLog,
      "\nTerminal 0 connected to {:s} = sim node {:d}"
      "\nTerminal 1 connected to {:s} = sim node {:d}",
      mTerminals[0]->node()->name(), mTerminals[0]->node()->matrixNodeIndex(),
      mTerminals[1]->node()->name(), mTerminals[1]->node()->matrixNodeIndex());
}

void EMT::Ph3::Transformer::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  // Ideal transformer equations
  if (terminalNotGrounded(0)) {
    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::A),
        mVirtualNodes[1]->matrixNodeIndex(PhaseType::A), -1.);
    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::B),
        mVirtualNodes[1]->matrixNodeIndex(PhaseType::B), -1.);
    Math::setMatrixElement(
        systemMatrix, mVirtualNodes[0]->matrixNodeIndex(PhaseType::C),
        mVirtualNodes[1]->matrixNodeIndex(PhaseType::C), -1.);

    Math::setMatrixElement(systemMatrix,
                           mVirtualNodes[1]->matrixNodeIndex(PhaseType::A),
                           mVirtualNodes[0]->matrixNodeIndex(PhaseType::A), 1.);
    Math::setMatrixElement(systemMatrix,
                           mVirtualNodes[1]->matrixNodeIndex(PhaseType::B),
                           mVirtualNodes[0]->matrixNodeIndex(PhaseType::B), 1.);
    Math::setMatrixElement(systemMatrix,
                           mVirtualNodes[1]->matrixNodeIndex(PhaseType::C),
                           mVirtualNodes[0]->matrixNodeIndex(PhaseType::C), 1.);
  }
  if (terminalNotGrounded(1)) {
    // loop through T matrix updating the system matrix 
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        if (mTurnsRatio(j, i) != 0.)
          Math::setMatrixElement(systemMatrix, matrixNodeIndex(1, i),
                                mVirtualNodes[1]->matrixNodeIndex(static_cast<PhaseType>(j)),
                                mTurnsRatio(j, i));
        if (mTurnsRatio(i,j) != 0.)
          Math::setMatrixElement(systemMatrix, 
                                mVirtualNodes[1]->matrixNodeIndex(static_cast<PhaseType>(i)),
                                matrixNodeIndex(1,j), -mTurnsRatio(i,j));
      }
    }
  }

  // Add subcomps to system matrix
  for (auto subcomp : mSubComponents)
    if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
      mnasubcomp->mnaApplySystemMatrixStamp(systemMatrix);

  if (terminalNotGrounded(0)) {
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0, 0)),
                       mVirtualNodes[0]->matrixNodeIndex(PhaseType::A),
                       mVirtualNodes[1]->matrixNodeIndex(PhaseType::A));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0, 0)),
                       mVirtualNodes[0]->matrixNodeIndex(PhaseType::B),
                       mVirtualNodes[1]->matrixNodeIndex(PhaseType::B));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(-1.0, 0)),
                       mVirtualNodes[0]->matrixNodeIndex(PhaseType::C),
                       mVirtualNodes[1]->matrixNodeIndex(PhaseType::C));

    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(1.0, 0)),
                       mVirtualNodes[1]->matrixNodeIndex(PhaseType::A),
                       mVirtualNodes[0]->matrixNodeIndex(PhaseType::A));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(1.0, 0)),
                       mVirtualNodes[1]->matrixNodeIndex(PhaseType::B),
                       mVirtualNodes[0]->matrixNodeIndex(PhaseType::B));
    SPDLOG_LOGGER_INFO(mSLog, "Add {:s} to system at ({:d},{:d})",
                       Logger::complexToString(Complex(1.0, 0)),
                       mVirtualNodes[1]->matrixNodeIndex(PhaseType::C),
                       mVirtualNodes[0]->matrixNodeIndex(PhaseType::C));
  }
}

void EMT::Ph3::Transformer::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::Transformer::mnaParentPreStep(Real time, Int timeStepCount) {
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::Transformer::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph3::Transformer::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateVoltage(**leftVector);
  mnaCompUpdateCurrent(**leftVector);
}

void EMT::Ph3::Transformer::mnaCompUpdateCurrent(const Matrix &leftVector) {
  **mIntfCurrent = mSubInductor->intfCurrent();
}

void EMT::Ph3::Transformer::mnaCompUpdateVoltage(const Matrix &leftVector) {
  // v1 - v0
  **mIntfVoltage = Matrix::Zero(3, 1);
  if (terminalNotGrounded(1)) {
    (**mIntfVoltage)(0, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
    (**mIntfVoltage)(1, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 1));
    (**mIntfVoltage)(2, 0) =
        Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 2));
  }
  if (terminalNotGrounded(0)) {
    (**mIntfVoltage)(0, 0) =
        (**mIntfVoltage)(0, 0) -
        Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
    (**mIntfVoltage)(1, 0) =
        (**mIntfVoltage)(1, 0) -
        Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
    (**mIntfVoltage)(2, 0) =
        (**mIntfVoltage)(2, 0) -
        Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
  }
}
