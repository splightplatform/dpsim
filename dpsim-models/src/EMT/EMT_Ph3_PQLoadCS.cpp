
#include <dpsim-models/EMT/EMT_Ph3_PQLoadCS.h>

using namespace CPS;

EMT::Ph3::PQLoadCS::PQLoadCS(String uid, String name, Logger::Level logLevel)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel),
      mActivePower(mAttributes->create<Real>("P", 0)),
      mReactivePower(mAttributes->create<Real>("Q", 0)),
      mNomVoltage(mAttributes->create<Real>("V_nom")) {
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(1);
  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);
}

EMT::Ph3::PQLoadCS::PQLoadCS(String name, Logger::Level logLevel)
    : PQLoadCS(name, name, logLevel) {}

EMT::Ph3::PQLoadCS::PQLoadCS(String name, Real activePower,
                             Real reactivePower, Real volt,
                             Logger::Level logLevel)
    : PQLoadCS(name, logLevel) {
  setParameters(activePower, reactivePower, volt);
}

void EMT::Ph3::PQLoadCS::setParameters(Real activePower, Real reactivePower,
                                       Real nomVolt) {
  **mActivePower = activePower;
  **mReactivePower = reactivePower;
  **mNomVoltage = nomVolt;
  mParametersSet = true;

  SPDLOG_LOGGER_INFO(mSLog,
                     "\nActive Power [W]: {:f}"
                     "\nReactive Power [VAr]: {:f}",
                     activePower, reactivePower);
  SPDLOG_LOGGER_INFO(mSLog, "Nominal Voltage={} [V]", nomVolt);
}

SimPowerComp<Real>::Ptr EMT::Ph3::PQLoadCS::clone(String name) {
  auto copy = PQLoadCS::make(name, mLogLevel);
  copy->setParameters(**mActivePower, **mReactivePower, **mNomVoltage);
  return copy;
}

void EMT::Ph3::PQLoadCS::createSubComponents() {
  if (mSubCompCreated)
    return;
  mSubCompCreated = true;

  mSubCurrentSource =
      std::make_shared<EMT::Ph3::CurrentSource>(**mName + "_cs", mLogLevel);
  // A positive power should result in a positive current to ground.
  mSubCurrentSource->connect({mTerminals[0]->node(), SimNode::GND});
  addMNASubComponent(mSubCurrentSource,
                     MNA_SUBCOMP_TASK_ORDER::TASK_AFTER_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
}

void EMT::Ph3::PQLoadCS::initializeParentFromNodesAndTerminals(
    Real frequency) {
  // Read power from terminals here, not in createSubComponents(): power-flow results aren't available that early.
  if (!mParametersSet) {
    **mActivePower = mTerminals[0]->singleActivePower();
    **mReactivePower = mTerminals[0]->singleReactivePower();
    **mNomVoltage = std::abs(mTerminals[0]->initialSingleVoltage());
  }

  mVoltageRefAngle = Math::phase(mTerminals[0]->initialSingleVoltage());

  // Balanced three-phase current at nominal voltage: the actual power
  // delivered will only match the setpoint exactly while the terminal
  // voltage stays close to V_nom<mVoltageRefAngle.
  Complex powerPerPhase = Complex(**mActivePower, **mReactivePower) / 3.;
  Complex vPhase = Math::polar(**mNomVoltage / sqrt(3.), mVoltageRefAngle);
  Complex iPhase =
      (vPhase != Complex(0, 0)) ? std::conj(powerPerPhase / vPhase) : 0.;

  mSubCurrentSource->setParameters(
      CPS::Math::singlePhaseVariableToThreePhase(iPhase), frequency);
  mSubCurrentSource->initializeFromNodesAndTerminals(frequency);
  updateIntfValues();

  SPDLOG_LOGGER_INFO(mSLog,
                     "\n--- Initialization from powerflow ---"
                     "\nVoltage across: {:s}"
                     "\nCurrent: {:s}"
                     "\nTerminal 0 voltage: {:s}"
                     "\nActive Power: {:f}"
                     "\nReactive Power: {:f}"
                     "\n--- Initialization from powerflow finished ---",
                     Logger::matrixToString(**mIntfVoltage),
                     Logger::matrixToString(**mIntfCurrent),
                     Logger::phasorToString(initialSingleVoltage(0)),
                     **mActivePower, **mReactivePower);
}

void EMT::Ph3::PQLoadCS::updateSetPoint() {
  // Phase-to-neutral RMS magnitude from the instantaneous, previous-step
  // solved terminal voltage. For a balanced three-phase set the sum of the
  // squared instantaneous samples is time-invariant (independent of the
  // actual phase angle), so this is exact and needs no PLL:
  //   va^2+vb^2+vc^2 = 3/2 * Vpeak^2  =>  Vrms = sqrt((va^2+vb^2+vc^2)/3)
  Real vMagRms = sqrt(((**mIntfVoltage)(0, 0) * (**mIntfVoltage)(0, 0) +
                       (**mIntfVoltage)(1, 0) * (**mIntfVoltage)(1, 0) +
                       (**mIntfVoltage)(2, 0) * (**mIntfVoltage)(2, 0)) /
                      3.);
  // Not measured yet (first step, before any MNA solve has run): fall back
  // to the nominal voltage used at initialization.
  if (vMagRms < 1e-6)
    vMagRms = **mNomVoltage / sqrt(3.);

  Complex powerPerPhase = Complex(**mActivePower, **mReactivePower) / 3.;
  Complex vPhase = Math::polar(vMagRms, mVoltageRefAngle);
  Complex iPhase =
      (vPhase != Complex(0, 0)) ? std::conj(powerPerPhase / vPhase) : 0.;

  **mSubCurrentSource->mCurrentRef =
      CPS::Math::singlePhaseVariableToThreePhase(iPhase);

  SPDLOG_LOGGER_DEBUG(mSLog,
                      "\n--- update set points ---"
                      "\nActive Power: {:f}"
                      "\nReactive Power: {:f}"
                      "\nCurrent: {:s}",
                      **mActivePower, **mReactivePower,
                      Logger::matrixCompToString(**mSubCurrentSource->mCurrentRef));
}

void EMT::Ph3::PQLoadCS::updateIntfValues() {
  **mIntfCurrent = mSubCurrentSource->intfCurrent();
  **mIntfVoltage = mSubCurrentSource->intfVoltage();
}

void EMT::Ph3::PQLoadCS::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  attributeDependencies.push_back(mActivePower);
  attributeDependencies.push_back(mReactivePower);
  attributeDependencies.push_back(mNomVoltage);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::PQLoadCS::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfCurrent);
  modifiedAttributes.push_back(mIntfVoltage);
}

void EMT::Ph3::PQLoadCS::mnaParentPreStep(Real time, Int timeStepCount) {
  updateSetPoint();
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::PQLoadCS::mnaParentPostStep(Real time, Int timeStepCount,
                                           Attribute<Matrix>::Ptr &leftVector) {
  updateIntfValues();
}
