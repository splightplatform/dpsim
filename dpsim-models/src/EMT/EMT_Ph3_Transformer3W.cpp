#include <dpsim-models/EMT/EMT_Ph3_Transformer3W.h>
#include <dpsim-models/MathUtils.h>

using namespace CPS;

using Winding = CPS::EMT::Ph3::Transformer3W::Winding;
using CPS::Base::Ph3::AllWindings3W;
using CPS::Base::Ph3::WindingConnection; 

namespace {

/// Short tag used in sub-component names and log lines
const char *windingTag(Winding w) {
  switch (w) {
  case Winding::Primary:
    return "p";
  case Winding::Secondary:
    return "s";
  default:
    return "t";
  }
}

const std::unordered_map<int, PhaseType> phaseMap = {
  {0, PhaseType::A}, 
  {1, PhaseType::B}, 
  {2, PhaseType::C}, 
};

UInt idx(Winding w) { return static_cast<UInt>(w); }

} // namespace

/// General

EMT::Ph3::Transformer3W::Transformer3W(String uid, String name,
                                       Logger::Level logLevel,
                                       Bool withResistiveLosses)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel),
      Base::Ph3::Transformer3W(mAttributes),
      mWithResistiveLosses(withResistiveLosses) {

  mPhaseType = PhaseType::ABC;
  setVirtualNodeNumber(virtualNodeCount(withResistiveLosses));
  setTerminalNumber(NumWindings);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);

  // One 3-phase (3x1) interface column per winding: row = phase (A,B,C),
  // column = winding (P,S,T)
  **mIntfVoltage = Matrix::Zero(3, NumWindings);
  **mIntfCurrent = Matrix::Zero(3, NumWindings);
}

UInt EMT::Ph3::Transformer3W::virtualNodeCount(Bool withResistiveLosses) {
  // TODO: only wye right now
  return 1 + (withResistiveLosses ? NumWindings : 0) + 2 * (NumWindings - 1);
}

void EMT::Ph3::Transformer3W::assignVirtualNodeSlots() {
  UInt next = 0;
  mVnStar = next++;

  for (auto w : AllWindings3W) {
    if (mWithResistiveLosses)
      mVnMid[idx(w)] = next++;

    if (!isReferenceWinding(w)) {
      mVnPreIdeal[idx(w)] = next++;
      mVnCurrent[idx(w)] = next++;
    }
  }

  // check we assign the correct number of virtual nodes
  if (next != mNumVirtualNodes) {
    SPDLOG_LOGGER_ERROR(mSLog,
                        "Virtual node slot assignment used {} slots but {} "
                        "were allocated; every allocated slot must be used or "
                        "the system matrix acquires an all-zero row.",
                        next, mNumVirtualNodes);
    throw InvalidArgumentException();
  }

  SPDLOG_LOGGER_INFO(mSLog,
                     "Virtual node layout (reference winding = {}): star = {}",
                     windingTag(mReferenceWinding), mVnStar);
  for (auto w : AllWindings3W)
    SPDLOG_LOGGER_INFO(mSLog, "  winding {}: mid = {}, preIdeal = {}, i = {}",
                       windingTag(w), mVnMid[idx(w)], mVnPreIdeal[idx(w)],
                       mVnCurrent[idx(w)]);
}

SimPowerComp<Real>::Ptr EMT::Ph3::Transformer3W::clone(String name) {
  auto copy = Transformer3W::make(name, mLogLevel);

  // Named with a suffix so they do not shadow the same-named accessors on
  // the base class, which the loop below still needs to call
  std::array<Real, NumWindings> nomVoltageArr{{0., 0., 0.}};
  std::array<Real, NumWindings> ratedPowerArr{{0., 0., 0.}};
  std::array<Real, NumWindings> ratioAbsArr{{0., 0., 0.}};
  std::array<Real, NumWindings> ratioPhaseArr{{0., 0., 0.}};
  std::array<Matrix, NumWindings> resistanceArr{
      {Matrix::Zero(3, 3), Matrix::Zero(3, 3), Matrix::Zero(3, 3)}};
  std::array<Matrix, NumWindings> inductanceArr{
      {Matrix::Zero(3, 3), Matrix::Zero(3, 3), Matrix::Zero(3, 3)}};

  for (auto w : AllWindings3W) {
    nomVoltageArr[idx(w)] = nominalVoltage(w);
    ratedPowerArr[idx(w)] = ratedPower(w);
    ratioAbsArr[idx(w)] = std::abs(ratio(w));
    ratioPhaseArr[idx(w)] = std::arg(ratio(w));
    resistanceArr[idx(w)] = resistance(w);
    inductanceArr[idx(w)] = inductance(w);
  }

  copy->setReferenceWinding(mReferenceWinding);
  for (auto w : AllWindings3W)
    copy->setWindingConnection(w, windingConnection(w), phaseShift(w));
  copy->setParameters(nomVoltageArr, ratedPowerArr, ratioAbsArr,
                      ratioPhaseArr, resistanceArr, inductanceArr);
  return copy;
}

void EMT::Ph3::Transformer3W::connectWinding(Winding w) {
  if (windingConnection(w) == WindingConnection::Delta){
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "Transformer3W {}: connectWinding() is not implemented yet for "
        "delta-leg topologies.",
        this->name());
    throw CPS::Exception();
  }
  const UInt i = idx(w);
  const String tag = windingTag(w);

  //create components 
  mSubInductor[i] = std::make_shared<EMT::Ph3::Inductor>(
      **mUID + "_ind_" + tag, **mName + "_ind_" + tag, Logger::Level::off);
  mSubInductor[i]->setParameters(inductance(w));
  addMNASubComponent(mSubInductor[i],
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  if (mWithResistiveLosses) {
    mSubResistor[i] = std::make_shared<EMT::Ph3::Resistor>(
        **mUID + "_res_" + tag, **mName + "_res_" + tag, Logger::Level::off);
    mSubResistor[i]->setParameters(resistance(w));
    addMNASubComponent(mSubResistor[i],
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }
  // for wye winding use same topology as SP_Ph1_Transformer3W.cpp
  auto starEnd = mVirtualNodes[mVnStar]; 
  auto outerEnd = isReferenceWinding(w) ? node(i) : mVirtualNodes[mVnPreIdeal[i]]; 

  if (mWithResistiveLosses){
    auto mid = mVirtualNodes[mVnMid[i]]; 
    if (isReferenceWinding(w)) { 
      mSubResistor[i]->connect({outerEnd, mid}); 
      mSubInductor[i]->connect({mid, starEnd}); 
    } else {
      mSubResistor[i]->connect({starEnd, mid}); 
      mSubInductor[i]->connect({mid, outerEnd}); 
    }
  } else {
    if (isReferenceWinding(w))
      mSubInductor[i]->connect({outerEnd, starEnd}); 
    else 
      mSubInductor[i]->connect({starEnd, outerEnd}); 
  }

}

void EMT::Ph3::Transformer3W::createSubComponents() {
  if (mSubCompCreated)
    return;
  mSubCompCreated = true;

  resolveReferenceWinding();
  assignVirtualNodeSlots(); // throws until Phase 3.0 topology is derived

  for (auto w : AllWindings3W)
    connectWinding(w);

  Bool snubbersEnabled =
      (mBehaviour == TopologicalPowerComp::Behaviour::Initialization ||
       mBehaviour == TopologicalPowerComp::Behaviour::MNASimulation);

  if (snubbersEnabled){
    for (auto w : AllWindings3W){
      if (ratedPower(w) <= 0) {
        SPDLOG_LOGGER_WARN(
            mSLog,
            "Rated power of winding {} is {} [VA]; snubbers disabled for this "
            "transformer (cannot be sized off non-positive power).",
            windingTag(w), ratedPower(w));
        snubbersEnabled = false;
        break;
      }
    }
  }

  if (!snubbersEnabled)
    return; 

  // build and connect snubbers (size later)
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    const String tag = windingTag(w);

    mSubSnubResistor[i] = std::make_shared<EMT::Ph3::Resistor>(
        **mName + "_snub_res_" + tag, mLogLevel);
    mSubSnubResistor[i]->connect({node(i), EMT::SimNode::GND});
    addMNASubComponent(mSubSnubResistor[i],
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

    // No capacitor on the reference (highest-voltage) terminal: the
    // two-winding model found an HV-side snubber capacitor made conditioning
    // worse
    if (isReferenceWinding(w))
      continue;

    mSubSnubCapacitor[i] = std::make_shared<EMT::Ph3::Capacitor>(
        **mName + "_snub_cap_" + tag, mLogLevel);
    mSubSnubCapacitor[i]->connect({node(i), EMT::SimNode::GND});
    addMNASubComponent(mSubSnubCapacitor[i],
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }

  // Diagnostic: star point -> GND, testing whether the star point is
  mSubSnubResistorStar = std::make_shared<EMT::Ph3::Resistor>(
      **mName + "_snub_res_star", mLogLevel);
  mSubSnubResistorStar->connect({mVirtualNodes[mVnStar], EMT::SimNode::GND});
  addMNASubComponent(mSubSnubResistorStar,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
}

void EMT::Ph3::Transformer3W::initializeParentFromNodesAndTerminals(
    Real frequency) {

  /* WORKAROUND SimPowerComp<Real>::initialize(), called earlier via
  SystemTopology::addComponent(), unconditionally resizes mIntfVoltage/
  mIntfCurrent to (3, mNumFreqs) -- (3, 1) for a normal single-frequency
  run
  This breakes the (3, NumWindings) shape set in the constructor
  This class uses columns to mean "winding," not "frequency" like the
  base class assumes
  To fix a memory error I'm re-establishing the real shape here, since nothing
  re-triggers initialize() after this point 
  TODO: a more permanent fix is need in the future. I've noted it */
  **mIntfVoltage = Matrix::Zero(3, NumWindings);
  **mIntfCurrent = Matrix::Zero(3, NumWindings);

  Real mNominalOmega = 2. * PI * frequency;

  for (auto w : AllWindings3W)
    SPDLOG_LOGGER_INFO(
        mSLog, "Winding {}: R = {} [Ohm], L = {} [H], X = {} [Ohm]",
        windingTag(w), Logger::matrixToString(resistance(w)),
        Logger::matrixToString(inductance(w)),
        Logger::matrixToString(mNominalOmega * inductance(w)));

  // size snubbers
  if (mSubSnubResistor[0]) {
    for (auto w : AllWindings3W) {
      const UInt i = idx(w);
      const Real pSnub = P_SNUB_TRANSFORMER * ratedPower(w);
      const Real qSnub = Q_SNUB_TRANSFORMER * ratedPower(w);
      const Real vNom = nominalVoltage(w);

      const Real snubberResistance = std::pow(std::abs(vNom), 2) / pSnub; 
      mSnubberResistance[i] = Math::singlePhaseParameterToThreePhase(snubberResistance); 
      mSubSnubResistor[i]->setParameters(mSnubberResistance[i]);
      SPDLOG_LOGGER_INFO(mSLog,
                         "Snubber resistance at winding {} ({}) = {} [Ohm]",
                         windingTag(w), node(i)->name(),
                         Logger::matrixToString(mSnubberResistance[i]));

      if (!mSubSnubCapacitor[i])
        continue;

      const Real snubberCapacitance = qSnub / std::pow(std::abs(vNom), 2) / mNominalOmega;
      mSnubberCapacitance[i] = Math::singlePhaseParameterToThreePhase(snubberCapacitance); 
      mSubSnubCapacitor[i]->setParameters(mSnubberCapacitance[i]);
      SPDLOG_LOGGER_INFO(mSLog,
                         "Snubber capacitance at winding {} ({}) = {} [F]",
                         windingTag(w), node(i)->name(),
                         Logger::matrixToString(mSnubberCapacitance[i]));
    }

    // Diagnostic star-point snubber. Sized off the reference winding's
    // rating
    const Real pSnubStar = P_SNUB_TRANSFORMER * ratedPower(mReferenceWinding);
    const Real vNomStar = nominalVoltage(mReferenceWinding);
    const Real snubberResistanceStar =
        std::pow(std::abs(vNomStar), 2) / pSnubStar;
    mSnubberResistanceStar =
        Math::singlePhaseParameterToThreePhase(snubberResistanceStar);
    mSubSnubResistorStar->setParameters(mSnubberResistanceStar);
    SPDLOG_LOGGER_INFO(mSLog, "Snubber resistance at star point = {} [Ohm]",
                       Logger::matrixToString(mSnubberResistanceStar));
  }

  std::array<MatrixComp, NumWindings> referredVoltage; // 3x1 vector (3ph) for each winding
  std::array<MatrixComp, NumWindings> windingImpedance;  // 3x3 vector (3ph w/ phase mutual inductance) for each winding
  MatrixComp sumCurrent = MatrixComp::Zero(3, 1);
  MatrixComp sumAdmittance = MatrixComp::Zero(3,3); 

  // calculate and star voltage and pre-ideal voltages from load-flow values
  for (auto w : AllWindings3W){
    if (windingConnection(w) == WindingConnection::Delta){
      SPDLOG_LOGGER_ERROR(
          mSLog,
          "Transformer3W {}: connectWinding() is not implemented yet for "
          "delta-leg topologies.",
          this->name());
      throw CPS::Exception();
    }
    
    const UInt i = idx(w);
    MatrixComp vABC = MatrixComp::Zero(3,1); 
    vABC(0,0) = initialSingleVoltage(i); 
    vABC(1,0) = initialSingleVoltage(i) * SHIFT_TO_PHASE_B; 
    vABC(2,0) = initialSingleVoltage(i) * SHIFT_TO_PHASE_C; 
    referredVoltage[i] = ratio(w) * vABC; 
    windingImpedance[i] = resistance(w).cast<Complex>() + Complex(0,mNominalOmega) * inductance(w).cast<Complex>();
  
    if (!isReferenceWinding(w))
      mVirtualNodes[mVnPreIdeal[i]]->setInitialVoltage(referredVoltage[i]);
  
    sumCurrent += windingImpedance[i].inverse() * referredVoltage[i];
    sumAdmittance += windingImpedance[i].inverse();
  }

  // calculate and set starVoltage
  const MatrixComp starVoltage = sumAdmittance.inverse() * sumCurrent;
  mVirtualNodes[mVnStar]->setInitialVoltage(starVoltage);

  SPDLOG_LOGGER_INFO(mSLog, "Star point initial voltage (3ph):\n{:s}",
                     Logger::matrixCompToString(starVoltage));

  // calculate and set midpoint voltages
  for (auto w : AllWindings3W){
    const UInt i = idx(w);
    const MatrixComp branchCurrent = windingImpedance[i].inverse() * (starVoltage - referredVoltage[i]);

    SPDLOG_LOGGER_INFO(mSLog, "Winding {}: branch current (3ph):\n{:s}",
                       windingTag(w), Logger::matrixCompToString(branchCurrent));

    // no midpoint node if not using resistors
    if (mWithResistiveLosses) {
      // if reference winding the voltage winding to midpoint is over inductor
      // else its over the resistor
      const MatrixComp compImpedance =
        isReferenceWinding(w) ? MatrixComp(Complex(0, mNominalOmega) * inductance(w).cast<Complex>())
                              : MatrixComp(resistance(w).cast<Complex>());
      const MatrixComp midVoltage = starVoltage - (compImpedance * branchCurrent);
      mVirtualNodes[mVnMid[i]]->setInitialVoltage(midVoltage);

      SPDLOG_LOGGER_INFO(mSLog, "Winding {}: midpoint initial voltage (3ph):\n{:s}",
                         windingTag(w), Logger::matrixCompToString(midVoltage));
    }

    // seed MNA current and voltages
    MatrixComp outerVoltage = MatrixComp::Zero(3,1); 
    Complex outerVoltageA = isReferenceWinding(w) ? initialSingleVoltage(i) : referredVoltage[i](0,0); 
    outerVoltage(0,0) = outerVoltageA; 
    outerVoltage(1,0) = outerVoltageA * SHIFT_TO_PHASE_B; 
    outerVoltage(2,0) = outerVoltageA * SHIFT_TO_PHASE_C; 

    for (auto &[phase_idx, phase] : phaseMap){
      (**mIntfVoltage)(phase_idx, i) = RMS3PH_TO_PEAK1PH * (outerVoltage(phase_idx,0) - starVoltage(phase_idx)).real(); 
      (**mIntfCurrent)(phase_idx, i) = RMS3PH_TO_PEAK1PH * (isReferenceWinding(w) ? -branchCurrent(phase_idx).real() :
                                      branchCurrent(phase_idx).real());
    }

    SPDLOG_LOGGER_INFO(mSLog,
                       "Winding {}: seeded intfVoltage = {:s}, intfCurrent = {:s}",
                       windingTag(w),
                       Logger::matrixToString((**mIntfVoltage).col(i)),
                       Logger::matrixToString((**mIntfCurrent).col(i)));
  }

  SPDLOG_LOGGER_INFO(mSLog,
                     "\n--- Initialization from powerflow ---"
                     "\nTerminal 0 (p) voltage: {:s}"
                     "\nTerminal 1 (s) voltage: {:s}"
                     "\nTerminal 2 (t) voltage: {:s}"
                     "\n--- Initialization from powerflow finished ---",
                     Logger::phasorToString(initialSingleVoltage(0)),
                     Logger::phasorToString(initialSingleVoltage(1)),
                     Logger::phasorToString(initialSingleVoltage(2)));
}

// MNA Section

void EMT::Ph3::Transformer3W::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  (void)omega;
  (void)leftVector;
  for (auto w : AllWindings3W)
    SPDLOG_LOGGER_INFO(mSLog,
                       "Terminal {} ({}) connected to {:s} = sim node {:d}",
                       idx(w), windingTag(w),
                       mTerminals[idx(w)]->node()->name(),
                       mTerminals[idx(w)]->node()->matrixNodeIndex());

  // Diagnostic: is the terminal snubber actually providing meaningful
  // damping relative to the branch's own trapezoidal conductance, or is it
  // negligible next to it?
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    const Matrix branchGeq = timeStep / 2. * inductance(w).inverse();
    SPDLOG_LOGGER_INFO(mSLog, "Winding {}: branch Geq (timeStep/2 * L^-1) = {} [S]",
                       windingTag(w), Logger::matrixToString(branchGeq));

    if (!mSubSnubResistor[i]) {
      SPDLOG_LOGGER_INFO(mSLog, "Winding {}: snubbers disabled, no comparison",
                         windingTag(w));
      continue;
    }

    const Matrix snubGeq = mSnubberResistance[i].inverse();
    const Real ratio = snubGeq(0, 0) / branchGeq(0, 0);
    SPDLOG_LOGGER_INFO(
        mSLog,
        "Winding {}: snubber Geq = {} [S], snubGeq/branchGeq = {} "
        "(near 0 means the snubber is negligible next to the branch)",
        windingTag(w), Logger::matrixToString(snubGeq), ratio);
  }

  // Diagnostic: same comparison for the star point
  if (mSubSnubResistorStar) {
    const Matrix starSnubGeq = mSnubberResistanceStar.inverse();
    SPDLOG_LOGGER_INFO(mSLog, "Star point: snubber Geq = {} [S]",
                       Logger::matrixToString(starSnubGeq));
    for (auto w : AllWindings3W) {
      const Matrix branchGeq = timeStep / 2. * inductance(w).inverse();
      const Real ratio = starSnubGeq(0, 0) / branchGeq(0, 0);
      SPDLOG_LOGGER_INFO(
          mSLog, "Star point: starSnubGeq/branchGeq[{}] = {}",
          windingTag(w), ratio);
    }
  }
}

void EMT::Ph3::Transformer3W::stampIdealTransformer(
    SparseMatrixRow &systemMatrix, Winding w) {

  const UInt i = idx(w); 
  if (isReferenceWinding(w))
      return; 
  if (!terminalNotGrounded(i))
      return; 
  if (windingConnection(w) == WindingConnection::Delta){
    SPDLOG_LOGGER_ERROR(
      mSLog,
      "Transformer3W {}: stampIdealTransformer() is not implemented yet for "
      "delta-leg topologies.",
      this->name());
      throw CPS::Exception();
    }

  // The ideal transformer imposes v_preIdeal = T * v_terminal and i_terminal = -T * i_branch
  //
  //          preIdeal   terminal   current
  //   preIdeal    .         .        -1
  //   terminal    .         .         T
  //   current     1        -T         0
  // This is stamped once per phase 
    
  Real T = ratio(w).real(); 
  for (auto &[phase_idx, phase] : phaseMap){

    const UInt preIdeal = mVirtualNodes[mVnPreIdeal[i]]->matrixNodeIndex(phase); 
    const UInt current = mVirtualNodes[mVnCurrent[i]]->matrixNodeIndex(phase);
    const UInt terminal = matrixNodeIndex(i, phase_idx); 
    
    Math::setMatrixElement(systemMatrix, preIdeal, current, -1.0); 
    Math::setMatrixElement(systemMatrix, current, preIdeal, 1.0); 
    Math::setMatrixElement(systemMatrix, terminal, current, T); 
    Math::setMatrixElement(systemMatrix, current, terminal, -T); 


    SPDLOG_LOGGER_INFO(mSLog,
                       "Ideal transformer for winding {}: preIdeal = {}, "
                       "terminal = {}, current unknown = {}",
                       windingTag(w), preIdeal, terminal, current);
  }
  
}

void EMT::Ph3::Transformer3W::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  for (auto w : AllWindings3W)
    stampIdealTransformer(systemMatrix, w); 

  for (auto subcomp : mSubComponents)
    if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
      mnasubcomp->mnaApplySystemMatrixStamp(systemMatrix);
}

void EMT::Ph3::Transformer3W::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::Transformer3W::mnaParentPreStep(Real time, Int timeStepCount) {
  (void)time;
  (void)timeStepCount;
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::Transformer3W::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph3::Transformer3W::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  (void)time;
  (void)timeStepCount;
  this->mnaCompUpdateVoltage(**leftVector);
  this->mnaCompUpdateCurrent(**leftVector);
}

void EMT::Ph3::Transformer3W::mnaCompUpdateCurrent(const Matrix &leftVector) {

  for (auto w : AllWindings3W) {
    if (windingConnection(w) == WindingConnection::Delta){
      SPDLOG_LOGGER_ERROR(
          mSLog,
          "Transformer3W {}: mnaCompUpdateCurrent() is not implemented yet for "
          "delta-leg topologies.",
          this->name());
      throw CPS::Exception();
    }
    const UInt i = idx(w);
    (**mIntfCurrent).col(i) = mSubInductor[i]->intfCurrent();
  }
}

void EMT::Ph3::Transformer3W::mnaCompUpdateVoltage(const Matrix &leftVector) {
  
  for (auto w : AllWindings3W) {
    if (windingConnection(w) == WindingConnection::Delta){
      SPDLOG_LOGGER_ERROR(
          mSLog,
          "Transformer3W {}: mnaCompUpdateVoltage() is not implemented yet for "
          "delta-leg topologies.",
          this->name());
      throw CPS::Exception();
    }
    const UInt i = idx(w);
    for (auto &[phase_idx, phase] : phaseMap){
      const Real vStar = Math::realFromVectorElement(
          leftVector, mVirtualNodes[mVnStar]->matrixNodeIndex(phase));
      const Real vOuter =
          isReferenceWinding(w)
              ? Math::realFromVectorElement(leftVector, matrixNodeIndex(i, phase_idx))
              : Math::realFromVectorElement(
                    leftVector, mVirtualNodes[mVnPreIdeal[i]]->matrixNodeIndex(phase));
      (**mIntfVoltage)(phase_idx, i) = vOuter - vStar;
    }
  }
}
