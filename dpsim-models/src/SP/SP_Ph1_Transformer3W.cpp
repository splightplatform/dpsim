#include <dpsim-models/MathUtils.h>
#include <dpsim-models/SP/SP_Ph1_Transformer3W.h>

using namespace CPS;

using Winding = CPS::Base::Ph1::Winding3W;
using CPS::Base::Ph1::AllWindings3W;

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

UInt idx(Winding w) { return static_cast<UInt>(w); }
} // namespace

/// General

SP::Ph1::Transformer3W::Transformer3W(String uid, String name,
                                      Logger::Level logLevel,
                                      Bool withResistiveLosses)
    : Base::Ph1::Transformer3W(mAttributes),
      CompositePowerComp<Complex>(uid, name, true, true, logLevel),
      mWithResistiveLosses(withResistiveLosses),
      mBaseVoltage(mAttributes->create<Real>("base_Voltage")),
      mCurrent(mAttributes->create<MatrixComp>("current_vector")),
      mActivePowerBranch(mAttributes->create<Matrix>("p_branch_vector")),
      mReactivePowerBranch(mAttributes->create<Matrix>("q_branch_vector")),
      mActivePowerInjection(mAttributes->create<Matrix>("p_inj_vector")),
      mReactivePowerInjection(mAttributes->create<Matrix>("q_inj_vector")) {

  setVirtualNodeNumber(virtualNodeCount(withResistiveLosses));
  setTerminalNumber(NumWindings);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);

  /// One interface entry per winding (3 x 1)
  **mIntfVoltage = MatrixComp::Zero(NumWindings, 1);
  **mIntfCurrent = MatrixComp::Zero(NumWindings, 1);
  **mCurrent = MatrixComp::Zero(NumWindings, 1);
  **mActivePowerBranch = Matrix::Zero(NumWindings, 1);
  **mReactivePowerBranch = Matrix::Zero(NumWindings, 1);
  **mActivePowerInjection = Matrix::Zero(NumWindings, 1);
  **mReactivePowerInjection = Matrix::Zero(NumWindings, 1);
}

UInt SP::Ph1::Transformer3W::virtualNodeCount(Bool withResistiveLosses) {
  // star point -> 1
  // + R-L midpoint per winding, if withResistveLosses -> 3
  // + pre-ideal xfmr node and current unknown, per non-ref winding -> 2 * 2
  return 1 + (withResistiveLosses ? NumWindings : 0) + 2 * (NumWindings - 1);
}

void SP::Ph1::Transformer3W::assignVirtualNodeSlots() {
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

SimPowerComp<Complex>::Ptr SP::Ph1::Transformer3W::clone(String name) {
  auto copy = Transformer3W::make(name, mLogLevel);

  // Named with a suffix so they do not shadow the same-named accessors on the
  // base class, which the loop below still needs to call
  std::array<Real, NumWindings> nomVoltageArr{{0., 0., 0.}};
  std::array<Real, NumWindings> ratedPowerArr{{0., 0., 0.}};
  std::array<Real, NumWindings> ratioAbsArr{{0., 0., 0.}};
  std::array<Real, NumWindings> ratioPhaseArr{{0., 0., 0.}};
  std::array<Real, NumWindings> resistanceArr{{0., 0., 0.}};
  std::array<Real, NumWindings> inductanceArr{{0., 0., 0.}};

  for (auto w : AllWindings3W) {
    nomVoltageArr[idx(w)] = nominalVoltage(w);
    ratedPowerArr[idx(w)] = ratedPower(w);
    ratioAbsArr[idx(w)] = std::abs(ratio(w));
    ratioPhaseArr[idx(w)] = std::arg(ratio(w));
    resistanceArr[idx(w)] = resistance(w);
    inductanceArr[idx(w)] = inductance(w);
  }

  copy->setReferenceWinding(mReferenceWinding);
  copy->setParameters(nomVoltageArr, ratedPowerArr, ratioAbsArr, ratioPhaseArr,
                      resistanceArr, inductanceArr);
  return copy;
}

void SP::Ph1::Transformer3W::connectWinding(Winding w) {
  const UInt i = idx(w);
  const String tag = windingTag(w);

  // create components
  mSubInductor[i] = std::make_shared<SP::Ph1::Inductor>(
      **mUID + "_ind_" + tag, **mName + "_ind_" + tag, Logger::Level::off);
  mSubInductor[i]->setParameters(inductance(w));
  addMNASubComponent(mSubInductor[i],
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  if (mWithResistiveLosses) {
    mSubResistor[i] = std::make_shared<SP::Ph1::Resistor>(
        **mUID + "_res_" + tag, **mName + "_res_" + tag, Logger::Level::off);
    mSubResistor[i]->setParameters(resistance(w));
    addMNASubComponent(mSubResistor[i],
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }

  // The reference winding runs terminal -> star with no ideal transformer,
  // because the star impedances are already on its side. Every other winding
  // runs star -> ideal transformer -> terminal
  // Sub-component orientation follows:
  //
  //   reference:      node(w) -[R]- mid -[L]- star
  //   non-reference:  star -[R]- mid -[L]- preIdeal = ideal = node(w)
  //
  // inductorcurrent is positive flowing towards the star point on the reference
  // winding and away from it on the others (assumed conventional flow)

  // wire up components
  auto starEnd = mVirtualNodes[mVnStar];
  auto outerEnd = isReferenceWinding(w) // before ideal xfmr
                      ? node(i)
                      : mVirtualNodes[mVnPreIdeal[i]];

  if (mWithResistiveLosses) {
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

void SP::Ph1::Transformer3W::createSubComponents() {
  if (mSubCompCreated)
    return;
  mSubCompCreated = true;

  resolveReferenceWinding();
  assignVirtualNodeSlots();

  for (auto w : AllWindings3W)
    connectWinding(w);

  // Snubbers are sized off the rated power, so rating required
  Bool snubbersEnabled =
      (mBehaviour == TopologicalPowerComp::Behaviour::Initialization ||
       mBehaviour == TopologicalPowerComp::Behaviour::MNASimulation);

  // ensure we can safely make the snubbers
  if (snubbersEnabled) {
    for (auto w : AllWindings3W) {
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

  // build and create snubbers (DPsim jargon to help simulation stability), 
  // size later in initializeParentFromNodesAndTerminals()
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    const String tag = windingTag(w);

    mSubSnubResistor[i] = std::make_shared<SP::Ph1::Resistor>(
        **mName + "_snub_res_" + tag, mLogLevel);
    mSubSnubResistor[i]->connect({node(i), SP::SimNode::GND});
    addMNASubComponent(mSubSnubResistor[i],
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

    // No capacitor on the reference (highest-voltage) terminal: the
    // two-winding model found an HV-side snubber capacitor made conditioning
    // worse, so we use that here
    if (isReferenceWinding(w))
      continue;

    mSubSnubCapacitor[i] = std::make_shared<SP::Ph1::Capacitor>(
        **mName + "_snub_cap_" + tag, mLogLevel);
    mSubSnubCapacitor[i]->connect({node(i), SP::SimNode::GND});
    addMNASubComponent(mSubSnubCapacitor[i],
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }
}

void SP::Ph1::Transformer3W::initializeParentFromNodesAndTerminals(
    Real frequency) {
  mNominalOmega = 2. * PI * frequency;

  for (auto w : AllWindings3W)
    SPDLOG_LOGGER_INFO(
        mSLog, "Winding {}: R = {} [Ohm], L = {} [H], X = {} [Ohm]",
        windingTag(w), resistance(w), inductance(w),
        mNominalOmega * inductance(w));

  //  Snubber sizing 
  if (mSubSnubResistor[0]) {
    for (auto w : AllWindings3W) {
      const UInt i = idx(w);
      const Real pSnub = P_SNUB_TRANSFORMER * ratedPower(w);
      const Real qSnub = Q_SNUB_TRANSFORMER * ratedPower(w);
      const Real vNom = nominalVoltage(w);

      mSnubberResistance[i] = std::pow(std::abs(vNom), 2) / pSnub;
      mSubSnubResistor[i]->setParameters(mSnubberResistance[i]);
      mSubSnubResistor[i]->setBaseVoltage(vNom);
      SPDLOG_LOGGER_INFO(mSLog,
                         "Snubber resistance at winding {} ({}) = {} [Ohm]",
                         windingTag(w), node(i)->name(),
                         Logger::realToString(mSnubberResistance[i]));

      if (!mSubSnubCapacitor[i])
        continue;

      mSnubberCapacitance[i] =
          qSnub / std::pow(std::abs(vNom), 2) / mNominalOmega;
      mSubSnubCapacitor[i]->setParameters(mSnubberCapacitance[i]);
      mSubSnubCapacitor[i]->setBaseVoltage(vNom);
      SPDLOG_LOGGER_INFO(mSLog,
                         "Snubber capacitance at winding {} ({}) = {} [F]",
                         windingTag(w), node(i)->name(),
                         Logger::realToString(mSnubberCapacitance[i]));
    }
  }

  std::array<Complex, NumWindings> referredVoltage{};
  std::array<Complex, NumWindings> windingImpedance{};
  Complex sumCurrent = 0;
  Complex sumAdmittance = 0;

  // calculate and set pre-ideal xfmr voltages
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    referredVoltage[i] = ratio(w) * initialSingleVoltage(i);
    windingImpedance[i] = {resistance(w), mNominalOmega * inductance(w)};

    if (!isReferenceWinding(w))
      mVirtualNodes[mVnPreIdeal[i]]->setInitialVoltage(referredVoltage[i]);

    sumCurrent += referredVoltage[i] / windingImpedance[i];
    sumAdmittance += 1.0 / windingImpedance[i];
  }

  // calculate and set starVoltage
  const Complex starVoltage = sumCurrent / sumAdmittance;
  mVirtualNodes[mVnStar]->setInitialVoltage(starVoltage);

  // calculate and set mid-voltages
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    const Complex branchCurrent = (starVoltage - referredVoltage[i]) / windingImpedance[i];

    if (mWithResistiveLosses) {
      const Complex compImpedance =
          isReferenceWinding(w) ? Complex(0, mNominalOmega * inductance(w))
                                : resistance(w);
      const Complex midVoltage = starVoltage - branchCurrent * compImpedance;
      mVirtualNodes[mVnMid[i]]->setInitialVoltage(midVoltage);
    }

    // seed MNA current and voltages
    const Complex outerVoltage =
        isReferenceWinding(w) ? initialSingleVoltage(i) : referredVoltage[i];
    (**mIntfVoltage)(i, 0) = outerVoltage - starVoltage;
    (**mIntfCurrent)(i, 0) =
        isReferenceWinding(w) ? -branchCurrent : branchCurrent; // DPsim convention
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

// #### Powerflow section ####

void SP::Ph1::Transformer3W::setBaseVoltage(Real baseVoltage) {
  **mBaseVoltage = baseVoltage;
}

void SP::Ph1::Transformer3W::calculatePerUnitParameters(Real baseApparentPower,
                                                        Real baseOmega) {
  SPDLOG_LOGGER_INFO(mSLog, "#### Calculate Per Unit Parameters for {}",
                     **mName);

  mBaseApparentPower = baseApparentPower;
  mBaseOmega = baseOmega;
  mBaseImpedance = **mBaseVoltage * **mBaseVoltage / mBaseApparentPower;
  mBaseAdmittance = 1.0 / mBaseImpedance;
  mBaseInductance = mBaseImpedance / mBaseOmega;
  // I_base = (S_3ph / 3) / (V_LL / sqrt(3))
  mBaseCurrent = baseApparentPower / (**mBaseVoltage * sqrt(3));

  SPDLOG_LOGGER_INFO(mSLog, "Base Power = {} [VA]  Base Omega = {} [1/s]",
                     baseApparentPower, baseOmega);
  SPDLOG_LOGGER_INFO(mSLog, "Base Voltage = {} [V]  Base Impedance = {} [Ohm]",
                     **mBaseVoltage, mBaseImpedance);

  for (auto w : AllWindings3W){
    const UInt i = idx(w); 
    mLeakagePerUnit[i] = Complex{resistance(w) / mBaseImpedance, inductance(w) / mBaseInductance};
    mRatioPerUnit[i] = ratio(w) * nominalVoltage(w) / (**mBaseVoltage); 
  }

  for (auto w : AllWindings3W)
    SPDLOG_LOGGER_INFO(mSLog, "Winding {}: Z = {} [pu], tap = {} [pu]",
                       windingTag(w),
                       Logger::complexToString(mLeakagePerUnit[idx(w)]),
                       Logger::complexToString(mRatioPerUnit[idx(w)]));

  // Sub-component per-unit values, so their own PF stamps land on the same
  // base as ours
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    if (mSubSnubResistor[i])
      mSubSnubResistor[i]->calculatePerUnitParameters(mBaseApparentPower);
    if (mSubSnubCapacitor[i])
      mSubSnubCapacitor[i]->calculatePerUnitParameters(mBaseApparentPower);
  }
}

void SP::Ph1::Transformer3W::pfApplyAdmittanceMatrixStamp(
    SparseMatrixCompRow &Y) {

  MatrixComp _mY_element = MatrixComp::Zero(NumWindings, NumWindings); // just Kron reduced w/o sum factor
  mY_element = MatrixComp::Zero(NumWindings, NumWindings); // Kron reduced w/ sum factor and tap
  Complex Yp = 1.0 / mLeakagePerUnit[idx(Winding::Primary)]; 
  Complex Ys = 1.0 / mLeakagePerUnit[idx(Winding::Secondary)]; 
  Complex Yt = 1.0 / mLeakagePerUnit[idx(Winding::Tertiary)]; 
  Complex sumY = Yp + Ys + Yt; 

  _mY_element(0, 0) = Yp * (Ys+Yt); 
  _mY_element(0, 1) = -Yp * Ys; 
  _mY_element(0, 2) = -Yp * Yt; 
  _mY_element(1, 0) = -Yp * Ys; 
  _mY_element(1, 1) = Ys * (Yp + Yt); 
  _mY_element(1, 2) = -Ys * Yt; 
  _mY_element(2, 0) = -Yp * Yt; 
  _mY_element(2, 1) = -Ys * Yt; 
  _mY_element(2, 2) = Yt * (Yp+Ys); 

  for (UInt i = 0 ; i < NumWindings; i++){
    for (UInt j = 0 ; j < NumWindings; j++){
    // see Grainger 9.6
    mY_element(i,j) = _mY_element(i,j) / sumY * std::conj(mRatioPerUnit[i]) * mRatioPerUnit[j]; 
    }
  }

  for (UInt i = 0; i < NumWindings; i++)
    for (UInt j = 0; j < NumWindings; j++)
      if (!Math::isFinite(mY_element.coeff(i, j))) {
        SPDLOG_LOGGER_ERROR(
            mSLog,
            "Transformer3W {}: non-finite per-unit admittance {} in element "
            "Y({},{})",
            this->name(), Logger::complexToString(mY_element.coeff(i, j)), i,
            j);
        throw InvalidArgumentException();
      }

  for (UInt i = 0; i < NumWindings; i++)
    for (UInt j = 0; j < NumWindings; j++)
      Y.coeffRef(this->matrixNodeIndex(i), this->matrixNodeIndex(j)) +=
          mY_element.coeff(i, j);

  SPDLOG_LOGGER_INFO(mSLog, "#### Y matrix stamping: {}", mY_element);

  // apply snubs to the admittance matrix
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    if (mSubSnubResistor[i])
      mSubSnubResistor[i]->pfApplyAdmittanceMatrixStamp(Y);
    if (mSubSnubCapacitor[i])
      mSubSnubCapacitor[i]->pfApplyAdmittanceMatrixStamp(Y);
  }
}

void SP::Ph1::Transformer3W::updateBranchFlow(VectorComp &current,
                                              VectorComp &powerflow) {
  **mCurrent = current * mBaseCurrent;
  **mActivePowerBranch = powerflow.real() * mBaseApparentPower;
  **mReactivePowerBranch = powerflow.imag() * mBaseApparentPower;
}

void SP::Ph1::Transformer3W::storeNodalInjection(Winding w,
                                                 Complex powerInjection) {
  const UInt i = idx(w);
  (**mActivePowerInjection)(i, 0) =
      std::real(powerInjection) * mBaseApparentPower;
  (**mReactivePowerInjection)(i, 0) =
      std::imag(powerInjection) * mBaseApparentPower;
}

MatrixComp SP::Ph1::Transformer3W::Y_element() { return mY_element; }

// MNA Section  

void SP::Ph1::Transformer3W::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  for (auto w : AllWindings3W)
    SPDLOG_LOGGER_INFO(mSLog, "Terminal {} ({}) connected to {:s} = sim node {:d}",
                       idx(w), windingTag(w),
                       mTerminals[idx(w)]->node()->name(),
                       mTerminals[idx(w)]->node()->matrixNodeIndex());
}

void SP::Ph1::Transformer3W::stampIdealTransformer(
    SparseMatrixRow &systemMatrix, Winding w) {
  if (isReferenceWinding(w))
    return;

  const UInt i = idx(w);
  const UInt preIdeal = mVirtualNodes[mVnPreIdeal[i]]->matrixNodeIndex();
  const UInt current = mVirtualNodes[mVnCurrent[i]]->matrixNodeIndex();
  const UInt terminal = matrixNodeIndex(i);

  // The ideal transformer imposes v_preIdeal = T * v_terminal and i_terminal = -conj(T) * i_branch
  //
  //          preIdeal   terminal   current
  //   preIdeal    .         .        -1
  //   terminal    .         .      coeffTerminalCurrent
  //   current     1   coeffCurrentTerminal    0

  Complex coeffCurrentTerminal = -ratio(w);
  Complex coeffTerminalCurrent = std::conj(ratio(w));

  Math::setMatrixElement(systemMatrix, preIdeal, current, Complex(-1.0, 0));
  Math::setMatrixElement(systemMatrix, current, preIdeal, Complex(1.0, 0));
  Math::setMatrixElement(systemMatrix, terminal, current,
                         coeffTerminalCurrent);
  Math::setMatrixElement(systemMatrix, current, terminal,
                         coeffCurrentTerminal);

  SPDLOG_LOGGER_INFO(mSLog,
                     "Ideal transformer for winding {}: preIdeal = {}, "
                     "terminal = {}, current unknown = {}",
                     windingTag(w), preIdeal, terminal, current);
}

void SP::Ph1::Transformer3W::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  for (auto w : AllWindings3W)
    stampIdealTransformer(systemMatrix, w);

  for (auto subcomp : mSubComponents)
    if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
      mnasubcomp->mnaApplySystemMatrixStamp(systemMatrix);
}

void SP::Ph1::Transformer3W::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mRightVector);
}

void SP::Ph1::Transformer3W::mnaParentPreStep(Real time, Int timeStepCount) {
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void SP::Ph1::Transformer3W::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void SP::Ph1::Transformer3W::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  this->mnaUpdateVoltage(**leftVector);
  this->mnaUpdateCurrent(**leftVector);
}

void SP::Ph1::Transformer3W::mnaCompUpdateCurrent(const Matrix &leftVector) {
  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    (**mIntfCurrent)(i, 0) = mSubInductor[i]->intfCurrent()(0, 0);
  }
}

void SP::Ph1::Transformer3W::mnaCompUpdateVoltage(const Matrix &leftVector) {
  // Voltage across each winding's branch: terminal side minus star side
  const Complex vStar = Math::complexFromVectorElement(
      leftVector, mVirtualNodes[mVnStar]->matrixNodeIndex());

  for (auto w : AllWindings3W) {
    const UInt i = idx(w);
    const Complex vOuter =
        isReferenceWinding(w)
            ? Math::complexFromVectorElement(leftVector, matrixNodeIndex(i))
            : Math::complexFromVectorElement(
                  leftVector, mVirtualNodes[mVnPreIdeal[i]]->matrixNodeIndex());
    (**mIntfVoltage)(i, 0) = vOuter - vStar;
  }
}
