#pragma once
#include <array>

#include <dpsim-models/Base/Base_Ph1_Transformer3W.h>
#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/SP/SP_Ph1_Capacitor.h>
#include <dpsim-models/SP/SP_Ph1_Inductor.h>
#include <dpsim-models/SP/SP_Ph1_Resistor.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Solver/PFSolverInterfaceBranch.h>

namespace CPS {
namespace SP {
namespace Ph1 {

// Three-winding transformer, static phasor.
//
// Star (wye) equivalent: three series R+L branches meeting at a fictitious
// star point (Chap 2.8 of Grainger) with an ideal transformer 
// between each non-reference winding's branch and its terminal
// The reference winding connects to the star point
// directly, since the star impedances are expressed on its side already
//
// Two representations of 3W xfmr
//
//   - Power flow (pfApplyAdmittanceMatrixStamp) eliminates the star point
//     entirely and hands the solver a dense 3x3 bus admittance block.
//   - MNA (mnaCompApplySystemMatrixStamp) keeps the star point as a real
//     node and stamps each ideal transformer as its own extra row/column.
//
class Transformer3W : public CompositePowerComp<Complex>,
                      public Base::Ph1::Transformer3W,
                      public SharedFactory<Transformer3W>,
                      public PFSolverInterfaceBranch {

public:
  using Winding = Base::Ph1::Winding3W;

private:
  // unused virtual node
  static constexpr UInt UNUSED_VN = static_cast<UInt>(-1);

  // DPsim Sub-components (one set per winding)

  // Series resistance of each star branch 
  // only used when mWithResistiveLosses is true (default for 3W xfmr) 
  std::array<std::shared_ptr<SP::Ph1::Resistor>, NumWindings> mSubResistor;
  // Series inductance of each star branch
  std::array<std::shared_ptr<SP::Ph1::Inductor>, NumWindings> mSubInductor;

  // Snubber conductance at each terminal, for numerical damping
  // copied from 2W xfmr case 
  std::array<std::shared_ptr<SP::Ph1::Resistor>, NumWindings> mSubSnubResistor;
  // Snubber capacitance at the non-reference terminals only
  std::array<std::shared_ptr<SP::Ph1::Capacitor>, NumWindings> mSubSnubCapacitor;

  // Snubber values for logging [Ohm] / [F].
  std::array<Real, NumWindings> mSnubberResistance{{0., 0., 0.}};
  std::array<Real, NumWindings> mSnubberCapacitance{{0., 0., 0.}};

  // Virtual-node slot assignment 
  // The layout depends on which winding is the reference and on whether the
  // resistances are split out, so the slots are assigned once in
  // createSubComponents() rather than hard-coded
  // The star point. Always allocated, always slot 0.
  UInt mVnStar = 0;
  // Node between R_i and L_i. UNUSED_VN when !mWithResistiveLosses.
  std::array<UInt, NumWindings> mVnMid{{UNUSED_VN, UNUSED_VN, UNUSED_VN}};
  // Node between L_i and the ideal transformer. UNUSED_VN for the ref winding
  std::array<UInt, NumWindings> mVnPreIdeal{{UNUSED_VN, UNUSED_VN, UNUSED_VN}};
  // Slot holding the ideal transformer's branch-current unknown. UNUSED_VN
  // for the ref winding. Carries no voltage, it exists to widen
  // the system by one row and one column.
  std::array<UInt, NumWindings> mVnCurrent{{UNUSED_VN, UNUSED_VN, UNUSED_VN}};

  // When true, each star branch is a separate Resistor plus Inductor
  // Default to true
  Bool mWithResistiveLosses;

  // Per-unit quantities (system base)

  Real mNominalOmega = 0;
  Real mBaseApparentPower = 0;
  Real mBaseOmega = 0;
  Real mBaseImpedance = 0;
  Real mBaseAdmittance = 0;
  Real mBaseInductance = 0;
  Real mBaseCurrent = 0;

  // Star-branch leakage impedance of each winding [pu on the system base]
  // Filled by calculatePerUnitParameters().
  std::array<Complex, NumWindings> mLeakagePerUnit{
      {Complex(0., 0.), Complex(0., 0.), Complex(0., 0.)}};
  // Off-nominal tap of each winding [pu]. 1.0 when the winding sits at its
  // nominal ratio; complex when phase shifting
  std::array<Complex, NumWindings> mRatioPerUnit{
      {Complex(1., 0.), Complex(1., 0.), Complex(1., 0.)}};

  // The reduced 3x3 bus admittance block handed to the power-flow solver
  MatrixComp mY_element;

  // Helpers 
  // Computes how many virtual nodes this configuration needs
  static UInt virtualNodeCount(Bool withResistiveLosses);
  // Fills mVnStar / mVnMid / mVnPreIdeal / mVnCurrent
  // Must run after the ref winding is resolved
  void assignVirtualNodeSlots();
  // Wires one winding's sub-components between the star point and its
  // terminal, honoring the virtual node assignment
  void connectWinding(Winding w);
  // Stamps one ideal transformer's extra row and column.
  void stampIdealTransformer(SparseMatrixRow &systemMatrix, Winding w);

public:
  // Base voltage [V], set to the reference winding's nominal voltage
  const Attribute<Real>::Ptr mBaseVoltage;

  // Power flow results (one entry per winding)
  // Branch current at each terminal [A], indexed by Winding3W
  const Attribute<MatrixComp>::Ptr mCurrent;
  // Branch active power at each terminal [W]
  const Attribute<Matrix>::Ptr mActivePowerBranch;
  // Branch reactive power at each terminal [Var]
  const Attribute<Matrix>::Ptr mReactivePowerBranch;
  // Nodal active power injection at each terminal [W], indexed by Winding3W
  const Attribute<Matrix>::Ptr mActivePowerInjection;
  // Nodal reactive power injection at each terminal [Var], indexed by Winding3W
  const Attribute<Matrix>::Ptr mReactivePowerInjection;

  // Defines UID, name and logging level.
  Transformer3W(String uid, String name,
                Logger::Level logLevel = Logger::Level::off,
                Bool withResistiveLosses = true);
  // Defines name and logging level.
  Transformer3W(String name, Logger::Level logLevel = Logger::Level::off)
      : Transformer3W(name, name, logLevel) {}

  SimPowerComp<Complex>::Ptr clone(String name) override;

  // General
  // Constructs and registers MNA sub-components
  void createSubComponents() override;
  // Initializes component from power flow data
  void initializeParentFromNodesAndTerminals(Real frequency) override;

  // Powerflow
  // Nominal voltage of a winding [V]. Used by PFSolver to resolve the base
  // voltage of whichever node it is asking about
  Real getNominalVoltage(Winding w) const { return nominalVoltage(w); }
  // Sets the base voltage (ref winding voltage)
  void setBaseVoltage(Real baseVoltage);
  // Computes per-unit parameters on the given system base.
  void calculatePerUnitParameters(Real baseApparentPower, Real baseOmega);
  // Stamps the reduced 3x3 admittance block into the bus admittance matrix.
  void pfApplyAdmittanceMatrixStamp(SparseMatrixCompRow &Y) override;
  // Stores branch current and power flow; inputs are per-unit, stored in SI.
  void updateBranchFlow(VectorComp &current, VectorComp &powerflow);
  // Stores nodal injection power for one terminal. Called once per terminal
  // by PFSolverPowerPolar::calculateNodalInjection()
  void storeNodalInjection(Winding w, Complex powerInjection);
  // The reduced 3x3 admittance block. Rows/columns are indexed by Winding3W.
  MatrixComp Y_element();

  // MNA Section

  // Initializes internal variables of the component
  void mnaParentInitialize(Real omega, Real timeStep,
                           Attribute<Matrix>::Ptr leftVector) override;
  // Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  // Updates internal current variable of the component
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
  // Updates internal voltage variable of the component
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  // MNA pre step operations
  void mnaParentPreStep(Real time, Int timeStepCount) override;
  // MNA post step operations
  void mnaParentPostStep(Real time, Int timeStepCount,
                         Attribute<Matrix>::Ptr &leftVector) override;
  // Add MNA pre step dependencies
  void mnaParentAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;
  // Add MNA post step dependencies
  void
  mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                   AttributeBase::List &attributeDependencies,
                                   AttributeBase::List &modifiedAttributes,
                                   Attribute<Matrix>::Ptr &leftVector) override;
};
} // namespace Ph1
} // namespace SP
} // namespace CPS
