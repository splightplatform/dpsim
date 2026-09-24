#pragma once
#include <array>

#include <dpsim-models/Base/Base_Ph3_Transformer3W.h>
#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/EMT/EMT_Ph3_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/Solver/MNAInterface.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

// EMT Three-winding transformer
class Transformer3W : public CompositePowerComp<Real>,
                      public Base::Ph3::Transformer3W,
                      public SharedFactory<Transformer3W> {

public:
  using Winding = Base::Ph3::Winding3W;

private:
  static constexpr UInt UNUSED_VN = static_cast<UInt>(-1);

  // --- Wye-branch sub-components 
  std::array<std::shared_ptr<EMT::Ph3::Resistor>, NumWindings> mSubResistor;
  std::array<std::shared_ptr<EMT::Ph3::Inductor>, NumWindings> mSubInductor;

  // --- Delta-leg sub-components -------------------------------------------

  // --- Snubbers
  std::array<std::shared_ptr<EMT::Ph3::Resistor>, NumWindings>
      mSubSnubResistor;
  std::array<std::shared_ptr<EMT::Ph3::Capacitor>, NumWindings>
      mSubSnubCapacitor;
  std::array<Matrix, NumWindings> mSnubberResistance{
      {Matrix::Zero(3, 3), Matrix::Zero(3, 3), Matrix::Zero(3, 3)}};
  std::array<Matrix, NumWindings> mSnubberCapacitance{
      {Matrix::Zero(3, 3), Matrix::Zero(3, 3), Matrix::Zero(3, 3)}};

  //virtual nodes
  UInt mVnStar = 0;
  std::array<UInt, NumWindings> mVnMid{{UNUSED_VN, UNUSED_VN, UNUSED_VN}};
  std::array<UInt, NumWindings> mVnPreIdeal{
      {UNUSED_VN, UNUSED_VN, UNUSED_VN}};
  std::array<UInt, NumWindings> mVnCurrent{{UNUSED_VN, UNUSED_VN, UNUSED_VN}};

  Bool mWithResistiveLosses;
  static UInt virtualNodeCount(Bool withResistiveLosses);

  void assignVirtualNodeSlots();
  void connectWinding(Winding w);
  void stampIdealTransformer(SparseMatrixRow &systemMatrix, Winding w);

public:
  Transformer3W(String uid, String name,
                Logger::Level logLevel = Logger::Level::off,
                Bool withResistiveLosses = true);
  Transformer3W(String name, Logger::Level logLevel = Logger::Level::off)
      : Transformer3W(name, name, logLevel) {}

  SimPowerComp<Real>::Ptr clone(String name) override;

  // General
  void createSubComponents() override;
  void initializeParentFromNodesAndTerminals(Real frequency) override;

  // MNA Section
  void mnaParentInitialize(Real omega, Real timeStep,
                           Attribute<Matrix>::Ptr leftVector) override;
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  void mnaParentPreStep(Real time, Int timeStepCount) override;
  void mnaParentPostStep(Real time, Int timeStepCount,
                         Attribute<Matrix>::Ptr &leftVector) override;
  void mnaParentAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;
  void
  mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                   AttributeBase::List &attributeDependencies,
                                   AttributeBase::List &modifiedAttributes,
                                   Attribute<Matrix>::Ptr &leftVector) override;
};
} // namespace Ph3
} // namespace EMT
} // namespace CPS
