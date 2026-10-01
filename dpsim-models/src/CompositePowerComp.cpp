// SPDX-License-Identifier: Apache-2.0

#include <dpsim-models/CompositePowerComp.h>

using namespace CPS;

template <typename VarType>
void CompositePowerComp<VarType>::initializeFromNodesAndTerminals(
    Real frequency) {
  this->createSubComponents();

  initializeParentFromNodesAndTerminals(frequency);

  for (auto subComp : this->mSubComponents) {
    subComp->initialize(this->mFrequencies);
    subComp->initializeFromNodesAndTerminals(frequency);
  }
}

template <typename VarType>
void CompositePowerComp<VarType>::addMNASubComponent(
    typename SimPowerComp<VarType>::Ptr subc,
    MNA_SUBCOMP_TASK_ORDER preStepOrder, MNA_SUBCOMP_TASK_ORDER postStepOrder,
    Bool contributeToRightVector) {
  this->mSubComponents.push_back(subc);
  if (auto mnasubcomp =
          std::dynamic_pointer_cast<MNASimPowerComp<VarType>>(subc)) {
    this->mSubcomponentsMNA.push_back(mnasubcomp);

    if (contributeToRightVector) {
      this->mRightVectorSubcomps.push_back(mnasubcomp); 
    }

    switch (preStepOrder) {
    case MNA_SUBCOMP_TASK_ORDER::NO_TASK:
      break;
    case MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT: {
      this->mSubcomponentsPreStepBeforeParent.push_back(mnasubcomp);
      break;
    }
    case MNA_SUBCOMP_TASK_ORDER::TASK_AFTER_PARENT: {
      this->mSubcomponentsPreStepAfterParent.push_back(mnasubcomp);
      break;
    }
    }
    switch (postStepOrder) {
    case MNA_SUBCOMP_TASK_ORDER::NO_TASK:
      break;
    case MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT: {
      this->mSubcomponentsPostStepBeforeParent.push_back(mnasubcomp);
      break;
    }
    case MNA_SUBCOMP_TASK_ORDER::TASK_AFTER_PARENT: {
      this->mSubcomponentsPostStepAfterParent.push_back(mnasubcomp);
      break;
    }
    }
  }
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  SimPowerComp<VarType>::updateMatrixNodeIndices();

  for (auto subComp : mSubcomponentsMNA) {
    subComp->mnaInitialize(omega, timeStep, leftVector);
  }

  **this->mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
  // for every MNA subcomponent (w/ R vector contribution)
  // add the subcomponents R vector rows to parent's list
  for (auto &sub : mRightVectorSubcomps)
    this->mRightVectorRows.insert(this->mRightVectorRows.end(), 
                                  sub->mRightVectorRows.begin(), 
                                  sub->mRightVectorRows.end()); 

  mnaParentInitialize(omega, timeStep, leftVector);
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  for (auto subComp : mSubcomponentsMNA) {
    subComp->mnaApplySystemMatrixStamp(systemMatrix);
  }
  mnaParentApplySystemMatrixStamp(systemMatrix);
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  // clear the composite's sums from last step only at rows that matter
  for (UInt r : this->mRightVectorRows)
    rightVector(r, 0) = 0; 
  
  // accumlate the contribution from each subcomponent's RightVectorRows
  for (auto &sub : mRightVectorSubcomps){ 
    const Matrix &stamp = **sub->mRightVector; 
    for (UInt r : sub->mRightVectorRows)
      rightVector(r, 0) += stamp(r, 0); 
  }
  mnaParentApplyRightSideVectorStamp(rightVector);
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompPreStep(Real time, Int timeStepCount) {
  for (auto subComp : mSubcomponentsPreStepBeforeParent) {
    subComp->mnaPreStep(time, timeStepCount);
  }
  mnaParentPreStep(time, timeStepCount);
  for (auto subComp : mSubcomponentsPreStepAfterParent) {
    subComp->mnaPreStep(time, timeStepCount);
  }
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  for (auto subComp : mSubcomponentsPostStepBeforeParent) {
    subComp->mnaPostStep(time, timeStepCount, leftVector);
  }
  mnaParentPostStep(time, timeStepCount, leftVector);
  for (auto subComp : mSubcomponentsPostStepAfterParent) {
    subComp->mnaPostStep(time, timeStepCount, leftVector);
  }
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  for (auto subComp : mSubcomponentsMNA) {
    subComp->mnaAddPreStepDependencies(
        prevStepDependencies, attributeDependencies, modifiedAttributes);
  }
  mnaParentAddPreStepDependencies(prevStepDependencies, attributeDependencies,
                                  modifiedAttributes);
}

template <typename VarType>
void CompositePowerComp<VarType>::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  for (auto subComp : mSubcomponentsMNA) {
    subComp->mnaAddPostStepDependencies(prevStepDependencies,
                                        attributeDependencies,
                                        modifiedAttributes, leftVector);
  }
  mnaParentAddPostStepDependencies(prevStepDependencies, attributeDependencies,
                                   modifiedAttributes, leftVector);
}

// Declare specializations to move definitions to .cpp
template class CPS::CompositePowerComp<Real>;
template class CPS::CompositePowerComp<Complex>;
