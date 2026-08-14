/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <algorithm>
#include <tuple>

#include <dpsim/PFSolverPowerPolar.h>

using namespace DPsim;
using namespace CPS;

PFSolverPowerPolar::PFSolverPowerPolar(CPS::String name,
                                       const CPS::SystemTopology &system,
                                       CPS::Real timeStep,
                                       CPS::Logger::Level logLevel)
    : PFSolver(name, system, timeStep, logLevel) {}

void PFSolverPowerPolar::generateInitialSolution(Real time,
                                                 bool keep_last_solution) {
  // Undo Q-limit PV<->PQ conversions left over from a previous solve.
  if (mEnforceReactiveLimits)
    resetToOriginalClassification();

  const UInt n = mSystem.mNodes.size();

  bool can_keep = keep_last_solution && mHasLastConvergedSolution &&
                  mLastConvergedV.size() == n && mLastConvergedD.size() == n;

  SPDLOG_LOGGER_INFO(mSLog,
                     "PF initialization: keep_last_solution={}, can_keep={}",
                     keep_last_solution, can_keep);

  resize_sol(n);
  resize_complex_sol(n);

  if (can_keep) {
    sol_V = mLastConvergedV;
    sol_D = mLastConvergedD;
  }

  if (!can_keep)
    mLmLambda = 1e-3;

  // update components
  for (auto comp : mSystem.mComponents) {
    if (auto load = std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
      if (load->use_profile)
        load->updatePQ(time);

      load->calculatePerUnitParameters(mBaseApparentPower,
                                       mSystem.mSystemOmega);
    }

    if (auto gen =
            std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(comp)) {
      gen->calculatePerUnitParameters(mBaseApparentPower, mSystem.mSystemOmega);
    }
  }

  // PQ buses
  for (auto pq : mPQBuses) {
    UInt idx = pq->matrixNodeIndex();

    sol_P(idx) = 0.0;
    sol_Q(idx) = 0.0;

    if (!can_keep) {
      if (auto v = pq->initialSingleVoltage(); std::abs(v) > 1e-6){
        Real base = mBaseVoltageAtNode[pq]; 
        sol_V(idx) = std::abs(v) / base; 
        sol_D(idx) = std::arg(v); 
      } else {
        sol_V(idx) = 1.0;
        sol_D(idx) = 0.0;
      }
    }

    for (auto comp : mSystem.mComponentsAtNode[pq]) {
      if (auto load = std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
        sol_P(idx) -= load->attributeTyped<CPS::Real>("P_pu")->get();
        sol_Q(idx) -= load->attributeTyped<CPS::Real>("Q_pu")->get();
      } else if (auto sst = std::dynamic_pointer_cast<
                     CPS::SP::Ph1::SolidStateTransformer>(comp)) {
        sol_P(idx) -= sst->getNodalInjection(pq).real();
        sol_Q(idx) -= sst->getNodalInjection(pq).imag();
      } else if (auto vsi = std::dynamic_pointer_cast<
                     CPS::SP::Ph1::AvVoltageSourceInverterDQ>(comp)) {
        sol_P(idx) +=
            vsi->attributeTyped<CPS::Real>("P_ref")->get() / mBaseApparentPower;
        sol_Q(idx) +=
            vsi->attributeTyped<CPS::Real>("Q_ref")->get() / mBaseApparentPower;
      } else if (auto gen =
                     std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(
                         comp)) {
        sol_P(idx) += gen->attributeTyped<CPS::Real>("P_set_pu")->get();
        sol_Q(idx) += gen->attributeTyped<CPS::Real>("Q_set_pu")->get();
      }
    }

    sol_S_complex(idx) = CPS::Complex(sol_P(idx), sol_Q(idx));
    sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
  }

  // PV buses
  for (auto pv : mPVBuses) {
    UInt idx = pv->matrixNodeIndex();

    sol_P(idx) = 0.0;
    sol_Q(idx) = 0.0;

    if (!can_keep) {
      if (auto v = pv->initialSingleVoltage(); std::abs(v) > 1e-6){
        Real base = mBaseVoltageAtNode[pv]; 
        sol_V(idx) = std::abs(v) / base; 
        sol_D(idx) = std::arg(v); 
      } else {
        sol_D(idx) = 0.0;
        sol_V(idx) = 1.0;
      }
    }

    for (auto comp : mSystem.mComponentsAtNode[pv]) {
      if (auto gen = std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(
              comp)) {
        sol_P(idx) += gen->attributeTyped<CPS::Real>("P_set_pu")->get();
        auto it = mLocalVSetOverride.find(idx);
        sol_V(idx) = (it != mLocalVSetOverride.end())
                      ? it->second
                      : gen->attributeTyped<CPS::Real>("V_set_pu")->get();
        // sol_P(idx) += gen->attributeTyped<CPS::Real>("P_set_pu")->get();
        // sol_V(idx) = gen->attributeTyped<CPS::Real>("V_set_pu")->get();
      } else if (auto load =
                     std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
        sol_P(idx) -= load->attributeTyped<CPS::Real>("P_pu")->get();
        sol_Q(idx) -= load->attributeTyped<CPS::Real>("Q_pu")->get();
      } else if (auto vsi = std::dynamic_pointer_cast<
                     CPS::SP::Ph1::AvVoltageSourceInverterDQ>(comp)) {
        sol_P(idx) +=
            vsi->attributeTyped<CPS::Real>("P_ref")->get() / mBaseApparentPower;
      } else if (auto extnet =
                     std::dynamic_pointer_cast<CPS::SP::Ph1::NetworkInjection>(
                         comp)) {
        sol_P(idx) += extnet->attributeTyped<CPS::Real>("p_inj")->get() /
                      mBaseApparentPower;
        sol_V(idx) = extnet->attributeTyped<CPS::Real>("V_set_pu")->get();
      }
    }

    sol_S_complex(idx) = CPS::Complex(sol_P(idx), sol_Q(idx));
    sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
  }

  // VD / slack buses
  for (auto vd : mVDBuses) {
    UInt idx = vd->matrixNodeIndex();

    sol_P(idx) = 0.0;
    sol_Q(idx) = 0.0;
    sol_D(idx) = 0.0;
    sol_V(idx) = 1.0;

    for (auto comp : mSystem.mComponentsAtNode[vd]) {
      if (auto extnet =
              std::dynamic_pointer_cast<CPS::SP::Ph1::NetworkInjection>(comp)) {
        sol_V(idx) = extnet->attributeTyped<CPS::Real>("V_set_pu")->get();
      } else if (auto load =
                     std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
        sol_P(idx) -= load->attributeTyped<CPS::Real>("P_pu")->get();
        sol_Q(idx) -= load->attributeTyped<CPS::Real>("Q_pu")->get();
      } else if (auto gen =
                     std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(
                         comp)) {
        sol_P(idx) += gen->attributeTyped<CPS::Real>("P_set_pu")->get();
        sol_Q(idx) += gen->attributeTyped<CPS::Real>("Q_set_pu")->get();
        sol_V(idx) = gen->attributeTyped<CPS::Real>("V_set_pu")->get();
      }
    }

    sol_S_complex(idx) = CPS::Complex(sol_P(idx), sol_Q(idx));
    sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
  }

  solutionInitialized = true;
  solutionComplexInitialized = true;

  Pesp = sol_P;
  Qesp = sol_Q;
}

void PFSolverPowerPolar::calculateMismatch() {
  UInt npqpv = mNumPQBuses + mNumPVBuses;
  UInt k;
  mF.setZero();

  for (UInt a = 0; a < npqpv; ++a) {
    // For PQ and PV buses calculate active power mismatch
    k = mPQPVBusIndices[a];
    mF(a) = Pesp.coeff(k) - P(k);

    //only for PQ buses calculate reactive power mismatch
    if (a < mNumPQBuses)
      mF(a + npqpv) = Qesp.coeff(k) - Q(k);
  }
}

void PFSolverPowerPolar::calculateJacobian() {
  UInt npqpv = mNumPQBuses + mNumPVBuses;
  double val;
  UInt k, j;
  UInt da, db;

  mJ.setZero();

  //J1
  for (UInt a = 0; a < npqpv; ++a) { //rows
    k = mPQPVBusIndices[a];
    //diagonal
    mJ.coeffRef(a, a) = -Q(k) - B(k, k) * sol_V.coeff(k) * sol_V.coeff(k);

    //non diagonal elements
    for (UInt b = 0; b < npqpv; ++b) {
      if (b != a) {
        j = mPQPVBusIndices[b];
        val = sol_V.coeff(k) * sol_V.coeff(j) *
              (G(k, j) * sin(sol_D.coeff(k) - sol_D.coeff(j)) -
               B(k, j) * cos(sol_D.coeff(k) - sol_D.coeff(j)));
        //if (val != 0.0)
        mJ.coeffRef(a, b) = val;
      }
    }
  }

  //J2
  da = 0;
  db = npqpv;
  for (UInt a = 0; a < npqpv; ++a) { //rows
    k = mPQPVBusIndices[a];
    //diagonal
    //std::cout << "J2D:" << (a + da) << "," << (a + db) << std::endl;
    if (a < mNumPQBuses)
      mJ.coeffRef(a + da, a + db) =
          P(k) + G(k, k) * sol_V.coeff(k) * sol_V.coeff(k);

    //non diagonal elements
    for (UInt b = 0; b < mNumPQBuses; ++b) {
      if (b != a) {
        j = mPQPVBusIndices[b];
        val = sol_V.coeff(k) * sol_V.coeff(j) *
              (G(k, j) * cos(sol_D.coeff(k) - sol_D.coeff(j)) +
               B(k, j) * sin(sol_D.coeff(k) - sol_D.coeff(j)));
        //if (val != 0.0)
        //std::cout << "J2ij:" << (a + da) << "," << (b + db) << std::endl;
        mJ.coeffRef(a + da, b + db) = val;
      }
    }
  }

  //J3
  da = npqpv;
  db = 0;
  for (UInt a = 0; a < mNumPQBuses; ++a) { //rows
    k = mPQPVBusIndices[a];
    //diagonal
    //std::cout << "J3:" << (a + da) << "," << (a + db) << std::endl;
    mJ.coeffRef(a + da, a + db) =
        P(k) - G(k, k) * sol_V.coeff(k) * sol_V.coeff(k);

    //non diagonal elements
    for (UInt b = 0; b < npqpv; ++b) {
      if (b != a) {
        j = mPQPVBusIndices[b];
        val = sol_V.coeff(k) * sol_V.coeff(j) *
              (G(k, j) * cos(sol_D.coeff(k) - sol_D.coeff(j)) +
               B(k, j) * sin(sol_D.coeff(k) - sol_D.coeff(j)));
        //if (val != 0.0)
        //std::cout << "J3:" << (a + da) << "," << (b + db) << std::endl;
        mJ.coeffRef(a + da, b + db) = -val;
      }
    }
  }

  //J4
  da = npqpv;
  db = npqpv;
  for (UInt a = 0; a < mNumPQBuses; ++a) { //rows
    k = mPQPVBusIndices[a];
    //diagonal
    //std::cout << "J4:" << (a + da) << "," << (a + db) << std::endl;
    mJ.coeffRef(a + da, a + db) =
        Q(k) - B(k, k) * sol_V.coeff(k) * sol_V.coeff(k);

    //non diagonal elements
    for (UInt b = 0; b < mNumPQBuses; ++b) {
      if (b != a) {
        j = mPQPVBusIndices[b];
        val = sol_V.coeff(k) * sol_V.coeff(j) *
              (G(k, j) * sin(sol_D.coeff(k) - sol_D.coeff(j)) -
               B(k, j) * cos(sol_D.coeff(k) - sol_D.coeff(j)));
        if (val != 0.0) {
          //std::cout << "J4:" << (a + da) << "," << (b + db) << std::endl;
          mJ.coeffRef(a + da, b + db) = val;
        }
      }
    }
  }
}

void PFSolverPowerPolar::updateSolution() {
  UInt npqpv = mNumPQBuses + mNumPVBuses;

  CPS::Vector sol_V_prev = sol_V;
  CPS::Vector sol_D_prev = sol_D;
  Real prevMismatchNorm = mF.norm();

  const double maxDVpu = 0.1;
  const double maxDThetaRad = 0.2;

  double baseScale = 1.0;
  for (UInt a = 0; a < npqpv; ++a) {
    double dTheta = std::abs(mX.coeff(a));
    if (dTheta > maxDThetaRad)
      baseScale = std::min(baseScale, maxDThetaRad / dTheta);
  }
  for (UInt b = 0; b < mNumPQBuses; ++b) {
    double dVrel = std::abs(mX.coeff(npqpv + b));
    if (dVrel > maxDVpu)
      baseScale = std::min(baseScale, maxDVpu / dVrel);
  }

  double scale = baseScale;
  const int maxBacktracks = 6;
  int attemptsUsed = 0;
  for (int attempt = 0; attempt <= maxBacktracks; ++attempt) {
    attemptsUsed = attempt;
    for (UInt a = 0; a < npqpv; ++a) {
      UInt k = mPQPVBusIndices[a];
      sol_D(k) = sol_D_prev(k) + scale * mX.coeff(a);
      if (a < mNumPQBuses)
        sol_V(k) = sol_V_prev(k) * (1.0 + scale * mX.coeff(a + npqpv));
    }
    for (auto node : mSystem.mNodes) {
      UInt idx = node->matrixNodeIndex();
      sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
    }

    calculateMismatch();

    if (mF.norm() <= 0.99 * prevMismatchNorm || attempt == maxBacktracks)
      break;
    scale *= 0.5;
  }

  if (attemptsUsed == 0) {
    mLmLambda = std::max(mLmLambda * 0.5, 1e-6);
  } else {
    mLmLambda = std::min(mLmLambda * 4.0, 1e4);
  }

  if (attemptsUsed >= maxBacktracks) {
    ++mStagnantIterations;
  } else {
    mStagnantIterations = 0;
  }

  //perturb if stuck (heuristic)
  if (mStagnantIterations >= 5) {
    mLmLambda = 1e-2;
    mStagnantIterations = 0;
    SPDLOG_LOGGER_WARN(mSLog, "LM stagnated at local minimum; perturbing state to escape");

    for (UInt a = 0; a < npqpv; ++a) {
      UInt k = mPQPVBusIndices[a];
      sol_D(k) += 0.01 * (2.0 * ((double)rand() / RAND_MAX) - 1.0);
      if (a < mNumPQBuses)
        sol_V(k) *= (1.0 + 0.01 * (2.0 * ((double)rand() / RAND_MAX) - 1.0));
    }
    for (auto node : mSystem.mNodes) {
      UInt idx = node->matrixNodeIndex();
      sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
    }
    calculateMismatch();
  }

}

void PFSolverPowerPolar::setSolution() {

  if (!isConverged) {
    SPDLOG_LOGGER_WARN(mSLog, "Not converged within {} iterations",
                       mIterations);

    SPDLOG_LOGGER_WARN(mSLog, "Writing last iterate to result state "
                              "(not stored as warm-start solution).");
  } else {
    calculatePAndQAtSlackBus();
    calculateQAtPVBuses();
    calculatePAndQInjectionPQBuses();

    mLastConvergedV = sol_V;
    mLastConvergedD = sol_D;
    mHasLastConvergedSolution = true;

    SPDLOG_LOGGER_INFO(mSLog, "Converged in {} iterations", mIterations);
  }

  SPDLOG_LOGGER_INFO(mSLog, "Solution written to result state:");

  SPDLOG_LOGGER_INFO(mSLog, "Name\tP\t\tQ\t\tV\t\tD");

  for (auto node : mSystem.mNodes) {
    UInt idx = node->matrixNodeIndex();

    SPDLOG_LOGGER_INFO(
        mSLog, "{}\t{}\t{}\t{}\t{}",
        std::dynamic_pointer_cast<CPS::SimNode<CPS::Complex>>(node)->name(),
        sol_P[idx], sol_Q[idx], sol_V[idx], sol_D[idx]);
  }

  mSLog->flush();

  for (UInt i = 0; i < mSystem.mNodes.size(); ++i) {
    sol_S_complex(i) = CPS::Complex(sol_P.coeff(i), sol_Q.coeff(i));
    sol_V_complex(i) = Math::polar(sol_V.coeff(i), sol_D.coeff(i));
  }

  for (auto node : mSystem.mNodes) {
    auto simNode = std::dynamic_pointer_cast<CPS::SimNode<CPS::Complex>>(node);

    UInt idx = node->matrixNodeIndex();

    simNode->setVoltage(sol_V_complex(idx) * mBaseVoltageAtNode[node]);

    simNode->setPower(sol_S_complex(idx) * mBaseApparentPower);
  }

  calculateBranchFlow();
  calculateNodalInjection();
}

void PFSolverPowerPolar::calculateBranchFlow() {
  for (auto line : mLines) {
    VectorComp v(2);
    v(0) = sol_V_complex.coeff(line->node(0)->matrixNodeIndex());
    v(1) = sol_V_complex.coeff(line->node(1)->matrixNodeIndex());
    /// I = Y * V
    VectorComp current = line->Y_element() * v;
    /// pf on branch [S_01; S_10] = [V_0 * conj(I_0); V_1 * conj(I_1)]
    VectorComp flow_on_branch =
        Math::elementwiseProduct(v, current.conjugate());
    line->updateBranchFlow(current, flow_on_branch);
  }
  for (auto trafo : mTransformers) {
    VectorComp v(2);
    v(0) = sol_V_complex.coeff(trafo->node(0)->matrixNodeIndex());
    v(1) = sol_V_complex.coeff(trafo->node(1)->matrixNodeIndex());
    /// I = Y * V
    VectorComp current = trafo->Y_element() * v;
    /// pf on branch [S_01; S_10] = [V_0 * conj(I_0); V_1 * conj(I_1)]
    VectorComp flow_on_branch =
        Math::elementwiseProduct(v, current.conjugate());
    trafo->updateBranchFlow(current, flow_on_branch);
  }
}

void PFSolverPowerPolar::calculateNodalInjection() {
  for (auto node : mSystem.mNodes) {
    std::list<std::shared_ptr<CPS::SP::Ph1::PiLine>> lines;
    for (auto comp : mSystem.mComponentsAtNode[node]) {
      if (std::shared_ptr<CPS::SP::Ph1::PiLine> line =
              std::dynamic_pointer_cast<CPS::SP::Ph1::PiLine>(comp)) {
        line->storeNodalInjection(sol_S_complex.coeff(node->matrixNodeIndex()));
        lines.push_back(line);
        break;
      }
    }
    if (lines.empty()) {
      for (auto comp : mSystem.mComponentsAtNode[node]) {
        if (std::shared_ptr<CPS::SP::Ph1::Transformer> trafo =
                std::dynamic_pointer_cast<CPS::SP::Ph1::Transformer>(comp)) {
          trafo->storeNodalInjection(
              sol_S_complex.coeff(node->matrixNodeIndex()));
          break;
        }
      }
    }
  }
}

Real PFSolverPowerPolar::P(UInt k) {
  Real val = 0.0;
  for (UInt j = 0; j < mSystem.mNodes.size(); ++j) {
    val += sol_V.coeff(j) * (G(k, j) * cos(sol_D.coeff(k) - sol_D.coeff(j)) +
                             B(k, j) * sin(sol_D.coeff(k) - sol_D.coeff(j)));
  }
  return sol_V.coeff(k) * val;
}

Real PFSolverPowerPolar::Q(UInt k) {
  Real val = 0.0;
  for (UInt j = 0; j < mSystem.mNodes.size(); ++j) {
    val += sol_V.coeff(j) * (G(k, j) * sin(sol_D.coeff(k) - sol_D.coeff(j)) -
                             B(k, j) * cos(sol_D.coeff(k) - sol_D.coeff(j)));
  }
  return sol_V.coeff(k) * val;
}

void PFSolverPowerPolar::calculatePAndQAtSlackBus() {
  for (auto topoNode : mVDBuses) {
    auto node_idx = topoNode->matrixNodeIndex();

    // Net nodal injection into the network: S_inj = V * conj(YV)
    CPS::Complex I(0.0, 0.0);
    for (UInt j = 0; j < mSystem.mNodes.size(); ++j)
      I += mY.coeff(node_idx, j) * sol_Vcx(j);

    CPS::Complex S = sol_Vcx(node_idx) * conj(I);

    // Generator/source power: S_gen = S_inj + S_load
    CPS::Complex Sgen = S;
    for (auto comp : mSystem.mComponentsAtNode[topoNode]) {
      if (auto loadPtr = std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
        Sgen += CPS::Complex(**(loadPtr->mActivePowerPerUnit),
                             **(loadPtr->mReactivePowerPerUnit));
      }
    }

    // Update connected VD source/generator with actual generated power
    for (auto comp : mSystem.mComponentsAtNode[topoNode]) {
      if (auto extnetPtr =
              std::dynamic_pointer_cast<CPS::SP::Ph1::NetworkInjection>(comp)) {
        extnetPtr->updatePowerInjection(Sgen * mBaseApparentPower);
        break;
      }

      if (auto sgPtr =
              std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(
                  comp)) {
        sgPtr->updatePowerInjection(Sgen * mBaseApparentPower);
        break;
      }
    }

    // Store net nodal injection, not generator power
    sol_P(node_idx) = S.real();
    sol_Q(node_idx) = S.imag();
  }
}

void PFSolverPowerPolar::calculateQAtPVBuses() {
  for (auto topoNode : mPVBuses) {
    auto node_idx = topoNode->matrixNodeIndex();

    // Net nodal injection into the network: S_inj = V * conj(YV)
    CPS::Complex I(0.0, 0.0);
    for (UInt j = 0; j < mSystem.mNodes.size(); ++j)
      I += mY.coeff(node_idx, j) * sol_Vcx(j);

    CPS::Complex S = sol_Vcx(node_idx) * conj(I);

    // Generator power: S_gen = S_inj + S_load
    CPS::Complex Sgen = S;
    for (auto comp : mSystem.mComponentsAtNode[topoNode]) {
      if (auto loadPtr = std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
        Sgen += CPS::Complex(**(loadPtr->mActivePowerPerUnit),
                             **(loadPtr->mReactivePowerPerUnit));
      }
    }

    // Update PV generator with actual generator Q
    for (auto comp : mSystem.mComponentsAtNode[topoNode]) {
      if (auto sgPtr =
              std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(
                  comp)) {
        sgPtr->updateReactivePowerInjection(Sgen * mBaseApparentPower);
        break;
      }
    }

    // Store net nodal Q injection, not generator Q
    sol_Q(node_idx) = S.imag();
  }
}

CPS::Real
PFSolverPowerPolar::loadReactivePowerPerUnit(CPS::TopologicalNode::Ptr node) {
  CPS::Real q = 0.0;
  for (auto comp : mSystem.mComponentsAtNode[node])
    if (auto load = std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp))
      q += load->attributeTyped<CPS::Real>("Q_pu")->get();
  return q;
}

CPS::Real PFSolverPowerPolar::generatorReactivePowerPerUnit(
    CPS::TopologicalNode::Ptr node) {
  UInt k = node->matrixNodeIndex();
  CPS::Complex I(0.0, 0.0);
  for (UInt j = 0; j < mSystem.mNodes.size(); ++j)
    I += mY.coeff(k, j) * sol_Vcx(j);
  // Net nodal injection S = generator - load; add load back for generator Q.
  CPS::Complex S = sol_Vcx(k) * conj(I);
  return S.imag() + loadReactivePowerPerUnit(node);
}

CPS::Bool PFSolverPowerPolar::enforceReactiveLimits() {
  static int callCount = 0;
  SPDLOG_LOGGER_INFO(mSLog, "=== enforceReactiveLimits() call #{} ===", ++callCount);

  // Returns false if the bus has no generator or no finite Q limit.
  auto busLimits = [&](CPS::TopologicalNode::Ptr node, CPS::Real &qMaxPU,
                       CPS::Real &qMinPU, CPS::Real &vSetPU) -> bool {
    qMaxPU = 0.0;
    qMinPU = 0.0;
    vSetPU = 0.0;
    bool hasGen = false, anyFinite = false;
    for (auto comp : mSystem.mComponentsAtNode[node]) {
      if (auto gen = std::dynamic_pointer_cast<CPS::SP::Ph1::SynchronGenerator>(
              comp)) {
        CPS::Real gMax = gen->attributeTyped<CPS::Real>("Q_max_pu")->get();
        CPS::Real gMin = gen->attributeTyped<CPS::Real>("Q_min_pu")->get();
        qMaxPU += gMax;
        qMinPU += gMin;
        vSetPU = gen->attributeTyped<CPS::Real>("V_set_pu")->get();
        hasGen = true;
        // isFinite is bit-pattern based, so it's immune to -Ofast/-ffast-math.
        if (CPS::Math::isFinite(gMax) || CPS::Math::isFinite(gMin))
          anyFinite = true;
      }
    }
    return hasGen && anyFinite;
  };

  auto frozen = [&](CPS::TopologicalNode::Ptr node) {
    auto it = mQLimitSwitchCount.find(node);
    return it != mQLimitSwitchCount.end() &&
           it->second >= mMaxQLimitSwitchesPerBus;
  };

  std::vector<std::tuple<CPS::TopologicalNode::Ptr, bool, CPS::Real>> toPQ;
  std::vector<CPS::TopologicalNode::Ptr> toPV;

  // PV -> PQ: voltage-controlling generators that hit a reactive limit.
  for (auto node : mPVBuses) {
    if (frozen(node))
      continue;
    CPS::Real qMaxPU, qMinPU, vSetPU;
    if (!busLimits(node, qMaxPU, qMinPU, vSetPU))
      continue;
    CPS::Real qGen = generatorReactivePowerPerUnit(node);
    if (CPS::Math::isFinite(qMaxPU) && qGen > qMaxPU)
      toPQ.emplace_back(node, true, qMaxPU);
    else if (CPS::Math::isFinite(qMinPU) && qGen < qMinPU)
      toPQ.emplace_back(node, false, qMinPU);
  }

  // PQ -> PV: pinned generators whose constraint is no longer binding.
  for (auto &kv : mQLimitConvertedAtMax) {
    auto node = kv.first;
    bool atMax = kv.second;
    if (frozen(node))
      continue;
    CPS::Real qMaxPU, qMinPU, vSetPU;
    if (!busLimits(node, qMaxPU, qMinPU, vSetPU))
      continue;
    CPS::Real v = sol_V(node->matrixNodeIndex());
    if (atMax && v > vSetPU)
      toPV.push_back(node);
    else if (!atMax && v < vSetPU)
      toPV.push_back(node);
  }

  SPDLOG_LOGGER_INFO(mSLog, "enforceReactiveLimits: toPQ.size()={} toPV.size()={}",
                     toPQ.size(), toPV.size());
  for (auto &c : toPQ)
    SPDLOG_LOGGER_INFO(mSLog, "  toPQ: bus {} pinning at {}",
                       std::get<0>(c)->name(), std::get<1>(c) ? "Qmax" : "Qmin");
  for (auto &n : toPV)
    SPDLOG_LOGGER_INFO(mSLog, "  toPV: bus {} relaxing back to PV", n->name());

  // Apply PV -> PQ switches (pin reactive injection at the limit via continuation).
  for (auto &c : toPQ) {
    auto node = std::get<0>(c);
    bool atMax = std::get<1>(c);
    CPS::Real qLimPU = std::get<2>(c);
    UInt idx = node->matrixNodeIndex();

    // Capture Q while still voltage-controlled, before switching bus type --
    // this is the continuation ramp's starting point.
    CPS::Real qStart = generatorReactivePowerPerUnit(node);

    mPVBuses.erase(std::remove(mPVBuses.begin(), mPVBuses.end(), node),
                  mPVBuses.end());
    mPQBuses.push_back(node);
    mQLimitConvertedAtMax[node] = atMax;
    if (++mQLimitSwitchCount[node] >= mMaxQLimitSwitchesPerBus)
      SPDLOG_LOGGER_WARN(
          mSLog, "Q-limit: bus {} frozen at {} after {} switches", node->name(),
          atMax ? "Qmax" : "Qmin", mQLimitSwitchCount[node]);
    else
      SPDLOG_LOGGER_INFO(mSLog, "Q-limit: PV bus {} -> PQ pinned at {}",
                        node->name(), atMax ? "Qmax" : "Qmin");

    // Rebuild index arrays now that this bus's classification changed, so the
    // continuation sub-steps below solve against the right unknown set.
    reclassifyBuses();

    CPS::Vector sol_V_prePin = sol_V;
    CPS::Vector sol_D_prePin = sol_D;

    // Pin Q directly at the limit.
    Qesp(idx) = qLimPU - loadReactivePowerPerUnit(node);
    sol_Q(idx) = Qesp(idx);

    SPDLOG_LOGGER_INFO(mSLog,
        "Q-pin: bus {} qStart={:.6f} -> qLim={:.6f}; attempting direct solve",
        node->name(), qStart, qLimPU);

     if (runNewtonRaphson(fmt::format("Q-pin direct bus {}", node->name()))) {
      SPDLOG_LOGGER_INFO(mSLog, "Q-pin: direct solve succeeded for bus {}",
                         node->name());
    } else {
      // The solution branch folds before Q reaches the limit, so continuing
      // from the pre-switch state cannot reach the Q-limited operating point.
      // Walk the load up instead, holding the pin fixed.
      SPDLOG_LOGGER_WARN(mSLog,
          "Q-pin: direct solve failed for bus {}; trying load homotopy",
          node->name());

      if (!solveWithLoadHomotopy(fmt::format("Q-pin bus {}", node->name()))) {
        SPDLOG_LOGGER_WARN(mSLog,
            "Q-pin: homotopy failed for bus {}; restoring pre-pin state",
            node->name());
        sol_V = sol_V_prePin;
        sol_D = sol_D_prePin;
        for (auto n : mSystem.mNodes) {
          UInt i = n->matrixNodeIndex();
          sol_V_complex(i) = Math::polar(sol_V(i), sol_D(i));
        }
      }
    }


    calculateMismatch();
  }

    // Apply PQ -> PV switches (restore voltage control).
  for (auto &node : toPV) {
    mPQBuses.erase(std::remove(mPQBuses.begin(), mPQBuses.end(), node),
                   mPQBuses.end());
    mPVBuses.push_back(node);
    CPS::Real qMaxPU, qMinPU, vSetPU;
    busLimits(node, qMaxPU, qMinPU, vSetPU);
    CPS::Real vBefore = sol_V(node->matrixNodeIndex());  
    sol_V(node->matrixNodeIndex()) = vSetPU;
    SPDLOG_LOGGER_WARN(mSLog,                           
        "Q-limit: PQ bus {} -> PV, voltage snapped {:.6f} -> {:.6f} ({:+.2f}%)",
        node->name(), vBefore, vSetPU, (vSetPU - vBefore) / vBefore * 100.0);
    mQLimitConvertedAtMax.erase(node);
    ++mQLimitSwitchCount[node];
    SPDLOG_LOGGER_INFO(mSLog, "Q-limit: PQ bus {} -> PV (constraint relaxed)",
                       node->name());
  }

  return !toPQ.empty() || !toPV.empty();
}

void PFSolverPowerPolar::computeLoadInjection(CPS::Vector &loadP,
                                              CPS::Vector &loadQ) {
  loadP = CPS::Vector::Zero(mSystem.mNodes.size());
  loadQ = CPS::Vector::Zero(mSystem.mNodes.size());
  for (auto node : mSystem.mNodes) {
    UInt idx = node->matrixNodeIndex();
    for (auto comp : mSystem.mComponentsAtNode[node]) {
      if (auto load = std::dynamic_pointer_cast<CPS::SP::Ph1::Load>(comp)) {
        loadP(idx) += load->attributeTyped<CPS::Real>("P_pu")->get();
        loadQ(idx) += load->attributeTyped<CPS::Real>("Q_pu")->get();
      }
    }
  }
}

Bool PFSolverPowerPolar::solveWithLoadHomotopy(const CPS::String &label) {
  CPS::Vector loadP, loadQ;
  computeLoadInjection(loadP, loadQ);

  // Full-load targets, restored before returning.
  CPS::Vector PespFull = Pesp;
  CPS::Vector QespFull = Qesp;

  const CPS::Real kLambdaStart = 0.5;
  const CPS::Real kMinStep     = 1e-3;
  const int       kMaxSteps    = 200;

  // Pesp = gen - load. At scale lambda: gen - lambda*load
  //                                   = Pesp_full + (1-lambda)*load
  auto applyLambda = [&](CPS::Real lambda) {
    for (auto node : mSystem.mNodes) {
      UInt i = node->matrixNodeIndex();
      if (std::abs(loadP(i)) < 1e-12 && std::abs(loadQ(i)) < 1e-12)
        continue;
      Pesp(i) = PespFull(i) + (1.0 - lambda) * loadP(i);
      Qesp(i) = QespFull(i) + (1.0 - lambda) * loadQ(i);
      sol_P(i) = Pesp(i);
      sol_Q(i) = Qesp(i);
    }
  };

  auto refreshComplex = [&]() {
    for (auto node : mSystem.mNodes) {
      UInt i = node->matrixNodeIndex();
      sol_V_complex(i) = Math::polar(sol_V(i), sol_D(i));
    }
  };

  auto restoreFullLoad = [&]() {
    Pesp = PespFull;
    Qesp = QespFull;
    sol_P = Pesp;
    sol_Q = Qesp;
  };

  // --- Phase 1: find a load level where the pinned problem is solvable.
  CPS::Real lambda = kLambdaStart;
  bool started = false;
  for (int t = 0; t < 8; ++t) {
    applyLambda(lambda);
    SPDLOG_LOGGER_INFO(mSLog, "Homotopy: seeking start at lambda={:.4f}", lambda);

    bool stepOk = runNewtonRaphson(fmt::format("{} homotopy start lambda={:.4f}",
                                               label, lambda));
    if (stepOk && !mRegulatedBusOfGen.empty())
      stepOk = resolveRemoteRegulation();

    if (stepOk) {
      started = true;
      break;
    }
    lambda *= 0.5;
  }
  if (!started) {
    SPDLOG_LOGGER_WARN(mSLog,
        "Homotopy: no solvable load level found down to lambda={:.4f}", lambda);
    restoreFullLoad();
    return false;
  }
  SPDLOG_LOGGER_INFO(mSLog, "Homotopy: started at lambda={:.4f}", lambda);

  // --- Phase 2: walk load up to 100%, warm-starting each step.
  CPS::Real step = (1.0 - lambda) / 4.0;
  int taken = 0;
  while (lambda < 1.0 - 1e-12 && taken < kMaxSteps) {
    CPS::Real trial = std::min(1.0, lambda + step);

    CPS::Vector Vsave = sol_V;
    CPS::Vector Dsave = sol_D;

    applyLambda(trial);
    SPDLOG_LOGGER_INFO(mSLog,
        "Homotopy: step {} lambda {:.4f} -> {:.4f} (dlambda={:.4f})",
        taken + 1, lambda, trial, step);
    
    bool stepOk = runNewtonRaphson(fmt::format("{} homotopy lambda={:.4f}", label, trial));
    
    if (stepOk && !mRegulatedBusOfGen.empty())
      stepOk = resolveRemoteRegulation();

    if (stepOk) {
      lambda = trial;
      ++taken;
      step *= 1.5;                      // grow cautiously after success

       // Are the remote-regulated buses being held where they should be?
      SPDLOG_LOGGER_INFO(mSLog,
          "Homotopy state @ lambda={:.4f}: SUB230-205={:.5f} (want 0.98000)  "
          "HYDRO-201={:.5f} (want 1.04000)  URBGEN-206={:.5f}  HYDRO_G-211={:.5f}  "
          "CATDOG-3008={:.5f}  CATDOG_G-3018={:.5f}",
          lambda, sol_V(10), sol_V(19), sol_V(18), sol_V(23),
          sol_V(6), sol_V(1));

    } else {
      sol_V = Vsave;                    // roll back to the last good point
      sol_D = Dsave;
      refreshComplex();
      step *= 0.5;
      SPDLOG_LOGGER_WARN(mSLog,
          "Homotopy: step failed, shrinking to {:.6f}", step);
      if (step < kMinStep) {
        SPDLOG_LOGGER_WARN(mSLog,
            "Homotopy: stalled at lambda={:.4f} (step below {:.1e}); "
            "load level appears to be a limit point",
            lambda, kMinStep);
        // restoreFullLoad();
        // return false;
        break; 
      }
    }
  }
  // --- Phase 2b: nose reached below full load. Try to cross to the lower
  // branch by overshooting, then re-solving at full load from the far side.
  if (lambda < 1.0 - 1e-12) {
    SPDLOG_LOGGER_INFO(mSLog,
        "Homotopy: nose at lambda={:.4f}; attempting lower-branch crossing", lambda);

    CPS::Vector Vnose = sol_V;
    CPS::Vector Dnose = sol_D;

    for (CPS::Real over : {1.05, 1.10, 1.20}) {
      sol_V = Vnose;
      sol_D = Dnose;
      applyLambda(over);
      refreshComplex();
      runNewtonRaphson(fmt::format("{} overshoot lambda={:.2f}", label, over));
      // Deliberately ignore the result: we want the diverged iterate, which
      // usually sits past the fold on the lower branch.

      restoreFullLoad();
      refreshComplex();
      if (runNewtonRaphson(fmt::format("{} lower-branch from overshoot {:.2f}",
                                       label, over))) {
        SPDLOG_LOGGER_INFO(mSLog,
            "Homotopy: reached a full-load solution via overshoot {:.2f}", over);
        return true;
      }
    }

    sol_V = Vnose;
    sol_D = Dnose;
    refreshComplex();
    restoreFullLoad();
    return false;
  }

  // --- Phase 3: land exactly on full load.
  restoreFullLoad();
  bool finalOk = runNewtonRaphson(fmt::format("{} homotopy final", label));
  if (finalOk && !mRegulatedBusOfGen.empty())
    finalOk = resolveRemoteRegulation();
  return finalOk;
}

void PFSolverPowerPolar::clearReactiveLimitState() {
  mQLimitConvertedAtMax.clear();
  mQLimitSwitchCount.clear();
}

void PFSolverPowerPolar::calculatePAndQInjectionPQBuses() {
  // calculates apparent power injection at PQ buses flowing to other nodes (i.e. S_inj_to_other = S_inj - S_shunt, with S_inj = S_gen - S_load)
  for (auto topoNode : mPQBuses) {
    auto node_idx = topoNode->matrixNodeIndex();

    // calculate power flowing out of the node into the admittance matrix (i.e. S_inj)
    CPS::Complex I(0.0, 0.0);
    for (UInt j = 0; j < mSystem.mNodes.size(); ++j)
      I += mY.coeff(node_idx, j) * sol_Vcx(j);
    CPS::Complex S = sol_Vcx(node_idx) * conj(I);

    // Subtracting shunt power to obtain power injection flowing from this node to the other nodes (i.e. S_inj_to_other)
    CPS::Real V = sol_V.coeff(node_idx);
    for (auto comp : mSystem.mComponentsAtNode[topoNode])
      if (auto shuntPtr = std::dynamic_pointer_cast<CPS::SP::Ph1::Shunt>(comp))
        // capacitive susceptance is positive --> q is injected into the node
        S += std::pow(V, 2) * Complex(-**(shuntPtr->mConductancePerUnit),
                                      **(shuntPtr->mSusceptancePerUnit));

    // TODO: check whether S_inj_to_other should be stored in sol_P and sol_Q or rather S_inj
    sol_P(node_idx) = S.real();
    sol_Q(node_idx) = S.imag();
  }
}

void PFSolverPowerPolar::resize_sol(Int n) {
  sol_P = CPS::Vector(n);
  sol_Q = CPS::Vector(n);
  sol_V = CPS::Vector(n);
  sol_D = CPS::Vector(n);
  sol_P.setZero(n);
  sol_Q.setZero(n);
  sol_V.setZero(n);
  sol_D.setZero(n);
}

void PFSolverPowerPolar::resize_complex_sol(Int n) {
  sol_S_complex = CPS::VectorComp(n);
  sol_V_complex = CPS::VectorComp(n);
  sol_S_complex.setZero(n);
  sol_V_complex.setZero(n);
}

CPS::Real PFSolverPowerPolar::sol_Vr(UInt k) {
  return sol_V.coeff(k) * cos(sol_D.coeff(k));
}

CPS::Real PFSolverPowerPolar::sol_Vi(UInt k) {
  return sol_V.coeff(k) * sin(sol_D.coeff(k));
}

CPS::Complex PFSolverPowerPolar::sol_Vcx(UInt k) {
  return CPS::Complex(sol_Vr(k), sol_Vi(k));
}
