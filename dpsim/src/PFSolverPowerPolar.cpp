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

#if defined(WITH_KLU)
#include <dpsim/KLUAdapter.h>
#elif defined(WITH_SPARSE)
#include <dpsim/SparseLUAdapter.h>
#else
#include <dpsim/DenseLUAdapter.h>
#endif

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

Bool PFSolverPowerPolar::updateSolution() {
  UInt npqpv = mNumPQBuses + mNumPVBuses;

  CPS::Vector sol_V_prev = sol_V;
  CPS::Vector sol_D_prev = sol_D;
  CPS::Vector mF_prev = mF; 
  Real prevMismatchNorm = mF.norm();
  Real bestNorm = prevMismatchNorm; 
  CPS::Vector step = mX; // vector of individual steps for buses
  
  const double maxDVpu = 0.1;
  const double maxDThetaRad = 0.2;
  UInt nClamped = 0; 
  
  // Use clampled scales rather than one base scale
  // Δθ_k ← clamp(Δθ_k, −0.2, +0.2)        for every PQ and PV bus
  for (UInt a = 0; a < npqpv; ++a) {
    step(a) = std::clamp(mX.coeff(a), -maxDThetaRad, maxDThetaRad); 
    if (step(a) != mX.coeff(a)) ++nClamped;
  }
  
  // ΔV_k/V_k ← clamp(ΔV_k/V_k, −0.1, +0.1) for every PQ bus
  for (UInt b = 0; b < mNumPQBuses; ++b) {
    step(npqpv + b) = std::clamp(mX.coeff(npqpv + b), -maxDVpu, maxDVpu); 
    if (step(npqpv + b) != mX.coeff(npqpv + b)) ++nClamped; 
  }
  
  double scale = 1.0; 
  const int maxBacktracks = 6;
  Bool improved = false;
  for (int attempt = 0; attempt <= maxBacktracks; ++attempt) {
    for (UInt a = 0; a < npqpv; ++a) {
      UInt k = mPQPVBusIndices[a];
      sol_D(k) = sol_D_prev(k) + scale * step(a);
      if (a < mNumPQBuses)
        sol_V(k) = sol_V_prev(k) * (1.0 + scale * step(a + npqpv));
      
    }
    for (auto node : mSystem.mNodes) {
      UInt idx = node->matrixNodeIndex();
      sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
    }

    calculateMismatch();
    bestNorm = std::min(bestNorm, mF.norm()); 

    if (mF.norm() <= 0.99 * prevMismatchNorm){
      improved = true; 
      break;
    }
    scale *= 0.5;
  }
  if (!improved){ 
    // restore old values 
    sol_V = sol_V_prev; 
    sol_D = sol_D_prev; 
    mF = mF_prev; 

    // rebuild sol_V_complex
    for (auto node : mSystem.mNodes) {
      UInt idx = node->matrixNodeIndex();
      sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
    }

    // warn 
    SPDLOG_LOGGER_WARN(mSLog, "line search: no step improved |F| (prev={}, best={}, minScale={}); restoring previous state",
                   prevMismatchNorm, bestNorm, 2 * scale);
    return false; 
  }

  SPDLOG_LOGGER_INFO(mSLog, "step: {} entries clamped, accepted scale {}", nClamped, scale);
  return true;
}

void PFSolverPowerPolar::buildFastDecoupledMatrices() {
  const UInt npqpv = mNumPQBuses + mNumPVBuses;
  std::vector<Eigen::Triplet<Real>> bpTriplets, bppTriplets;
  bpTriplets.reserve(mY.nonZeros());
  bppTriplets.reserve(mY.nonZeros());

  // PQ buses come first in mPQPVBusIndices, so a PQ bus has the same row in B' and B''.
  for (UInt a = 0; a < npqpv; ++a) {
    const UInt k = mPQPVBusIndices[a];
    const bool kIsPQ = a < mNumPQBuses;
    Real bpDiag = 0.0;

    for (CPS::SparseMatrixCompRow::InnerIterator it(mY, k); it; ++it) {
      const UInt j = static_cast<UInt>(it.col());
      const Real negB = -it.value().imag();
      if (j == k) {
        if (kIsPQ)
          bppTriplets.emplace_back(a, a, negB);
        continue;
      }
      // Off-diagonals to the slack still count towards the shunt-free diagonal
      bpDiag -= negB;
      const Int b = mPQPVPositionOfBus[j];
      if (b < 0)
        continue;
      bpTriplets.emplace_back(a, b, negB);
      if (kIsPQ && static_cast<UInt>(b) < mNumPQBuses)
        bppTriplets.emplace_back(a, b, negB);
    }
    bpTriplets.emplace_back(a, a, bpDiag);
  }

  mBp.resize(npqpv, npqpv);
  mBp.setFromTriplets(bpTriplets.begin(), bpTriplets.end());
  mBp.makeCompressed();
  mBpp.resize(mNumPQBuses, mNumPQBuses);
  mBpp.setFromTriplets(bppTriplets.begin(), bppTriplets.end());
  mBpp.makeCompressed();

  std::vector<std::pair<UInt, UInt>> noVaryingEntries;
  auto factorizeOnce = [&](CPS::SparseMatrixRow &matrix,
                           std::shared_ptr<DirectLinearSolver> &solver) {
    if (matrix.rows() == 0) {
      solver.reset();
      return;
    }
#if defined(WITH_KLU)
    solver = std::make_shared<KLUAdapter>(mSLog);
#elif defined(WITH_SPARSE)
    solver = std::make_shared<SparseLUAdapter>(mSLog);
#else
    solver = std::make_shared<DenseLUAdapter>(mSLog);
#endif
    solver->preprocessing(matrix, noVaryingEntries);
    solver->factorize(matrix);
  };
  factorizeOnce(mBp, mBpSolver);
  factorizeOnce(mBpp, mBppSolver);

  mFastDecoupledBusIndices = mPQPVBusIndices;
  mFastDecoupledNumPQBuses = mNumPQBuses;
}

CPS::Bool PFSolverPowerPolar::runFastDecoupled(CPS::UInt maxIter,
                                               CPS::Real handoffTol) {
  const CPS::String label = "fast-decoupled";
  const UInt npqpv = mNumPQBuses + mNumPVBuses;
  if (npqpv == 0)
    return true;

  if (!mBpSolver || mFastDecoupledBusIndices != mPQPVBusIndices ||
      mFastDecoupledNumPQBuses != mNumPQBuses)
    buildFastDecoupledMatrices();

  SPDLOG_LOGGER_INFO(mSLog, ">>> fast-decoupled warm start (max {} iters, handoff at {} pu)",
                     maxIter, handoffTol);

  auto refreshComplexVoltages = [&]() {
    for (auto node : mSystem.mNodes) {
      UInt idx = node->matrixNodeIndex();
      sol_V_complex(idx) = Math::polar(sol_V(idx), sol_D(idx));
    }
  };

  calculateMismatch();
  reportProgress(label, 0);

  CPS::Vector bestV = sol_V;
  CPS::Vector bestD = sol_D;
  Real bestNorm = mF.norm();
  Real prevNorm = bestNorm;
  UInt growingIterations = 0;
  UInt iter = 0;
  Bool reachedHandoff = false;

  Matrix rhsP(npqpv, 1);
  Matrix rhsQ(mNumPQBuses, 1);

  while (true) {
    if (mF.cwiseAbs().maxCoeff() < handoffTol) {
      reachedHandoff = true;
      break;
    }
    if (iter >= maxIter)
      break;
    ++iter;

    // Angle half-step: B' dTheta = dP / V for PQ and PV buses
    for (UInt a = 0; a < npqpv; ++a)
      rhsP(a, 0) = mF(a) / sol_V(mPQPVBusIndices[a]);
    Matrix dTheta = mBpSolver->solve(rhsP);
    for (UInt a = 0; a < npqpv; ++a)
      sol_D(mPQPVBusIndices[a]) += dTheta(a, 0);

    calculateMismatch();

    // Magnitude half-step: B'' dV = dQ / V for PQ buses
    if (mNumPQBuses > 0) {
      for (UInt b = 0; b < mNumPQBuses; ++b)
        rhsQ(b, 0) = mF(npqpv + b) / sol_V(mPQPVBusIndices[b]);
      Matrix dV = mBppSolver->solve(rhsQ);
      for (UInt b = 0; b < mNumPQBuses; ++b)
        sol_V(mPQPVBusIndices[b]) += dV(b, 0);
    }

    calculateMismatch();
    reportProgress(label, iter);

    const Real norm = mF.norm();
    SPDLOG_LOGGER_INFO(mSLog, "fast-decoupled iter {} max|mismatch|={} norm={}",
                       iter, mF.cwiseAbs().maxCoeff(), norm);

    if (!Math::isFinite(norm)) {
      SPDLOG_LOGGER_WARN(mSLog, "fast-decoupled: mismatch not finite at iter {}; "
                                "handing best state to Newton", iter);
      break;
    }
    if (norm < bestNorm) {
      bestNorm = norm;
      bestV = sol_V;
      bestD = sol_D;
    }
    growingIterations = norm > prevNorm ? growingIterations + 1 : 0;
    prevNorm = norm;
    if (growingIterations >= 2) {
      SPDLOG_LOGGER_WARN(mSLog, "fast-decoupled: mismatch grew for 2 iterations; "
                                "handing best state to Newton");
      break;
    }
  }

  // Hand over the best state seen, which is not necessarily the last one
  if (!(mF.norm() <= bestNorm)) {
    sol_V = bestV;
    sol_D = bestD;
    calculateMismatch();
  }
  refreshComplexVoltages();

  SPDLOG_LOGGER_INFO(mSLog, "<<< fast-decoupled end: iters={} reached_handoff={} "
                            "max|mismatch|={} norm={}",
                     iter, reachedHandoff, mF.cwiseAbs().maxCoeff(), mF.norm());
  return reachedHandoff;
}

void PFSolverPowerPolar::setSolution() {

  if (!mSolveSucceeded) {
    if (isConverged)
      SPDLOG_LOGGER_WARN(mSLog, "Last Newton-Raphson run converged, but the "
                                "Q-limit or remote-regulation loop did not succeed");
    else
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
  for (auto xfmr : mTransformers3W) {
    VectorComp v(3); 
    v(0) = sol_V_complex.coeff(xfmr->node(0)->matrixNodeIndex()); 
    v(1) = sol_V_complex.coeff(xfmr->node(1)->matrixNodeIndex()); 
    v(2) = sol_V_complex.coeff(xfmr->node(2)->matrixNodeIndex()); 
    
    VectorComp current = xfmr->Y_element() * v; 
    VectorComp flow_on_branch = Math::elementwiseProduct(v, current.conjugate()); 
    xfmr->updateBranchFlow(current, flow_on_branch); 
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
        if (std::shared_ptr<CPS::SP::Ph1::Transformer3W> xfmr = 
                std::dynamic_pointer_cast<CPS::SP::Ph1::Transformer3W>(comp)) {
          if (xfmr->terminal(0)->node()->name() == node->name()){
            xfmr->storeNodalInjection(Base::Ph1::Winding3W::Primary, 
                                      sol_S_complex.coeff(node->matrixNodeIndex()));
            break; 
          }
          if (xfmr->terminal(1)->node()->name() == node->name()){
            xfmr->storeNodalInjection(Base::Ph1::Winding3W::Secondary, 
                                      sol_S_complex.coeff(node->matrixNodeIndex()));
            break; 
          }
          if (xfmr->terminal(2)->node()->name() == node->name()){
            xfmr->storeNodalInjection(Base::Ph1::Winding3W::Tertiary, 
                                      sol_S_complex.coeff(node->matrixNodeIndex()));
            break; 
          }
        }
      }
    }
  }
}

Real PFSolverPowerPolar::P(UInt k) {
  // P(k) = Vk · Σj Vj · (G·cos d + B·sin d)
  // d = θk − θj
  // use sparse matrix inner-iterator 
  // (it(mY, k) = iterator to first entry in row k
  // it.col() is j, it.value() is admittance, it.row() is k 
  Real val = 0.0;
  for (CPS::SparseMatrixCompRow::InnerIterator it(mY, k); it; ++it){
    const Real d = sol_D.coeff(k) - sol_D.coeff(it.col()); 
    val += sol_V.coeff(it.col()) * (it.value().real() * cos(d) + it.value().imag() * sin(d));
  }
  return sol_V.coeff(k) * val; 
}

Real PFSolverPowerPolar::Q(UInt k) {
  // Q(k) = Vk · Σj Vj · (G·sin d − B·cos d)
  Real val = 0.0;
  for (CPS::SparseMatrixCompRow::InnerIterator it(mY, k); it; ++it){
    const Real d = sol_D.coeff(k) - sol_D.coeff(it.col()); 
    val += sol_V.coeff(it.col()) * (it.value().real() * sin(d) - it.value().imag() * cos(d)); 
  }
  return sol_V.coeff(k) * val;
}

void PFSolverPowerPolar::calculatePAndQAtSlackBus() {
  for (auto topoNode : mVDBuses) {
    auto node_idx = topoNode->matrixNodeIndex();

    // Net nodal injection into the network: S_inj = V * conj(YV)
    CPS::Complex I(0.0, 0.0);
    // use inner iterator
    for (CPS::SparseMatrixCompRow::InnerIterator it(mY, node_idx); it; ++it)
      I += it.value() * sol_Vcx(it.col()); 

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
    // use inner iterator
    for (CPS::SparseMatrixCompRow::InnerIterator it(mY, node_idx); it; ++it)
      I += it.value() * sol_Vcx(it.col());

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
  // user inner iterator
  for (CPS::SparseMatrixCompRow::InnerIterator it(mY, k); it; ++it)
    I += it.value() * sol_Vcx(it.col());
  // Net nodal injection S = generator - load; add load back for generator Q.
  CPS::Complex S = sol_Vcx(k) * conj(I);
  return S.imag() + loadReactivePowerPerUnit(node);
}

CPS::Bool PFSolverPowerPolar::enforceReactiveLimits() {
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

  // node, pinned at Qmax?, Q limit [pu], violation beyond the limit [pu]
  std::vector<std::tuple<CPS::TopologicalNode::Ptr, bool, CPS::Real, CPS::Real>> toPQ;
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
      toPQ.emplace_back(node, true, qMaxPU, qGen - qMaxPU);
    else if (CPS::Math::isFinite(qMinPU) && qGen < qMinPU)
      toPQ.emplace_back(node, false, qMinPU, qMinPU - qGen);
  }

  // optionally switch only the worst violators this pass
  if (mMaxQLimitSwitchesPerPass > 0 && toPQ.size() > mMaxQLimitSwitchesPerPass) {
    std::sort(toPQ.begin(), toPQ.end(), [](const auto &l, const auto &r) {
      return std::get<3>(l) > std::get<3>(r);
    });
    toPQ.resize(mMaxQLimitSwitchesPerPass);
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

  // apply every switch of this pass at once

  // Apply PV -> PQ switches (pin reactive injection at the limit).
  for (auto &c : toPQ) {
    auto node = std::get<0>(c);
    bool atMax = std::get<1>(c);
    CPS::Real qLimPU = std::get<2>(c);
    UInt idx = node->matrixNodeIndex();

    mPVBuses.erase(std::remove(mPVBuses.begin(), mPVBuses.end(), node),
                   mPVBuses.end());
    mPQBuses.push_back(node);
    mQLimitConvertedAtMax[node] = atMax;
    if (++mQLimitSwitchCount[node] >= mMaxQLimitSwitchesPerBus)
      SPDLOG_LOGGER_WARN(
          mSLog, "Q-limit: bus {} frozen at {} after {} switches", node->name(),
          atMax ? "Qmax" : "Qmin", mQLimitSwitchCount[node]);
    else
      SPDLOG_LOGGER_INFO(mSLog, "Q-limit: PV bus {} -> PQ pinned at {} (violation {:.4f} pu)",
                         node->name(), atMax ? "Qmax" : "Qmin", std::get<3>(c));

    // Pin Q directly at the limit.
    Qesp(idx) = qLimPU - loadReactivePowerPerUnit(node);
    sol_Q(idx) = Qesp(idx);
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
    SPDLOG_LOGGER_INFO(mSLog,
        "Q-limit: PQ bus {} -> PV, voltage snapped {:.6f} -> {:.6f} ({:+.2f}%)",
        node->name(), vBefore, vSetPU, (vSetPU - vBefore) / vBefore * 100.0);
    mQLimitConvertedAtMax.erase(node);
    ++mQLimitSwitchCount[node];
    SPDLOG_LOGGER_INFO(mSLog, "Q-limit: PQ bus {} -> PV (constraint relaxed)",
                       node->name());
  }

  for (auto n : mSystem.mNodes) {
    UInt i = n->matrixNodeIndex();
    sol_V_complex(i) = Math::polar(sol_V(i), sol_D(i));
  }

  SPDLOG_LOGGER_INFO(mSLog, "Q-limit pass: {} PV -> PQ, {} PQ -> PV",
                     toPQ.size(), toPV.size());
  return !toPQ.empty() || !toPV.empty();
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
    for (CPS::SparseMatrixCompRow::InnerIterator it(mY, node_idx); it; ++it)
      I += it.value() * sol_Vcx(it.col());
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
