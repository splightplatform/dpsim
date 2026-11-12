// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim/PFSolverPowerPolarSparse.h>

#if defined(WITH_KLU)
#include <dpsim/KLUAdapter.h>
#elif defined(WITH_SPARSE)
#include <dpsim/SparseLUAdapter.h>
#else
#include <dpsim/DenseLUAdapter.h>
#endif

using namespace DPsim;
using namespace CPS;

PFSolverPowerPolarSparse::PFSolverPowerPolarSparse(
    CPS::String name, const CPS::SystemTopology &system, CPS::Real timeStep,
    CPS::Logger::Level logLevel)
    : PFSolverPowerPolar(name, system, timeStep, logLevel) {}

void PFSolverPowerPolarSparse::buildJacobianPattern() {
  /* 
  This function just creates the entries of the Jacobian matrix
  We are now using the reverse lookup table instead of scanning
  the whole array looking for connection
  */
  UInt npqpv = mNumPQBuses + mNumPVBuses;
  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.reserve(4 * mY.nonZeros());

  // Rows: P at a, Q at a + npqpv (PQ only).
  // Cols: theta at b, V at b + npqpv (PQ only)
  // Every stored entry of row k of mY is a connection to bus k.
  for (UInt a = 0; a < npqpv; ++a) {
    const UInt k = mPQPVBusIndices[a];
    const bool kIsPQ = a < mNumPQBuses;

    // diagonals: bus k's own rows and columns
    triplets.emplace_back(a, a, 0.0);
    if (kIsPQ) {
      triplets.emplace_back(a, a + npqpv, 0.0);
      triplets.emplace_back(a + npqpv, a, 0.0);
      triplets.emplace_back(a + npqpv, a + npqpv, 0.0);
    }

    // off-diagonals: one set per neighbour j
    // use inner-iterator to find connected
    for (CPS::SparseMatrixCompRow::InnerIterator it(mY, k); it; ++it) {
      const UInt j = static_cast<UInt>(it.col());
      const Int b = mPQPVPositionOfBus[j];
      if (j == k || b < 0) // diagonal handled above & slack has no column
        continue;
      const bool jIsPQ = static_cast<UInt>(b) < mNumPQBuses;

      triplets.emplace_back(a, b, 0.0); // J1
      if (jIsPQ)
        triplets.emplace_back(a, b + npqpv, 0.0); // J2
      if (kIsPQ)
        triplets.emplace_back(a + npqpv, b, 0.0); // J3
      if (kIsPQ && jIsPQ)
        triplets.emplace_back(a + npqpv, b + npqpv, 0.0); // J4
    }
  }

  mJsparse.resize(mNumUnknowns, mNumUnknowns);
  mJsparse.setFromTriplets(triplets.begin(), triplets.end());
  mJsparse.makeCompressed();
}

void PFSolverPowerPolarSparse::setUpJacobianStorage() {
  if (mNumUnknowns == 0)
    return;

  buildJacobianPattern();

#if defined(WITH_KLU)
  mLinearSolver = std::make_shared<KLUAdapter>(mSLog);
#elif defined(WITH_SPARSE)
  mLinearSolver = std::make_shared<SparseLUAdapter>(mSLog);
#else
  mLinearSolver = std::make_shared<DenseLUAdapter>(mSLog);
#endif

  // Analyze the fixed pattern once; values are refactorized each iteration.
  mLinearSolver->preprocessing(mJsparse, mVariableSystemMatrixEntries);
}

void PFSolverPowerPolarSparse::calculateJacobian() {
  /* 
  This function fills the entries of the Jacobian matrix
  We are now using the reverse lookup table instead of scanning
  the whole array looking for connection
  Same entries as the previous block-by-block fill, computed in one pass over each row of mY
  */
  const UInt npqpv = mNumPQBuses + mNumPVBuses;
  // Pattern is fixed; clear values and refill the existing entries in place.
  mJsparse.coeffs().setZero();

  for (UInt a = 0; a < npqpv; ++a) {
    const UInt k = mPQPVBusIndices[a];
    const bool kIsPQ = a < mNumPQBuses;
    const Real Vk = sol_V.coeff(k);
    Real Gkk = 0.0, Bkk = 0.0;

    for (CPS::SparseMatrixCompRow::InnerIterator it(mY, k); it; ++it) {
      const UInt j = static_cast<UInt>(it.col());
      if (j == k) {
        Gkk = it.value().real();
        Bkk = it.value().imag();
        continue;
      }
      const Int b = mPQPVPositionOfBus[j]; 
      if (b < 0) // slack 
        continue; 
      const bool jIsPQ = static_cast<UInt>(b) < mNumPQBuses; 

      const Real G = it.value().real(); 
      const Real B = it.value().imag(); 
      const Real d = sol_D.coeff(k) - sol_D.coeff(j);
      const Real VkVj = Vk * sol_V.coeff(j); 
      const Real s = VkVj * (G * sin(d) - B * cos(d)); 
      const Real c = VkVj * (G * cos(d) + B * sin(d)); 

      mJsparse.coeffRef(a, b) = s; // J1 dP/dTheta
      if (jIsPQ)
        mJsparse.coeffRef(a, b + npqpv) = c; // J2 dP/dV
      if (kIsPQ)
        mJsparse.coeffRef(a + npqpv, b) = -c; // J3 dQ/dTheta
      if (kIsPQ && jIsPQ)
        mJsparse.coeffRef(a + npqpv, b + npqpv) = s; // J4 dQ/dV
    }
      
    // diagonals, with P(k) and Q(k) computed once
    const Real Pk = P(k);
    const Real Qk = Q(k);
    mJsparse.coeffRef(a, a) = -Qk - Bkk * Vk * Vk; // J1
    if (kIsPQ) {
      mJsparse.coeffRef(a, a + npqpv) = Pk + Gkk * Vk * Vk;         // J2
      mJsparse.coeffRef(a + npqpv, a) = Pk - Gkk * Vk * Vk;         // J3
      mJsparse.coeffRef(a + npqpv, a + npqpv) = Qk - Bkk * Vk * Vk; // J4
    }
  }

  // Current-balance rows
  if (mCurrentBalanceRows) {
    for (UInt a = 0; a < mNumPQBuses; ++a) {
      const UInt k = mPQPVBusIndices[a];
      const Real voltage = sol_V.coeff(k);
      // same condition as calculateMismatch, so rows and Jacobian stay consistent
      if (voltage <= 0.0)
        continue;
      for (const UInt row : {a, a + npqpv})
        for (CPS::SparseMatrixRow::InnerIterator it(mJsparse, row); it; ++it)
          it.valueRef() /= voltage;
      mJsparse.coeffRef(a, a + npqpv) += mF(a);
      mJsparse.coeffRef(a + npqpv, a + npqpv) += mF(a + npqpv);
    }
  }
}

void PFSolverPowerPolarSparse::solveJacobianSystem() {
  // Full numeric factorization with pivoting on the first iteration of each run;
  // refactorize (reuse ordering + symbolic) afterwards. mIterations is 0 only on
  // the first iteration (set after the solve in PFSolver::solvePowerflow).
  if (mIterations == 0)
    mLinearSolver->factorize(mJsparse);
  else
    mLinearSolver->refactorize(mJsparse);

  Matrix rhs = mF;
  mX = mLinearSolver->solve(rhs);
}

CPS::Bool PFSolverPowerPolarSparse::remoteRegulationSensitivities(
      const std::vector<CPS::UInt> &inputs, const std::vector<CPS::UInt> &vOut,
      const std::vector<CPS::UInt> &qOut, CPS::Matrix &dV, CPS::Matrix &dQ) {
  
  // inputs is a list of bus indices (PV buses) whose terminal setpoints the Newton loop
  // will move this pass. An empty call is the capability probe.
  if (inputs.empty())
    return true;

  calculateMismatch();
  calculateJacobian();
  mLinearSolver->factorize(mJsparse);

  // build rhs matrix (one column per input in 64 size chunks)
  const UInt n = inputs.size();
  const UInt npqpv = mNumPQBuses + mNumPVBuses;
  const Int nPQ = static_cast<Int>(mNumPQBuses);

  dV = Matrix::Zero(vOut.size(), n);
  dQ = Matrix::Zero(qOut.size(), n);
  // P and Q at each Q output, needed for the diagonal terms of every column
  std::vector<Real> pOut(qOut.size()), qOutVal(qOut.size());
  for (UInt i = 0; i < qOut.size(); ++i) {
    pOut[i] = P(qOut[i]);
    qOutVal[i] = Q(qOut[i]);
  }

  for (UInt start = 0; start < n; start += 64){
    UInt cols = std::min<UInt>(64, n - start);
    Matrix rhs = Matrix::Zero(mNumUnknowns, cols);
    std::vector<bool> valid(cols, false);
    for (UInt c = 0; c < cols; ++c){
      UInt g = inputs[start + c];
      const Int pg = mPQPVPositionOfBus[g]; // position index for bus g
      if (pg < (Int)mNumPQBuses) {
        // PQ (pinned at a Q limit) or slack: its setpoint isn't a free input
        SPDLOG_LOGGER_WARN(mSLog,
                           "Remote regulation sensitivities: input bus idx {} is "
                           "not PV (it is {}); its column is left at zero",
                           g, pg < 0 ? "the slack" : "PQ");
        continue;
      }
      valid[c] = true;
      const Real Vg = sol_V(g);
      rhs(pg, c) = -(P(g) / Vg + mY.coeff(g, g).real() * Vg); // own P row
      for (CPS::SparseMatrixCompRow::InnerIterator it(mY, g); it; ++it) { // gathers (x,y) pairs only when they exist for g
        const UInt k = it.col();
        if (k == g) continue; // own bus done above
        const Int pk = mPQPVPositionOfBus[k]; // position index for bus k
        if (pk < 0) continue; // slack: no rows
        const bool kIsPQ = pk < (Int)mNumPQBuses;
        const Complex Ykg = mY.coeff(k, g); // row k: P_k depends on bus k's row of Y
        const Real Vk = sol_V(k), d = sol_D(k) - sol_D(g);
        Real dPk = -Vk * (Ykg.real() * cos(d) + Ykg.imag() * sin(d));
        Real dQk = -Vk * (Ykg.real() * sin(d) - Ykg.imag() * cos(d));
        if (kIsPQ && mCurrentBalanceRows && Vk > 0.0) { 
          dPk /= Vk; 
          dQk /= Vk; 
        }
        rhs(pk, c) = dPk; // B: neighbour P row
        if (kIsPQ) rhs(pk + npqpv, c) = dQk;  // C: neighbour Q row
      }
    }
    // x column c: change of every theta and dV/V per +1 pu of input c's setpoint
    Matrix x = mLinearSolver->solve(rhs);

    for (UInt c = 0; c < cols; ++c) {
      if (!valid[c])
        continue;
      const UInt j = start + c;
      const UInt g = inputs[j];
      // angle change at a bus; the slack holds its angle
      auto dTh = [&](UInt b) -> Real {
        const Int p = mPQPVPositionOfBus[b];
        return p < 0 ? 0.0 : x(p, c);
      };
      // |V| change at a bus: 1 at the moved input, x*V at PQ buses,
      // 0 at other PV buses and the slack (they hold their voltage)
      auto dVm = [&](UInt b) -> Real {
        if (b == g)
          return 1.0;
        const Int p = mPQPVPositionOfBus[b];
        return (p >= 0 && p < nPQ) ? x(p + npqpv, c) * sol_V(b) : 0.0;
      };

      for (UInt i = 0; i < vOut.size(); ++i)
        dV(i, j) = dVm(vOut[i]);

      // Q_h = V_h sum_j V_j (G sin d - B cos d), d = theta_h - theta_j.
      // Loads are constant power, so dQ of the generator = dQ of the net injection.
      for (UInt i = 0; i < qOut.size(); ++i) {
        const UInt h = qOut[i];
        const Real Vh = sol_V(h);
        Real Ghh = 0.0, Bhh = 0.0, sum = 0.0;
        for (CPS::SparseMatrixCompRow::InnerIterator it(mY, h); it; ++it) {
          const UInt k = static_cast<UInt>(it.col());
          const Real G = it.value().real(), B = it.value().imag();
          if (k == h) {
            Ghh = G;
            Bhh = B;
            continue;
          }
          const Real d = sol_D(h) - sol_D(k);
          const Real Vk = sol_V(k);
          sum += -Vh * Vk * (G * cos(d) + B * sin(d)) * dTh(k) +
                 Vh * (G * sin(d) - B * cos(d)) * dVm(k);
        }
        sum += (pOut[i] - Ghh * Vh * Vh) * dTh(h) +
               (qOutVal[i] / Vh - Bhh * Vh) * dVm(h);
        dQ(i, j) = sum;
      }
    }
  }
  return true;
}