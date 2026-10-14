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
