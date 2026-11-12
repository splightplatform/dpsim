/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <map>
#include <limits>

#include <dpsim/PFSolver.h>

namespace DPsim {
/// Powerflow solver class considering power mismatch and voltages in polar coordinates.
class PFSolverPowerPolar : public PFSolver {
protected:
  /// Solution vector of active power
  CPS::Vector sol_P;
  /// Solution vector of reactive power
  CPS::Vector sol_Q;
  /// Solution vector of voltage magnitude
  CPS::Vector sol_V;
  /// Solution vector of voltage angle
  CPS::Vector sol_D;
  /// Solution vector of representing sol_V and sol_D as complex quantity
  CPS::VectorComp sol_V_complex;
  /// Solution vector of representing sol_P and sol_Q as complex quantity
  CPS::VectorComp sol_S_complex;

  CPS::Vector Pesp;
  CPS::Vector Qesp;

  // Cache for last converged solution
  CPS::Vector mLastConvergedV;
  CPS::Vector mLastConvergedD;
  bool mHasLastConvergedSolution = false;

  // Cache best iterate of the run for reload if non-convergence
  CPS::Vector mBestRunV; 
  CPS::Vector mBestRunD;
  CPS::Real mBestRunNorm = std::numeric_limits<CPS::Real>::infinity(); 

  // Solution state of the saved operating point
  CPS::Vector mSavedV, mSavedD, mSavedP, mSavedQ, mSavedPesp, mSavedQesp;
  std::map<CPS::TopologicalNode::Ptr, bool> mSavedQLimitConvertedAtMax;
  std::map<CPS::TopologicalNode::Ptr, CPS::UInt> mSavedQLimitSwitchCount;

  // counter & limits for iteration improvement tracking
  CPS::UInt mNumNonImprovingSteps = 0; 
  CPS::UInt mMaxNonImprovingSteps = 8; 
  CPS::Real mBlowUpFactor = 10.0; 

  /// Divide each PQ bus's mismatch rows by |V| (current balance) instead of power balance
  CPS::Bool mCurrentBalanceRows = true;
  /// Largest angle step [rad] after whole-step scaling
  CPS::Real mDThetaLimit = 0.2;
  /// Largest relative voltage step dV/V after whole-step scaling
  CPS::Real mDVLimit = 0.1;

  // Remote bus regulation
  CPS::Real busVoltageMagnitude(CPS::UInt busIdx) override { return sol_V(busIdx); }
  void setBusVoltageMagnitude(CPS::UInt busIdx, CPS::Real value) override { sol_V(busIdx) = value; }
  CPS::Bool isQLimitPinned(CPS::UInt busIdx) override {
    for (auto &node : mSystem.mNodes)
      if (node->matrixNodeIndex() == busIdx)
        return mQLimitConvertedAtMax.count(node) > 0;
    return false;
  }

  // Core methods
  /// Generate initial solution for current time step
  void generateInitialSolution(Real time,
                               bool keep_last_solution = false) override;
  /// Calculate the Jacobian
  void calculateJacobian() override;
  /// Update solution in each iteration
  Bool updateSolution() override;
  /// Set final solution
  void setSolution() override;
  // Set the best values at the beginning of the Newton run
  void startNewtonRun() override;
  // Check if Newton run converged, if not, reset state to best values
  void finishNewtonRun() override;
  /// Calculate mismatch
  void calculateMismatch() override;
  /// Save / restore the operating point, solution and Q-limit state included
  void saveOperatingPoint() override;
  void restoreOperatingPoint() override;

  // Helper methods
  /// Resize solution vector
  void resize_sol(CPS::Int n);
  /// Resize complex solution vector
  void resize_complex_sol(CPS::Int n);
  /// Calculate real part of voltage from sol_V and sol_D
  CPS::Real sol_Vr(CPS::UInt k);
  /// Calculate imaginary part of voltage from sol_V and sol_D
  CPS::Real sol_Vi(CPS::UInt k);
  /// Calculate complex voltage from sol_V and sol_D
  CPS::Complex sol_Vcx(CPS::UInt k);
  /// Calculate active power at a bus from current solution
  CPS::Real P(CPS::UInt k);
  /// Calculate the reactive power at a bus from current solution
  CPS::Real Q(CPS::UInt k);
  /// Calculate P and Q at slack bus from current solution
  void calculatePAndQAtSlackBus();
  /// Calculate the reactive power at all PV buses from current solution
  void calculateQAtPVBuses();
  /// Generator reactive-power injection at a bus [pu], used by the Q-limit check
  CPS::Real generatorReactivePowerPerUnit(CPS::TopologicalNode::Ptr node);
  /// Total local load reactive power at a bus [pu]
  CPS::Real loadReactivePowerPerUnit(CPS::TopologicalNode::Ptr node);
  /// Q-limit PV<->PQ switching pass (overrides the base no-op)
  CPS::Bool enforceReactiveLimits() override;
  /// Clear the Q-limit conversion bookkeeping between solves
  void clearReactiveLimitState() override;
  /// Buses switched PV->PQ by the Q-limit loop -> pinned at Qmax (true) or Qmin (false)
  std::map<CPS::TopologicalNode::Ptr, bool> mQLimitConvertedAtMax;
  /// Total PV<->PQ switches per bus over the solve (anti-oscillation cap)
  std::map<CPS::TopologicalNode::Ptr, CPS::UInt> mQLimitSwitchCount;
  /// Calculate complex power flowing from this node to the other nodes
  void calculatePAndQInjectionPQBuses();
  /// Calculate branch flows from current solution and store them in line and transformer components
  void calculateBranchFlow();
  /// Calculate nodal power injections and store them in first line or transformer (in case no line is connected), so that lower level classes (Node, TopologicalTerminal) can stay untouched
  void calculateNodalInjection();

public:
  /// Constructor to be used in simulation examples.
  PFSolverPowerPolar(CPS::String name, const CPS::SystemTopology &system,
                     CPS::Real timeStep, CPS::Logger::Level logLevel);
  ///
  virtual ~PFSolverPowerPolar(){};

  /// Use current-balance rows for PQ buses (true) or plain power-balance rows (false)
  void setCurrentBalanceRows(CPS::Bool value) { mCurrentBalanceRows = value; }
  /// Largest angle step [rad] per Newton iteration
  void setDThetaLimit(CPS::Real value) { mDThetaLimit = value; }
  /// Largest relative voltage step dV/V per Newton iteration
  void setDVLimit(CPS::Real value) { mDVLimit = value; }
};
} // namespace DPsim
