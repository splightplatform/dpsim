/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <cmath>
#include <iterator>
#include <map>

#include "dpsim-models/Components.h"
#include "dpsim-models/SystemTopology.h"
#include <dpsim/Scheduler.h>
#include <dpsim/Solver.h>
#include <functional>

namespace DPsim {
/// Snapshot passed to the progress callback after each mismatch evaluation
struct PFProgress {
  CPS::String phase;          // label of the current Newton run
  CPS::UInt iteration;        // within this run; 0 = before the first step
  CPS::UInt totalIterations;  // across the whole solve
  CPS::Real maxMismatch;      // max |mismatch|, per unit
  CPS::Real mismatchNorm;     // 2-norm of the mismatch, per unit
};

using PFProgressCallback = std::function<void(const PFProgress &)>;
/// Solver class using the nonlinear powerflow (PF) formulation.
class PFSolver : public Solver {
protected:
  PFProgressCallback mProgressCallback;
  void reportProgress(const CPS::String &phase, CPS::UInt iteration);

  // Levenberg-Marquardt regularization for the Newton-Raphson solve
  Real mLmLambda = 1e-3;

  // Remote regulation
  std::map<CPS::UInt, CPS::UInt> mRegulatedBusOfGen; // genBusIdx -> targetBusIdx (only entries where they differ)
  std::map<CPS::UInt, CPS::Real> mLocalVSetOverride; // genBusIdx -> current local V_set_pu being tuned
  std::map<CPS::UInt, CPS::Real> mRegulatedVSetPU;   // genBusIdx -> target pu at the remote bus
  CPS::Real mRemoteRegTolerance = 1e-5;
  CPS::UInt mMaxRemoteRegIterations = 30;

  /// Drive local PV setpoints so remote-regulated buses reach their target voltage
  Bool resolveRemoteRegulation();
  /// Current solved |V| at a bus, in the solver's own representation.
  virtual CPS::Real busVoltageMagnitude(CPS::UInt busIdx) = 0;
  virtual void setBusVoltageMagnitude(CPS::UInt busIdx, CPS::Real value) = 0;
  virtual CPS::Bool isQLimitPinned(CPS::UInt busIdx) { return false; } // base: no Q-limit concept

  /// Number of PQ nodes
  UInt mNumPQBuses = 0;
  /// Number of PV nodes
  UInt mNumPVBuses = 0;
  /// Number of PV nodes
  UInt mNumVDBuses = 0;
  /// Number of unknowns, defining system dimension
  UInt mNumUnknowns = 0;
  /// Vector of nodes characterized as PQ buses
  CPS::TopologicalNode::List mPQBuses;
  /// Vector of nodes characterized as PV buses
  CPS::TopologicalNode::List mPVBuses;
  /// Vector of nodes characterized as VD buses
  CPS::TopologicalNode::List mVDBuses;
  /// Original PQ/PV classification (snapshot before Q-limit switching)
  CPS::TopologicalNode::List mPQBusesOrig;
  CPS::TopologicalNode::List mPVBusesOrig;
  /// Vector with indices of PQ buses
  std::vector<CPS::UInt> mPQBusIndices;
  /// Vector with indices of PV buses
  std::vector<CPS::UInt> mPVBusIndices;
  /// Vector with indices of VD buses
  std::vector<CPS::UInt> mVDBusIndices;
  /// Vector with indices of both PQ and PV buses
  std::vector<CPS::UInt> mPQPVBusIndices;
  /// Position of each bus in mPQPVBusIndices for reverse lookup, -1 for VD buses
  std::vector<Int> mPQPVPositionOfBus; 

  /// Admittance matrix
  CPS::SparseMatrixCompRow mY;

  /// Jacobian matrix
  CPS::Matrix mJ;
  /// Solution vector
  CPS::Vector mX;
  /// Vector of mismatch values
  CPS::Vector mF;

  /// System list
  CPS::SystemTopology mSystem;
  /// Vectors of transformer components
  std::vector<std::shared_ptr<CPS::SP::Ph1::Transformer>> mTransformers;
  std::vector<std::shared_ptr<CPS::SP::Ph1::Transformer3W>> mTransformers3W; 
  /// Vector of solid state transformer components
  std::vector<std::shared_ptr<CPS::SP::Ph1::SolidStateTransformer>>
      mSolidStateTransformers;
  /// Vector of synchronous generator components
  std::vector<std::shared_ptr<CPS::SP::Ph1::SynchronGenerator>>
      mSynchronGenerators;
  /// Vector of load components
  std::vector<std::shared_ptr<CPS::SP::Ph1::Load>> mLoads;
  /// Vector of line components
  std::vector<std::shared_ptr<CPS::SP::Ph1::PiLine>> mLines;
  /// Vector of shunt components
  std::vector<std::shared_ptr<CPS::SP::Ph1::Shunt>> mShunts;
  /// Vector of external grid components
  std::vector<std::shared_ptr<CPS::SP::Ph1::NetworkInjection>> mExternalGrids;
  /// Vector of average voltage source inverters
  std::vector<std::shared_ptr<CPS::SP::Ph1::AvVoltageSourceInverterDQ>>
      mAverageVoltageSourceInverters;
  /// Map providing determined base voltages for each node
  std::map<CPS::TopologicalNode::Ptr, CPS::Real> mBaseVoltageAtNode;

  /// Solver tolerance
  Real mTolerance = 1e-8;
  /// Maximum number of iterations
  CPS::UInt mMaxIterations = 20;
  /// Actual number of iterations of the last Newton-Raphson run
  CPS::UInt mIterations = 0;
  /// Newton-Raphson iterations summed over every run of the last solve
  CPS::UInt mTotalIterations = 0;
  /// Whether the Q-limit and remote-regulation outer loops of the last solve
  /// reached their fixed point (true when neither loop is enabled)
  CPS::Bool mOuterLoopsSettled = false;
  /// Enforce generator reactive-power limits via PV<->PQ outer-loop switching
  CPS::Bool mEnforceReactiveLimits = false;
  /// Maximum number of Q-limit outer iterations
  CPS::UInt mMaxOuterIterations = 10;
  /// Maximum number of PV<->PQ switches per bus before it is frozen (anti-oscillation)
  CPS::UInt mMaxQLimitSwitchesPerBus = 2;
  /// Run fast-decoupled iterations before the initial Newton-Raphson solve
  CPS::Bool mUseFastDecoupledWarmStart = false;
  /// Iteration cap for the fast-decoupled warm start
  CPS::UInt mFastDecoupledMaxIterations = 20;
  /// Max |mismatch| [pu] at which the warm start hands over to Newton-Raphson
  CPS::Real mFastDecoupledHandoffTol = 1e-2;
  /// Relative tolerance for non-authoritative (e.g. Load) base-voltage candidates vs. the zone's rating
  CPS::Real mBaseVoltageLooseTolerance = 0.1;
  /// Relative tolerance between authoritative base-voltage sources (generator/transformer/network-injection/VSI) within a zone
  CPS::Real mBaseVoltageStrictTolerance = 0.01;
  /// Base power of per-unit system
  CPS::Real mBaseApparentPower;
  /// Fallback base apparent power if no generator or transformer rating is found
  CPS::Real mBaseApparentPowerFallback = 100e6;
  /// Convergence flag
  CPS::Bool isConverged = false;
  /// Flag whether solution vectors are initialized
  CPS::Bool solutionInitialized = false;
  /// Flag whether complex solution vectors are initialized
  CPS::Bool solutionComplexInitialized = false;

  /// Use last converged solution as initial guess
  CPS::Bool mKeepLastSolution = false;

  /// Generate initial solution for current time step
  virtual void generateInitialSolution(Real time,
                                       bool keep_last_solution = false) = 0;
  /// Calculate mismatch
  virtual void calculateMismatch() = 0;
  /// Calculate the Jacobian
  virtual void calculateJacobian() = 0;
  /// Update solution in each iteration
  // true means a step was accepted (went correct direction)
  // false means the previous state was restored and the run should stop (divergence)
  virtual Bool updateSolution() = 0;
  /// Set final solution
  virtual void setSolution() = 0;

  /// Initialization of the solver
  void initialize() override;
  /// Initialization of individual components
  void initializeComponents();
  /// Assignment of matrix indices for nodes
  void assignMatrixNodeIndices();
  /// Set apparent base power of per-unit system
  void setBaseApparentPower();
  /// Determine bus type for all buses
  void determinePFBusType();
  /// Rebuild index vectors + counts from the PQ/PV/VD node lists
  void rebuildBusIndexAggregates();
  /// Re-derive index vectors and resize the system after PV<->PQ switching
  void reclassifyBuses();
  /// Restore the original PV/PQ classification before a fresh solve
  void resetToOriginalClassification();
  /// Clear Q-limit bookkeeping; overridden by PFSolverPowerPolar
  virtual void clearReactiveLimitState() {}
  /// Base voltage a single component reports for `node`, or 0 if unknown
  CPS::Real componentBaseVoltage(CPS::TopologicalPowerComp::Ptr comp,
                                 CPS::TopologicalNode::Ptr node);
  /// Determine, verify and propagate each node's base voltage per electrical zone
  void propagateAndVerifyBaseVoltage();

  /// Compose admittance matrix
  void composeAdmittanceMatrix();
  /// Gets the real part of admittance matrix element
  CPS::Real G(int i, int j);
  /// Gets the imaginary part of admittance matrix element
  CPS::Real B(int i, int j);
  /// Solves the powerflow problem
  Bool solvePowerflow();
  /// Run a single Newton-Raphson solve with the current bus classification
  Bool runNewtonRaphson(const CPS::String &label = "");
  /// Switch generators violating their Q limits between PV/PQ; base impl is a no-op
  virtual CPS::Bool enforceReactiveLimits() { return false; }
  /// Fast-decoupled warm start ahead of Newton-Raphson; true if it reached
  /// handoffTol. Base impl is a no-op
  virtual CPS::Bool runFastDecoupled(CPS::UInt maxIter, CPS::Real handoffTol) {
    return false;
  }
  /// Allocate Jacobian storage; dense by default, sparse subclass overrides
  virtual void setUpJacobianStorage();
  /// Solve the linearized system mJ*mX = mF into mX; sparse subclass overrides
  virtual void solveJacobianSystem();
  /// Check whether below tolerance
  CPS::Bool checkConvergence();
  /// Logging for integer vectors
  CPS::String logVector(std::vector<CPS::UInt> indexVector) {
    std::stringstream result;
    std::copy(indexVector.begin(), indexVector.end(),
              std::ostream_iterator<CPS::UInt>(result, " "));
    return result.str();
  };
  ///
  CPS::Task::List getTasks() override;
  // determines power flow bus type for each node according to the components attached to it.
public:
  void setProgressCallback(PFProgressCallback cb) { mProgressCallback = std::move(cb); }
  // Accessors for more observability during load-flow solved
  /// Whether the last Newton-Raphson run converged
  CPS::Bool getConverged() const { return isConverged; }
  /// Iterations of the last Newton-Raphson run
  CPS::UInt getIterations() const { return mIterations; }
  /// Iterations summed over every Newton-Raphson run of the last solve,
  /// including the Q-limit and remote-regulation outer loops
  CPS::UInt getTotalIterations() const { return mTotalIterations; }
  /// Whether the Q-limit and remote-regulation outer loops of the last solve
  /// settled; false if either gave up or was skipped
  CPS::Bool getOuterLoopsSettled() const { return mOuterLoopsSettled; }
  /// Per-unit bus admittance matrix, indexed by matrixNodeIndex()
  const CPS::SparseMatrixCompRow &getAdmittanceMatrix() const { return mY; }
  /// Base power of the per-unit system [VA]
  CPS::Real getBaseApparentPower() const { return mBaseApparentPower; }
  /// Base voltage the solver assigned to a node [V]
  CPS::Real getBaseVoltageAtNode(const CPS::TopologicalNode::Ptr &node) const {
    auto it = mBaseVoltageAtNode.find(node);
    if (it == mBaseVoltageAtNode.end())
      throw std::out_of_range("No base voltage for node " + node->name());
    return it->second;
  }
  /// Constructor to be used in simulation examples.
  PFSolver(CPS::String name, CPS::SystemTopology system, Real timeStep,
           CPS::Logger::Level logLevel);
  ///
  virtual ~PFSolver(){};

  /// Set a node to VD using its name
  void setVDNode(CPS::String name);
  /// Allows to modify the powerflow bus type of a specific component
  void modifyPowerFlowBusComponent(CPS::String name,
                                   CPS::PowerflowBusType powerFlowBusType);
  /// set solver and component to initialization or simulation behaviour
  void setSolverAndComponentBehaviour(Solver::Behaviour behaviour) override;

  void setKeepLastSolution(CPS::Bool keepLastSolution) {
    mKeepLastSolution = keepLastSolution;
  }

  /// Enable generator reactive-limit enforcement (PV<->PQ outer loop)
  void setEnforceReactiveLimits(CPS::Bool value) {
    mEnforceReactiveLimits = value;
  }

  /// Run a fast-decoupled warm start before the initial Newton-Raphson solve
  void setUseFastDecoupledWarmStart(CPS::Bool value) {
    mUseFastDecoupledWarmStart = value;
  }

  CPS::Bool getUseFastDecoupledWarmStart() const {
    return mUseFastDecoupledWarmStart;
  }

  void setFastDecoupledMaxIterations(CPS::UInt maxIterations) {
    mFastDecoupledMaxIterations = maxIterations;
  }

  /// Max |mismatch| [pu] at which the warm start hands over to Newton-Raphson
  void setFastDecoupledHandoffTolerance(CPS::Real tolerance) {
    mFastDecoupledHandoffTol = tolerance;
  }

  /// Raise for grids with legitimate large voltage drop (e.g. untapped feeders)
  void setBaseVoltageLooseTolerance(CPS::Real tolerance) {
    mBaseVoltageLooseTolerance = tolerance;
  }

  /// Override the tolerance between authoritative base-voltage sources within a zone
  void setBaseVoltageStrictTolerance(CPS::Real tolerance) {
    mBaseVoltageStrictTolerance = tolerance;
  }

  CPS::Bool getKeepLastSolution() const { return mKeepLastSolution; }

  void setBaseApparentPowerFallback(CPS::Real baseApparentPowerFallback) {
    mBaseApparentPowerFallback = baseApparentPowerFallback;
  }

  CPS::Real getBaseApparentPowerFallback() const {
    return mBaseApparentPowerFallback;
  }

  void setMaxIterations(CPS::UInt maxIterations) {
    mMaxIterations = maxIterations;
  }

  CPS::UInt getMaxIterations() const { return mMaxIterations; }

  class SolveTask : public CPS::Task {
  public:
    SolveTask(PFSolver &solver)
        : Task(solver.mName + ".Solve"), mSolver(solver) {
      mModifiedAttributes.push_back(Scheduler::external);
    }

    void execute(Real time, Int timeStepCount);

  private:
    PFSolver &mSolver;
  };
};
} // namespace DPsim
