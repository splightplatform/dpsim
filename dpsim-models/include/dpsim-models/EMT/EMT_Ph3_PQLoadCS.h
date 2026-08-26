/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/EMT/EMT_Ph3_CurrentSource.h>
#include <dpsim-models/Solver/MNAInterface.h>

namespace CPS {
namespace EMT {
namespace Ph3 {
/// \brief 3-Phase PQ load represented by a current source
///
/// Unlike RXLoad, this load recomputes its current injection from mActivePower
/// and mReactivePower on every time step, so the setpoint can be updated at
/// runtime (e.g. **load->mActivePower = newP) without triggering a system
/// matrix refactorization.
class PQLoadCS : public CompositePowerComp<Real>,
                 public SharedFactory<PQLoadCS> {
protected:
  /// Internal current source
  std::shared_ptr<EMT::Ph3::CurrentSource> mSubCurrentSource;
  Real mVoltageRefAngle = 0;
  /// Recomputes the current reference from the present P/Q setpoint
  void updateSetPoint();
  void updateIntfValues();

public:
  /// Active power [Watt]
  const Attribute<Real>::Ptr mActivePower;
  /// Reactive power [VAr]
  const Attribute<Real>::Ptr mReactivePower;
  /// Nominal voltage [V]
  const Attribute<Real>::Ptr mNomVoltage;

  /// Defines UID, name and logging level
  PQLoadCS(String uid, String name,
           Logger::Level logLevel = Logger::Level::off);
  /// Defines name and logging level
  PQLoadCS(String name, Logger::Level logLevel = Logger::Level::off);
  /// Defines name, component parameters and logging level
  PQLoadCS(String name, Real activePower, Real reactivePower, Real volt,
           Logger::Level logLevel = Logger::Level::off);

  SimPowerComp<Real>::Ptr clone(String name) override;

  // information print
  virtual String description() override {
    return fmt::format("Active: {}MW, Reactive: {}MVAr, Voltage: {}kV",
                       **mActivePower * 1e-6, **mReactivePower * 1e-6,
                       **mNomVoltage * 1e-3);
  };
  ///
  void setParameters(Real activePower, Real reactivePower, Real nomVolt);
  /// Constructs and registers MNA subcomponents; idempotent.
  void createSubComponents() override;
  /// Derives values from power flow data and pushes them to subcomponents
  void initializeParentFromNodesAndTerminals(Real frequency) override;

  // #### MNA section ####
  /// MNA pre and post step operations
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
