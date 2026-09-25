#pragma once
#include <array>
#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Base/Base_Ph1_Transformer3W.h>
#include <dpsim-models/Definitions.h>

namespace CPS {
namespace Base {
namespace Ph3 {


using Winding3W = CPS::Base::Ph1::Winding3W;
using CPS::Base::Ph1::AllWindings3W;

// How a winding is physically connected to its three line terminals.
//
// Wye:   the star-equivalent topology SP::Ph1::Transformer3W already
//        implements (phase shift encoded in turns ratio for PF)
// Delta: three legs (AB, BC, CA) each span two of
//        this winding's own terminals, and terminal currents are the
//        superposition of two leg currents. This creates a phase shift
//        in EMT
//
// phaseShift [deg] is only meaningful when connection == Delta 
// it selects which two lines a given leg bridges and with what polarity

enum class WindingConnection { Wye, Delta };

class Transformer3W {
public:
  static constexpr UInt NumWindings = 3; // don't change

protected:
  // Nominal (rated) voltage of each winding [V], indexed by Winding3W
  std::array<Real, NumWindings> mNominalVoltage{{0., 0., 0.}};

  // Winding the star-branch impedances are referred to.
  // Resolved in resolveReferenceWinding().
  Winding3W mReferenceWinding = Winding3W::Primary;
  Bool mReferenceWindingSet = false;

  // Wye by default; set per winding
  std::array<WindingConnection, NumWindings> mWindingConnection{
      {WindingConnection::Wye, WindingConnection::Wye,
       WindingConnection::Wye}};
  // Phase shift [deg] relative to the reference winding
  // use when corresponding winding's connection is Delta
  std::array<Real, NumWindings> mPhaseShift{{0., 0., 0.}};

  // Star-branch (wye) or leg (delta) resistance/inductance of each winding,
  // per phase [Ohm]/[H], referred to the reference winding
  std::array<Matrix, NumWindings> mResistance{
      {Matrix::Zero(3, 3), Matrix::Zero(3, 3), Matrix::Zero(3, 3)}};
  std::array<Matrix, NumWindings> mInductance{
      {Matrix::Zero(3, 3), Matrix::Zero(3, 3), Matrix::Zero(3, 3)}};

  void resolveReferenceWinding() {
    if (mReferenceWindingSet)
      return;
    UInt best = 0;
    for (UInt i = 1; i < NumWindings; i++)
      if (mNominalVoltage[i] > mNominalVoltage[best])
        best = i;
    mReferenceWinding = static_cast<Winding3W>(best);
  }

public:
  // Rated apparent power of each winding [VA]
  const Attribute<Real>::Ptr mRatedPowerP, mRatedPowerS, mRatedPowerT;
  const Attribute<Complex>::Ptr mRatioP, mRatioS, mRatioT;

  explicit Transformer3W(CPS::AttributeList::Ptr attributeList)
      : mRatedPowerP(attributeList->create<Real>("S_p")),
        mRatedPowerS(attributeList->create<Real>("S_s")),
        mRatedPowerT(attributeList->create<Real>("S_t")),
        mRatioP(attributeList->create<Complex>("ratio_p")),
        mRatioS(attributeList->create<Complex>("ratio_s")),
        mRatioT(attributeList->create<Complex>("ratio_t")){};

  // Accessors

  Real nominalVoltage(Winding3W w) const {
    return mNominalVoltage[static_cast<UInt>(w)];
  }
  Winding3W referenceWinding() const { return mReferenceWinding; }
  Bool isReferenceWinding(Winding3W w) const {
    return w == mReferenceWinding;
  }

  Real ratedPower(Winding3W w) const {
    switch (w) {
    case Winding3W::Primary:
      return **mRatedPowerP;
    case Winding3W::Secondary:
      return **mRatedPowerS;
    default:
      return **mRatedPowerT;
    }
  }
  Complex ratio(Winding3W w) const {
    switch (w) {
    case Winding3W::Primary:
      return **mRatioP;
    case Winding3W::Secondary:
      return **mRatioS;
    default:
      return **mRatioT;
    }
  }
  const Matrix &resistance(Winding3W w) const {
    return mResistance[static_cast<UInt>(w)];
  }
  const Matrix &inductance(Winding3W w) const {
    return mInductance[static_cast<UInt>(w)];
  }
  WindingConnection windingConnection(Winding3W w) const {
    return mWindingConnection[static_cast<UInt>(w)];
  }
  Real phaseShift(Winding3W w) const {
    return mPhaseShift[static_cast<UInt>(w)];
  }


  void setReferenceWinding(Winding3W w) {
    mReferenceWinding = w;
    mReferenceWindingSet = true;
  }

  // Marks a winding Delta-connected at the given phase-shift [deg] relative
  // to the reference winding Leaving a winding defaults to wye 
  void setWindingConnection(Winding3W w, WindingConnection connection,
                            Real phaseShiftDeg = 0.) {
    mWindingConnection[static_cast<UInt>(w)] = connection;
    mPhaseShift[static_cast<UInt>(w)] = phaseShiftDeg;
  }

  // Set parameters directly from the star equivalent. R, L are per phase,
  // already referred to the reference winding 
  void setParameters(const std::array<Real, NumWindings> &nomVoltage,
                     const std::array<Real, NumWindings> &ratedPower,
                     const std::array<Real, NumWindings> &ratioAbs,
                     const std::array<Real, NumWindings> &ratioPhase,
                     const std::array<Matrix, NumWindings> &resistance,
                     const std::array<Matrix, NumWindings> &inductance) {
    mNominalVoltage = nomVoltage;
    resolveReferenceWinding();

    **mRatedPowerP = ratedPower[0];
    **mRatedPowerS = ratedPower[1];
    **mRatedPowerT = ratedPower[2];

    **mRatioP = std::polar<Real>(ratioAbs[0], ratioPhase[0]);
    **mRatioS = std::polar<Real>(ratioAbs[1], ratioPhase[1]);
    **mRatioT = std::polar<Real>(ratioAbs[2], ratioPhase[2]);

    mResistance = resistance;
    mInductance = inductance;
  }
};
} // namespace Ph3
} // namespace Base
} // namespace CPS
