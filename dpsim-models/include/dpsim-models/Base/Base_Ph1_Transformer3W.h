#pragma once
#include <array>
#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Definitions.h>

namespace CPS {
namespace Base {
namespace Ph1 {

enum class Winding3W : UInt { Primary = 0, Secondary = 1, Tertiary = 2 };

/// Iteration helper so the winding loops below read as loops 
inline constexpr std::array<Winding3W, 3> AllWindings3W{
    {Winding3W::Primary, Winding3W::Secondary, Winding3W::Tertiary}};

// this class holds the  parameters shared by every domain's three-winding xfmr
// holds Z_p, Z_s, Z_t, and turns ratio and rating per winding

class Transformer3W {
public:
  static constexpr UInt NumWindings = 3; // don't change

protected:
  // Nominal (rated) voltage of each winding [V], indexed by Winding3W
  std::array<Real, NumWindings> mNominalVoltage{{0., 0., 0.}};

  // Winding the star-point impedances are referred to. Every Z_i, and the
  // star point's own voltage level is expressed on this windings side
  // Resolved in resolveReferenceWinding()
  Winding3W mReferenceWinding = Winding3W::Primary;

  // True once setReferenceWinding() is called
  Bool mReferenceWindingSet = false;

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
  // Complex turns ratio of each winding
  const Attribute<Complex>::Ptr mRatioP, mRatioS, mRatioT;
  // Star-branch resistance of each winding [Ohm], referred to the ref winding
  const Attribute<Real>::Ptr mResistanceP, mResistanceS, mResistanceT;
  // Star-branch inductance of each winding [H], referred to the ref winidng
  const Attribute<Real>::Ptr mInductanceP, mInductanceS, mInductanceT;

  // attributeList is DPsim jargon 
  explicit Transformer3W(CPS::AttributeList::Ptr attributeList)
      : mRatedPowerP(attributeList->create<Real>("S_p")),
        mRatedPowerS(attributeList->create<Real>("S_s")),
        mRatedPowerT(attributeList->create<Real>("S_t")),
        mRatioP(attributeList->create<Complex>("ratio_p")),
        mRatioS(attributeList->create<Complex>("ratio_s")),
        mRatioT(attributeList->create<Complex>("ratio_t")),
        mResistanceP(attributeList->create<Real>("R_p")),
        mResistanceS(attributeList->create<Real>("R_s")),
        mResistanceT(attributeList->create<Real>("R_t")),
        mInductanceP(attributeList->create<Real>("L_p")),
        mInductanceS(attributeList->create<Real>("L_s")),
        mInductanceT(attributeList->create<Real>("L_t")){};

  // Some accessors

  Real nominalVoltage(Winding3W w) const {
    return mNominalVoltage[static_cast<UInt>(w)];
  }
  Winding3W referenceWinding() const { return mReferenceWinding; }
  Bool isReferenceWinding(Winding3W w) const { return w == mReferenceWinding; }

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
  Real resistance(Winding3W w) const {
    switch (w) {
    case Winding3W::Primary:
      return **mResistanceP;
    case Winding3W::Secondary:
      return **mResistanceS;
    default:
      return **mResistanceT;
    }
  }
  Real inductance(Winding3W w) const {
    switch (w) {
    case Winding3W::Primary:
      return **mInductanceP;
    case Winding3W::Secondary:
      return **mInductanceS;
    default:
      return **mInductanceT;
    }
  }

  // Force the reference winding instead of taking the highest-nominal-voltage
  // default. Call before setParameters()
  void setReferenceWinding(Winding3W w) {
    mReferenceWinding = w;
    mReferenceWindingSet = true;
  }

  // Set parameters directly from the star equivalent
  // Arrays are indexed by Winding3W: [Primary, Secondary, Tertiary]
  void setParameters(const std::array<Real, NumWindings> &nomVoltage,
                     const std::array<Real, NumWindings> &ratedPower,
                     const std::array<Real, NumWindings> &ratioAbs,
                     const std::array<Real, NumWindings> &ratioPhase,
                     const std::array<Real, NumWindings> &resistance,
                     const std::array<Real, NumWindings> &inductance) {
    mNominalVoltage = nomVoltage;
    resolveReferenceWinding();

    **mRatedPowerP = ratedPower[0];
    **mRatedPowerS = ratedPower[1];
    **mRatedPowerT = ratedPower[2];

    **mRatioP = std::polar<Real>(ratioAbs[0], ratioPhase[0]);
    **mRatioS = std::polar<Real>(ratioAbs[1], ratioPhase[1]);
    **mRatioT = std::polar<Real>(ratioAbs[2], ratioPhase[2]);

    **mResistanceP = resistance[0];
    **mResistanceS = resistance[1];
    **mResistanceT = resistance[2];

    **mInductanceP = inductance[0];
    **mInductanceS = inductance[1];
    **mInductanceT = inductance[2];
  }
};
} // namespace Ph1
} // namespace Base
} // namespace CPS
