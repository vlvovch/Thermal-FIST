/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/ParticleDecay.h"

#include <cmath>

using namespace std;

namespace thermalfist {

  double ParticleDecayChannel::ModifiedWidth(double m) const
  {
    if (m < mM0) return 0.;
    if (mM0 >= mPole)
      return mBratio * pow(1. - (mM0 / m)*(mM0 / m), mL + 1. / 2.);
    return mBratio * pow(1. - (mM0 / m)*(mM0 / m), mL + 1. / 2.) / pow(1. - (mM0 / mPole)*(mM0 / mPole), mL + 1. / 2.);
  }

  bool ParticleDecayChannel::operator==(const ParticleDecayChannel & rhs) const
  {
    bool ret = true;

    ret &= mBratio == rhs.mBratio;
    ret &= mDaughters == rhs.mDaughters;
    ret &= mM0 == rhs.mM0;
    ret &= mPole == rhs.mPole;
    ret &= mL == rhs.mL;
    ret &= mBratioVsM == rhs.mBratioVsM;
    ret &= mBratioAverage == rhs.mBratioAverage;
    ret &= mChannelName == rhs.mChannelName;

    return ret;
  }

} // namespace thermalfist