#include "genericDetector.h"

double genericDetector::etaMin() {
  switch (mType) {
    case kBarrel:
      return -log(tan(atan2(mRadius, -mLength) * 0.5));
    case kForward:
      return -log(tan(atan2(mRadiusOut, mPositionZ) * 0.5));
  }
  return 0.;
}

double genericDetector::etaMax() {
  switch (mType) {
    case kBarrel:
      return -log(tan(atan2(mRadius, mLength) * 0.5));
    case kForward:
      return -log(tan(atan2(mRadiusIn, mPositionZ) * 0.5));
  }
  return 0.;
}

double genericDetector::ptMin() {
  switch (mType) {
    case kBarrel:
      return 0.003 * mMagneticField * 0.5 * mRadius;
    case kForward:
      return 0.003 * mMagneticField * 0.5 * mRadiusIn;
  }
  return 0.;
}

double genericDetector::trackLength(double eta) {
  auto theta = 2. * atan(exp(-eta));

  switch (mType) {
    case kBarrel:
      return mRadius / sin(theta);
    case kForward:
      return mPositionZ / cos(theta);
  }
  return 0.;
}