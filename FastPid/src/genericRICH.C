#include "genericRICH.h"

double genericRICH::cherenkovAngleSigma(double eta, double p, double m) const {
  auto theta = 2. * atan(exp(-eta)) * 57.295780;
  auto chromatic = mChromaticSigma ? mChromaticSigma->Eval(theta) : 0.;
  auto position = mPositionSigma ? mPositionSigma->Eval(theta) : 0.;
  auto emission = mEmissionSigma ? mEmissionSigma->Eval(theta) : 0.;
  auto field = mFieldSigma ? mFieldSigma->Eval(theta) : 0.;
  auto tracking = mTrackingSigma ? mTrackingSigma->Eval(theta) : 0.;

  // contributions that scale with number of detected photons
  auto ndet = numberOfDetectedPhotons(cherenkovAngle(p, m));
  auto sigma1 = sqrt(chromatic * chromatic + position * position +
                     emission * emission + field * field + tracking * tracking);
  // contributions that do not
  auto sigma2 = 0.;
  //
  return sqrt(sigma1 * sigma1 / ndet + sigma2 * sigma2);
};

double genericRICH::numSigma(double eta, double p, PID::type PID) {
  double mass1, mass2;
  switch (PID) {
    case e_pi:
      mass1 = mMassElectron;
      mass2 = mMassPion;
      break;
    case pi_k:
      mass1 = mMassPion;
      mass2 = mMassKaon;
      break;
    case k_p:
      mass1 = mMassKaon;
      mass2 = mMassProton;
      break;
  }

  double thr1 = cherenkovThreshold(mass1);
  double thr2 = cherenkovThreshold(mass2);
  // cout << thr1 << " thres " << thr2 << " p " << p<<endl;

  /** both particles are above threshold **/
  if (p > thr1 && p > thr2) {
    return (cherenkovAngle(p, mass1) - cherenkovAngle(p, mass2)) /
           cherenkovAngleSigma(eta, p, mass1) / sqrt(2.);
  }

  /** lightest particle above threshold **/
  if (mThresholdMode && p > thr1) {
    return (cherenkovAngle(thr2 + 0.001, mass1) -
            cherenkovAngle(thr2 + 0.001, mass2)) /
           cherenkovAngleSigma(eta, thr2 + 0.001, mass1);
  }

  /** none above threshold **/
  return 0.;
}

double genericRICH::maxP(double eta, double nsigma, PID::type PID) {
  double mass1, mass2;
  switch (PID) {
    case e_pi:
      mass1 = mMassElectron;
      mass2 = mMassPion;
      break;
    case pi_k:
      mass1 = mMassPion;
      mass2 = mMassKaon;
      break;
    case k_p:
      mass1 = mMassKaon;
      mass2 = mMassProton;
      break;
  }

  /** let's do it numerically, starting from the threshold **/
  double p = minP(eta, nsigma, PID) + 0.001;
  while (numSigma(eta, p, PID) > nsigma) p += 0.001;
  return p;
}

double genericRICH::minP(double eta, double nsigma, PID::type PID) {
  double mass;
  switch (PID) {
    case e_pi:
      mass = mThresholdMode ? mMassElectron : mMassPion;
      break;
    case pi_k:
      mass = mThresholdMode ? mMassPion : mMassKaon;
      break;
    case k_p:
      mass = mThresholdMode ? mMassKaon : mMassProton;
      break;
  }

  /** let's do it numerically, starting from the threshold **/
  double p = cherenkovThreshold(mass);

  while (numberOfPhotons(cherenkovAngle(p, mass)) * mEfficiency < mMinPhotons)
    p += 0.001;
  return std::max(p, pMin(eta));
}