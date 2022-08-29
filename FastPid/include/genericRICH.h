/// \author R+Preghenella
/// \email  preghenella@bo.infn.it
/// \date   March 2020

#ifndef __GENERICRICH_H__
#define __GENERICRICH_H__

#include "genericDetector.h"

class genericRICH : public genericDetector {
 public:
  genericRICH() = default;
  virtual ~genericRICH() = default;

  /** setters **/
  void setIndex(double val) { mIndex = val; };
  void setEfficiency(double val) { mEfficiency = val; };
  void setMinPhotons(double val) { mMinPhotons = val; };
  void setThresholdMode(bool val) { mThresholdMode = val; };

  void setChromaticSigma(int n, double *valx, double *valy) {
    if (mChromaticSigma) delete mChromaticSigma;
    mChromaticSigma = new TGraph(n, valx, valy);
  }
  void setPositionSigma(int n, double *valx, double *valy) {
    if (mPositionSigma) delete mPositionSigma;
    mPositionSigma = new TGraph(n, valx, valy);
  }
  void setEmissionSigma(int n, double *valx, double *valy) {
    if (mEmissionSigma) delete mEmissionSigma;
    mEmissionSigma = new TGraph(n, valx, valy);
  }
  void setFieldSigma(int n, double *valx, double *valy) {
    if (mFieldSigma) delete mFieldSigma;
    mFieldSigma = new TGraph(n, valx, valy);
  }
  void setTrackingSigma(int n, double *valx, double *valy) {
    if (mTrackingSigma) delete mTrackingSigma;
    mTrackingSigma = new TGraph(n, valx, valy);
  }

  /** methods to override **/
  double numSigma(double eta, double p, PID::type PID);
  double maxP(double eta, double nsigma, PID::type PID);
  double minP(double eta, double nsigma, PID::type PID);

  double cherenkovAngle(double p, double m) const {
    double a = sqrt(m * m + p * p) / (mIndex * p);
    double b = acos(a);
    return b;
  };
  double cherenkovThreshold(double m) const {
    return m / sqrt(mIndex * mIndex - 1.);
  };
  double numberOfPhotons(double angle) const {
    return 490. * sin(angle) * sin(angle) * mLength;
  };
  double numberOfDetectedPhotons(double angle) const {
    return numberOfPhotons(angle) * mEfficiency;
  };
  double cherenkovAngleSigma(double eta, double p, double m) const;

 protected:
  // RICH parameters
  double mIndex = 1.0014;     // refractive index
  double mEfficiency = 0.25;  // overall photon detection efficiency
  double mMinPhotons = 3.;    // minimum number of detected photons

  // contributions to resolution
  TGraph *mChromaticSigma =
      nullptr;  // chromatic resolution vs. polar angle [rad]
  TGraph *mPositionSigma =
      nullptr;  // position resolution vs. polar angle [rad]
  TGraph *mEmissionSigma =
      nullptr;                    // emission resolution vs. polar angle [rad]
  TGraph *mFieldSigma = nullptr;  // field resolution vs. polar angle [rad]
  TGraph *mTrackingSigma =
      nullptr;  // tracking resolution vs. polar angle [rad]

  // threshold mode
  bool mThresholdMode = true;
};

#endif /* __GENERICRICH_H__ */