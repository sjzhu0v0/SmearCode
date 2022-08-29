// Fast reconstruction for EIC Barrel DIRC
// original author: r.dzhygadlo at gsi.de

#ifndef DrcPidFast_h
#define DrcPidFast_h 1

#include "TFile.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TF1.h"

#include <iostream>

// probability - normalized to 1 probability for e,mu,pi,k,p
// sigma - deviation of the determined Cherenkov angle from expected in terms of Cherenkov track
// resolution cangle - Cherenkov angle cctr -  combined Cherenkov track resolution
struct DrcPidInfo
{
  double probability[5];
  double sigma[5];
  double cangle;
  double cctr;
};

class DrcPidFast
{

public:
  // barid = 0 for 17 mm thickness of the radiator
  // barid = 1 for 10 mm thickness of the radiator
  DrcPidFast(int barid = 0);
  ~DrcPidFast() {}

  // read Cherenkov track resolution map from a file
  void ReadMap(TString name);

  // pdg - Particle Data Group code of the particle
  // mom - 3-momentum of the particle [GeV/c]
  // track_err - error assosiated with track direction [mrad]
  DrcPidInfo GetInfo(int pdg, TVector3 mom, double track_err = 0);

  // p - momentum of the particle [GeV/c]
  // theta - polar angle of the particle [deg]
  DrcPidInfo GetInfo(int pdg, double p, double theta, double track_err = 0);
  TH2F *GetTrrMap() { return fTrrMap; }
  //additional function

  double GetCctr(int pid, double theta, double p, double track_err);
  double ePi_maxP(double NSigma, double theta, double track_err);
  double MomentumThreshold(int pid);
  double MomentumThreshold(double mass)
  {
    return mass / sqrt(1.46907 * 1.46907 - 1.);
  }

  double GetAng(int pid, double p);
  double GetAng(double mass, double p);

private:
  int get_pid(int pdg);
  TH2F *fTrrMap;
  double fMass[5];
  TF1 *fMs_mom, *fMs_thickness;
  double fMs_thickness_max;
};

#endif
