#ifndef __MRICH_H__
#define __MRICH_H__
	
//
//  Hello mRICH Fans:
//
//  This is an example class that inherits from the PID base class.
//  We shall implement all the necessary functions as required by the
//  base class.
//
//  The UNIQUE identifiers for this class are radius, eta extent, and 
//  time resolution. 
//
//  Note that in keeping with any well written base class, we shall
//  give default values to all the constructor arguments.
//
//  This routine assumes units of cm for distances and picoSeconds for time.
//

#include "PID.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include "TFile.h"
#include "string"

using namespace std;

	
class mRICH: public PID
{
public:
  mRICH(double trackResolution=0.5, double timePrecision=1.0, double pixS=0.5, double p=5.0);
  virtual ~mRICH() {}
	
  //bool   valid   (double eta, double p) {return (((eta>-8&&eta<-1) || (eta>1&&eta<8)) && (p>pLow && p<pHigh));}
  bool   valid   (double eta) {return (eta>-8&&eta<8);}
  double maxP    (double eta, double numSigma, PID::type PID);
  double minP    (double eta, double numSigma, PID::type PID) {return 0;}
  string name    () {return myName;}
  void   description ();
  void ReadMap(TString name);
  double getAng(double mass);
  double getdAng(double mass);
  double getNgamma(double t, double mass);
  double getdAngTrk(double mass);
  double T_Aer(double t, double lam);
  double T_QE(double lam);

  double SetLensFocalLength(double f);
  double SetAerogelThickness(double t);
  double GetmRICHParameters(); 

  //additional function
  double maxP_epi_separation(double Sigma);
  double MomentumThreshold(double mass);
  double getAng(double mass,double _mom);
  double getdAng(double mass,double _mom);
  double getdAngTrk(double mass,double _mom);
  double getNgamma(double t, double mass,double _mom);
 protected:
  std::string myName;

  double Nsigma;
  double Nsigma_epi;

  // Physical constants (should come from elsewhere!)
  double mPion;    // GeV/c^2
  double mKaon;    // GeV/c^2
  double mProton;  // GeV/c^2
  double mom;      // GeV/c
  double mElectron;
  double c;        // cm/picosecond;
  double n;
  double f;        //mm
  double a;        //mm
  double N_gam;
  double pi;
  double alpha;
  double L;
  double th0;

  double fTrackResolution;
  double fTimePrecision;

  double pLow;
  double pHigh;

  TH2F *fpixMap;

};
	
#endif /* __PID_H__ */
