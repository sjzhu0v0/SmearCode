#include "string"
#include "TRandom3.h"

#ifndef __TOFSMEAR_H__
#define __TOFSMEAR_H__
	
//
//  Hello tofBarrel Fans:
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
using namespace std;
	
class tofSmear
{
public:
  enum typeDet {endcap,barrel};
  tofSmear(double r, double eL, double eH, double sT,double sT0,typeDet dettype);
  
 //tofSmear(double radius=100, double etaLow=-1.0, double etaHigh=1.0, double sigmaT=10,double signmaT0=10,PID::typeDet detType=barrel); 
  virtual ~tofSmear() {}
	
  bool   valid   (double eta) {return (eta>etaLow && eta<etaHigh);}
  double GetTof  (double eta,double p,double m);
  double GetRecipBeta();
  string name    () {return myName;}
  void   description ();

  //addition function
  double get_1_beta(double p,double m,double eta);
  double get_1_beta_sigma(double eta);
		
protected:
  std::string myName;

  // utility function
  double tof(double L, double p, double m);
  // TOF wall parameters
  TRandom3 RndGen;
  typeDet detType;
  double radius;   // cm
  double etaLow;
  double etaHigh;
  double sigmaT;   // picosecond
  double sigmaT0;
  double smearedTof;
  double recip_beta;
  // Physical constants (should come from elsewhere!)
  double mPion;    // GeV/c^2
  double mKaon;    // GeV/c^2
  double mProton;  // GeV/c^2
  double mElec;    // GeV/c^2
  double c;        // cm/picosecond;
};
	
#endif /* __PID_H__ */
