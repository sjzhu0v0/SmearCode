#include "tofSmear.h"
#include "math.h"
#include "iostream"

using namespace std;

tofSmear::tofSmear(double r, double eL, double eH, double sT,double sT0,typeDet dettype)
{
    detType = dettype;
    radius  = r;
    etaLow  = eL;
    etaHigh = eH;
    sigmaT  = sT;
    sigmaT0 = sT0;

    char nameit[500];
    sprintf(nameit,"TOF barrel R=%d dT=%d",int(radius),int(sigmaT));
    myName = nameit;

    //  Initialize a few constants:  Wikipedia...should be redone in a better manner...
    c       = 0.0299792458;  // cm/picosecond
    mPion   = 0.13957018;    //GeV/c^2
    mKaon   = 0.493677;      //GeV/c^2
    mElec   = 5.11e-4;     //GeV/C^2
    mProton = 0.93827208816; //GeV/c^2
  }


double tofSmear::tof(double L, double p, double m)
{
  double p2_m2 = p*p/(m*m);
  double beta  = sqrt(p2_m2/(1+p2_m2));
  return (L/(c*beta));
}

double tofSmear::GetTof(double eta,double p,double m)
{ 
  double theta = 2.0*atan( exp(-eta) );
  double L = radius/sin(theta);
  if(detType==endcap)
  {
    L = abs(radius/cos(theta));
  }
  smearedTof = tof(L,p,m)+RndGen.Gaus(0,sigmaT)+RndGen.Gaus(0,sigmaT0);
  recip_beta = smearedTof/L*c;
  return smearedTof;
}

double tofSmear::GetRecipBeta()
{
  return recip_beta;
}

double tofSmear::get_1_beta(double p,double m,double eta)
{
  double theta = 2.0*atan( exp(-eta) );
  double L = radius/sin(theta);
  if(detType==endcap)
  {
    L = abs(radius/cos(theta));
  }
  return tof(L,p,m)/L*c;
}

double tofSmear::get_1_beta_sigma(double eta)
{
  double theta = 2.0*atan( exp(-eta) );
  double L = radius/sin(theta);
  if(detType==endcap)
  {
    L = abs(radius/cos(theta));
  }
  return sqrt(sigmaT0*sigmaT0+sigmaT*sigmaT)/L*c;
}

void tofSmear::description() 
{
  //  Here one should describe what this detector is and what assumptions are 
  //  made in calculating the performance specifications.

  cout << "My name is \"" << myName << "\" and I am described as follows:" <<endl;
  cout << "    I am a Time-of-Flight barrel." <<endl;
  cout << "    I extend from eta =  " << etaLow << " until eta = " << etaHigh << " ." <<endl;
  cout << "    I am located at radius R= " << radius << " cm." <<endl;
  cout << "    I have a time resolution of " << sigmaT << " picoseconds." <<endl;
  cout << "    My calculations have assumed perfect momentum resolution and pointing." <<endl;
  cout << "    My calculations have assumed purely Gaussian time response." <<endl;
  cout << endl;

}
