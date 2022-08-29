#ifndef __MDEFINITION_H__
#define __MDEFINITION_H__

#include "EicC_Mvd_DP/EicC_Mvd_DP.h"
#include "Math/ProbFuncMathCore.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObject.h"
#include "TPad.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "iostream"
#include "vector"

#include "dualRICH.h"
#include "mRICH.h"
#include "TOF.h"
#include "DrcPidFast.h"

using namespace std;

// Detector Eta Region  0 backward 1 barrel 2 endcap
extern double min_tof_eta_region[];
extern double max_tof_eta_region[];
extern double min_cherenkov_eta_region[];
extern double max_cherenkov_eta_region[];
extern double min_emc_eta_region[];
extern double max_emc_eta_region[];
extern double min_momentum_applying_emc[];
extern double p_track_becal[];
extern double pi_rej_becal[];
extern double Nsigma_becal[];
extern double p_track_femc[];
extern double pi_rej_femc[];
extern double Nsigma_femc[];
extern double p_track_eemc[];
extern double pi_rej_eemc[];
extern double Nsigma_eemc[];

extern TLorentzVector p4BeamHadron;
extern TLorentzVector p4BeamElectron;

extern EicC_Mvd_DP mvd_dp;
extern dualRICH_C2F6 Mdual_rich_c2f6;
extern dualRICH_aerogel Mdual_rich_aerogel;
extern DrcPidFast MdrcPidFast;
extern mRICH MmRICH;
extern TOF MTof;

extern TVector3 nz;
extern TRandom3 MRndgen_tof;
extern TRandom3 MRndgen_cherenkov;
extern TRandom3 MRndgen_emc;
extern TRandom3 MRndgen_eff;
extern TRandom3 MRndgen;
extern TRandom3 rndgenXY, rndgenZ, rndgenP, rndgenPVX, rndgenPVY, rndgenPVZ;

extern TH2F *map_elec_eop_eff;
extern TH2F *map_pion_eop_eff;
extern TGraph *gr_emc_log_pi_rej[3];  // eemc becal femc
extern TGraph *gr_emc_Nsigma_pi_rej[3];

// TH2F *map_elec_matching_efficiency;
// TH2F *map_pion_matching_efficiency;

void ReadMap(string ver_det);

void p_smear(TVector3 *p, TRandom3 *rndgenP, int pid);

void dca_smear(double px, double py, double pz, double vx, double vy, double vz,
               int pid, double *dcaz, double *dcaxy, TRandom3 *rndgenXY,
               TRandom3 *rndgenZ);

double eta_to_theta(double eta);

void coor_cal(TVector3 p, double dcaxy, double dcaz, TVector3 *coor);

double dca_cal(TVector3 p, TVector3 c, TVector3 *c_, TVector3 _c);

double trackdca_cal(TVector3 p1, TVector3 p2, TVector3 coor1, TVector3 coor2,
                    TVector3 *coor3);
#endif