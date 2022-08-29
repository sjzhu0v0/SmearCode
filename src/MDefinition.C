#include "MDefinition.h"

double min_tof_eta_region[3] = {-3.7, -1.7, 1.31};
double max_tof_eta_region[3] = {-1.7, 1.31, 4.21};
double min_cherenkov_eta_region[3] = {-3.5, -1.42, 1.42};
double max_cherenkov_eta_region[3] = {-1.42, 1.42, 3.9};
double min_emc_eta_region[3] = {-3.76, -1.71, 1.31};
double max_emc_eta_region[3] = {-1.81, 1.27, 4.21};
double min_momentum_applying_emc[3] = {4., 1., 6.};
// double min_momentum_applying_emc[3] = {4, 1., 5.8};

double p_track_becal[9] = {0.848343, 1.2452,  1.74209, 2.4879, 3.4807,
                           4.97083,  7.97601, 12.452,  17.4209};
double pi_rej_becal[9] = {194.353, 555.861, 420.883, 440.855, 216.556,
                          147.158, 83.0736, 57.331,  44.085};
double Nsigma_becal[9] = {2.7978,  3.12155, 3.03871, 3.05265, 2.83256,
                          2.70671, 2.51104, 2.37724, 2.27872};
double p_track_femc[14] = {0.842549, 1.2452,  1.74209, 2.4879,  3.4807,
                           4.97083,  7.97601, 12.3669, 17.4209, 22.2953,
                           27.3842,  32.2796, 37.0214, 44.5461};
double pi_rej_femc[14] = {16.6524, 62.901,  140.492, 203.575, 203.575,
                          124.154, 108.033, 92.5642, 70.087,  57.331,
                          49.1219, 44.772,  43.4094, 36.6234};
double Nsigma_femc[14] = {1.88042, 2.41125, 2.69128, 2.81274, 2.81274,
                          2.64978, 2.60244, 2.549,   2.45045, 2.37724,
                          2.31969, 2.28462, 2.27283, 2.20713};
double p_track_eemc[14] = {0.842549, 1.2452,  1.74209, 2.4879,  3.4807,
                           4.97083,  7.92154, 12.452,  17.4209, 22.2953,
                           27.3842,  32.2796, 37.2759, 44.5461};
double pi_rej_eemc[14] = {194.353, 582.238, 883.708, 793.102, 414.429,
                          290.461, 144.902, 122.25,  101.557, 78.094,
                          67.954,  62.901,  53.8944, 47.6269};
double Nsigma_eemc[14] = {2.7978,  3.13518, 3.25558, 3.22474, 3.03405,
                          2.92516, 2.70157, 2.64455, 2.58117, 2.48914,
                          2.4393,  2.41125, 2.35436, 2.30805};

TLorentzVector p4BeamHadron(0., 0., 25., 25.01760082);
TLorentzVector p4BeamElectron(0., 0., -5., 5.);

EicC_Mvd_DP mvd_dp;
dualRICH_C2F6 Mdual_rich_c2f6;
dualRICH_aerogel Mdual_rich_aerogel;
DrcPidFast MdrcPidFast;
mRICH MmRICH;
TOF MTof;
TVector3 nz(0., 0., 1.);

TRandom3 MRndgen_tof;
TRandom3 MRndgen_cherenkov;
TRandom3 MRndgen_emc;
TRandom3 MRndgen_eff;
TRandom3 MRndgen;

TRandom3 rndgenXY, rndgenZ, rndgenP, rndgenPVX, rndgenPVY, rndgenPVZ;

// EMC response map
TH2F *map_elec_eop_eff;
TH2F *map_pion_eop_eff;

TGraph *gr_emc_log_pi_rej[3]; // eemc becal femc
TGraph *gr_emc_Nsigma_pi_rej[3];

// TH2F *map_elec_matching_efficiency;
// TH2F *map_pion_matching_efficiency;

void ReadMap(string ver_det) {
  mvd_dp.Init_para(ver_det);
  // TFile *file_emcresponse = new TFile("~/Share/FastPid/EmcPid.root");
  // map_elec_eop_eff = (TH2F *)file_emcresponse->Get("h2_elec_eop_eff");
  // map_pion_eop_eff = (TH2F *)file_emcresponse->Get("h2_pion_eop_eff");

  // matching efficiency isnt used.
  // map_elec_matching_efficiency = new TH2F("h2_elec_matching_eff",
  // "h2_elec_matching_eff", 1, 0, 40, 1, -10, 10); map_pion_matching_efficiency
  // = new TH2F("h2_pion_matching_eff", "h2_pion_matching_eff", 1, 0, 40, 1,
  // -10, 10); map_elec_matching_efficiency->Fill(0.,0.);
  // map_pion_matching_efficiency->Fill(0.,0.);

  // matching efficiency is used.
  // map_elec_matching_efficiency = (TH2F
  // *)file_emcresponse->Get("h2_elec_matching_eff");
  // map_pion_matching_efficiency = (TH2F
  // *)file_emcresponse->Get("h2_pion_matching_eff");

  double log_p_track_becal[9];
  double log_pi_rej_becal[9];
  double log_p_track_femc[14];
  double log_pi_rej_femc[14];
  double log_p_track_eemc[14];
  double log_pi_rej_eemc[14];
  for (int i = 0; i < 9; i++) {
    log_p_track_becal[i] = log(p_track_becal[i]);
    log_pi_rej_becal[i] = log(pi_rej_becal[i]);
  }
  for (int i = 0; i < 14; i++) {
    log_p_track_femc[i] = log(p_track_femc[i]);
    log_pi_rej_femc[i] = log(pi_rej_femc[i]);
    log_p_track_eemc[i] = log(p_track_eemc[i]);
    log_pi_rej_eemc[i] = log(pi_rej_eemc[i]);
  }
  gr_emc_log_pi_rej[0] = new TGraph(14, log_p_track_eemc, log_pi_rej_eemc);
  gr_emc_log_pi_rej[1] = new TGraph(9, log_p_track_becal, log_pi_rej_becal);
  gr_emc_log_pi_rej[2] = new TGraph(14, log_p_track_femc, log_pi_rej_femc);

  gr_emc_Nsigma_pi_rej[0] = new TGraph(14, log_p_track_eemc, Nsigma_eemc);
  gr_emc_Nsigma_pi_rej[1] = new TGraph(9, log_p_track_becal, Nsigma_becal);
  gr_emc_Nsigma_pi_rej[2] = new TGraph(14, log_p_track_femc, Nsigma_femc);
  MdrcPidFast.ReadMap("~/Share/FastPid/ctr_map_p1_0.95.root");
}

void p_smear(TVector3 *p, TRandom3 *rndgenP, int pid) {
  double sigma_p = mvd_dp.GetResP(p->Mag(), p->PseudoRapidity(), pid);
  *p = p->Unit() * p->Mag() * (1 + rndgenP->Gaus(0, sigma_p));
}

void dca_smear(double px, double py, double pz, double vx, double vy, double vz,
               int pid, double *dcaz, double *dcaxy, TRandom3 *rndgenXY,
               TRandom3 *rndgenZ) {
  double sigma_dcaz = 0;
  double sigma_dcaxy = 0;
  TVector3 momentum(px, py, pz);
  TVector3 coor(vx, vy, vz);
  double sign = vx * py - vy * px > 0 ? +1. : -1.;
  sigma_dcaz =
      mvd_dp.GetResDCAz(momentum.Perp(), momentum.PseudoRapidity(), pid) /
      10000.;
  sigma_dcaxy =
      mvd_dp.GetResDCArp(momentum.Perp(), momentum.PseudoRapidity(), pid) /
      10000.;
  *dcaz = (coor -
           (coor.Dot(momentum.Cross(nz).Unit())) * (momentum.Cross(nz).Unit()))
              .Z() +
          rndgenZ->Gaus(0, sigma_dcaz);
  *dcaxy = sign * ((coor.Dot(momentum.Cross(nz).Unit())) *
                   (momentum.Cross(nz).Unit()))
                      .Mag() +
           rndgenXY->Gaus(0, sigma_dcaxy);
}

double eta_to_theta(double eta) { return 2. * atan(exp(-1. * eta)); }

void coor_cal(TVector3 p, double dcaxy, double dcaz, TVector3 *coor) {
  *coor = dcaxy * (p.Cross(nz).Unit());
  coor->SetZ(dcaz);
}

double dca_cal(TVector3 p, TVector3 c, TVector3 *c_, TVector3 _c) {
  /*p,c and _c are inputs and c_(output) is coordinate of DCA from the track
      constructed with p and c to _c. _c is the original point by default.*/
  TVector3 a = (c - _c) - p.Unit() * (c - _c) * p.Unit();
  *c_ = _c + a;
  return a.Mag();
}

double trackdca_cal(TVector3 p1, TVector3 p2, TVector3 coor1, TVector3 coor2,
                    TVector3 *coor3) {
  /*p1 is the momentum of the first daughter.coor1 is the coordinate of the
     first daughter. p2 is ... . coor2 is ... . coor3 is the middle point of DCA
     of the tracks of the two daughters.*/
  TVector3 n = p1.Cross(p2);
  n = n.Unit();
  TVector3 delta = coor2 - coor1 + (coor2 - coor1) * n * n;
  double sol = ((delta * p1) * (p2 * p2) - (p1 * p2) * (p2 * delta)) /
               ((p1 * p1) * (p2 * p2) - (p1 * p2) * (p2 * p1));
  /*sol is the solution of a_1 of
   * equation{p1*delta=p1*(a_1*p1+a_2*p2);p2*delta=p2*(a_1*p1+a_2*p2)}*/
  *coor3 = coor1 + sol * p1 + 0.5 * delta;
  return abs(n * (coor1 - coor2));
}