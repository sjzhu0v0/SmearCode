#include "MDefinition.h"

#ifndef __PIDMAP_H__
#define __PIDMAP_H__

const int n_particle_type = 4;
const int n_det_type = 3;

class PidMap {
public:
  enum PID { electron = 0, pion = 1, kaon = 2, proton = 3 };
  enum DET { TOF = 0, cherenkov = 1, EMC = 2 };
  enum PID_SYSTEM { ECCE_PID = 0, EICC_ACCEPTANCE = 1 };

  int PDGID[4] = {11, 211, 321, 2212};
  PidMap::PID all_particle_type[n_particle_type] = {
      PidMap::electron, PidMap::pion, PidMap::kaon, PidMap::proton};
  PidMap(PidMap::PID pid, double eta, double p);
  PidMap(int pid, double eta, double p);
  ~PidMap() {}

  void is_hit();
  int in_which_region(int iDet);
  void smear_particle_mcinfo();
  void cal_pidmap();

  double get_probability_to_particle(PidMap::PID toWhichParticle) {
    return bayes_probability_to_particle[toWhichParticle];
  }

  void GetPidMap(PID_SYSTEM pid_system, double *pidmap);

  int get_expected_pid_int(PID_SYSTEM pid_system);

  double probability_to_particle_tof(PidMap::PID toWhichParticle);
  double probability_to_particle_cherenkov(PidMap::PID toWhichParticle);
  double probability_to_particle_emc(PidMap::PID toWhichParticle);

  double sigma_to_particle_tof();

  double Nsigma_to_particle_tof(PidMap::PID toWhichParticle);
  double Nsigma_to_particle_cherenkov(PidMap::PID toWhichParticle);
  double Nsigma_to_particle_emc(PidMap::PID toWhichParticle);

  // cherenkov detector
  double probability_to_particle_dRICH_C2F6(PidMap::PID toWhichParticle);
  double probability_to_particle_dRICH_aerogel(PidMap::PID toWhichParticle);
  double probability_to_particle_hpDIRC(PidMap::PID toWhichParticle);
  double probability_to_particle_mRICH(PidMap::PID toWhichParticle);

  double sigma_to_particle_dRICH_C2F6(PidMap::PID toWhichParticle);
  double sigma_to_particle_dRICH_aerogel(PidMap::PID toWhichParticle);
  double sigma_to_particle_hpDIRC(PidMap::PID toWhichParticle);
  double sigma_to_particle_mRICH(PidMap::PID toWhichParticle);

  double Nsigma_to_particle_dRICH_C2F6(PidMap::PID toWhichParticle);
  double Nsigma_to_particle_dRICH_aerogel(PidMap::PID toWhichParticle);
  double Nsigma_to_particle_hpDIRC(PidMap::PID toWhichParticle);
  double Nsigma_to_particle_mRICH(PidMap::PID toWhichParticle);

  double momentum_threshold_dRICH(
      int type,
      PidMap::PID
          toWhichParticle); // forward //type1: threshold of c2f6 is used.
                            // type0: threshold of aerogel is used.
  double momentum_threshold_hpDIRC(PidMap::PID toWhichParticle); // barrel
  double momentum_threshold_mRICH(PidMap::PID toWhichParticle);  // backward

private:
  PidMap::PID pdg_id;
  bool isEmcUsed = false;
  double momentum;
  double pseduorapidity;

  // EicC PID acceptance
  double min_pid_eta_region[3] = {-3.5, -1, 1};
  double max_pid_eta_region[3] = {-1, 1, 3.5};
  double PID_acceptance[4][3] = {
      {30, 30, 20},
      {4, 6, 15},
      {4, 6, 15},
      {4, 6, 15}}; // [identification type][detector region]

  // random number
  double rnd_tof = MRndgen_tof.Gaus();
  double rnd_cherenkov = MRndgen_cherenkov.Gaus();
  double rnd_emc = MRndgen.Uniform(0, 1);
  double rnd_emcMatching = MRndgen.Uniform(0, 1);

  // electron  //pion
  double fMass[n_particle_type] = {5.11e-4, 0.13957, 0.493677, 0.938272};
  double prob_sum = 0;
  double bayes_probability_to_particle[n_particle_type] = {1., 1., 1., 1.};
  double Nsigma_to_particle[n_particle_type][n_det_type] = {{1.}};
  double TotalSigma_to_particle[n_particle_type] = {0.};
  // TOF  //cherenkov //EMC
  int det_hit[n_det_type] = {
      0}; // 0 not covered, 1 backward, 2 barrel, 3 forward, -1 backward but no
          // hit, -2 ..., -3 ...
};

#endif