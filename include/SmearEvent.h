#ifndef __SmearEvent_h__
#define __SmearEvent_h__

#include "PidMap.h"
#include "eicsmear/erhic/EventMC.h"
#include "eicsmear/functions.h"

using namespace std;

class SmearEvent {
public:
  enum USED_PID_DET {
    ECCE_PID,
    ECCE_TOF_EMC,
    ECCE_CHERENKOV_EMC,
    EICC_ACCEPTANCE,
    PERFECT_PID,
    NO_PID
  };

  typedef struct {
    TVector3 p_mc;
    bool isDetected = false;
    TVector3 p_smeared;
    double smeared_dcaxy = 0;
    double smeared_dcaz = 0;
    int pid;
    int expected_pid[4]; // 0 ecce_pid 1 ecce_tof_emc 2 ecce_cherenkov_emc 3
                         // eicc_acceptance
    double n_sigma_pid_ecce[n_particle_type]
                           [n_det_type]; // the second index is pid system, the
                                         // first index is pid:0 electron 1 pion
                                         // 2 kaon 3 proton
    double n_sigma_pid_eicc[n_particle_type]
                           [n_det_type]; // the second index is pid system, the
                                         // first index is pid:0 electron 1 pion
                                         // 2 kaon 3 proton
  } SmearParticleStr;

  typedef struct {
    double dca_proton = 0;
    double dca_kaon = 0;
    double dca_pion = 0;
    double cos_theta = 0;
    double decay_length_lambda_c = 0;
    double maxDca_daughters = 0;
    double dca_lambda_c2pv = 0;
    double inv_mass = 0;
    double pT_lambda_c = 0;
    double rapidity_lambda_c = 0;
    int Lambda_c_Sign = 0;
    double p3_mc[3] = {0, 0, 0};
    double p3_elec[3] = {0, 0, 0};
    int num_final_particles = 0;
    int num_mc_final_charged_particles = 0;
  } CutStr_lambda_c;

  typedef struct {
    double dca_kaon = 0;
    double dca_pion = 0;
    double cos_theta = 0;
    double dca_12 = 0;
    double decay_length_D0 = 0;
    double inv_mass = 0;
    double pT_D0 = 0;
    double rapidity_D0 = 0;
    int D0_Sign = 0;
    double p3_mc[3] = {0, 0, 0};
    double p3_elec[3] = {0, 0, 0};
    int num_final_particles = 0;
    int num_mc_final_charged_particles = 0;
  } CutStr_D0;

  typedef struct {
    double dca_electron = 0;
    double dca_positron = 0;
    double cos_theta = 0;
    double decay_length_Jpsi = 0;
    double dca_12 = 0;
    double inv_mass = 0;
    double pT_Jpsi = 0;
    double rapidity_Jpsi = 0;
    double p3_mc[3] = {0, 0, 0};
    double p3_elec[3] = {0, 0, 0};
    int num_final_particles = 0;
    int num_mc_final_charged_particles = 0;
  } CutStr_Jpsi;

  typedef struct {
    TVector3 p3_pion;
    TVector3 p3_kaon;
    TVector3 p3_proton;
    TVector3 c3_pion;
    TVector3 c3_kaon;
    TVector3 c3_proton;
  } McEventStr;

  int event_id = 0;
  int num_final_particles = 0;
  double true_px[40];
  double true_py[40];
  double true_pz[40];
  double smeared_px[40];
  double smeared_py[40];
  double smeared_pz[40];
  double smeared_dcaxy[40];
  double smeared_dcaz[40];
  double n_sigma_pid_ecce[40][n_particle_type][n_det_type] = {{{0}}};
  double n_sigma_pid_eicc[40][n_particle_type][n_det_type] = {{{0}}};
  int true_pid[40];

  double smeared_px_scattered_lepton = 0;
  double smeared_py_scattered_lepton = 0;
  double smeared_pz_scattered_lepton = 0;
  double true_px_scattered_lepton = 0;
  double true_py_scattered_lepton = 0;
  double true_pz_scattered_lepton = 0;

  double prime_vx = 0.;
  double prime_vy = 0.;
  double prime_vz = 0.;

  SmearEvent(bool isRead, string fileName);
  ~SmearEvent();

  TVector3 GetTrueP3(int iTrack) {
    return TVector3(true_px[iTrack], true_py[iTrack], true_pz[iTrack]);
  }
  TVector3 GetSmearedP3(int iTrack) {
    return TVector3(smeared_px[iTrack], smeared_py[iTrack], smeared_pz[iTrack]);
  }
  TVector3 GetPrimeVertex() {
    return TVector3(prime_vx / 10000., prime_vy / 10000., prime_vz / 10000.);
  };
  TVector3 GetPV3() {
    return TVector3(prime_vx / 10000., prime_vy / 10000., prime_vz / 10000.);
  } // micrometer 2 centimeter
  bool IsTriggered() { return is_triggered; }
  bool isScatteredElectronDetected() { return is_scattered_electron_detected; };
  void GetEntry(long iEvent);
  long GetEntries() { return tree_event->GetEntries(); }
  void Fill() { tree_event->Fill(); }
  void Save();
  void EventInput(erhic::EventMC *mcEvent,
                  int index); // void GetEventData(SmearEvent *event);
  void McEventInput(McEventStr mcEvent);

  void McTreeInitializing(TTree *tree);
  CutStr_lambda_c
  CalCut_Lambda_c(int no1, int no2,
                  int no3); // calculate cuts no1: pion, no2:kaon, no3:proton
  bool bool4D0Cut(USED_PID_DET used_pid_det, double n_sigma, int no1, int no2);
  bool bool4JpsiCut(USED_PID_DET used_pid_det, double n_sigma, int no1,
                    int no2); // no1 electron;no2 positron
  bool bool4LowMomentum(int no);
  bool bool4Lambda_cCut(USED_PID_DET used_pid_det, double n_sigma, int no1,
                        int no2, int no3);
  CutStr_D0 CalCut_D0(int no1, int no2);
  CutStr_Jpsi CalCut_Jpsi(int no1, int no2);

  SmearParticleStr SmearFunction(const erhic::ParticleMC *mcParticle);
  SmearParticleStr SmearFunction(TVector3 p3McParticle, TVector3 v3McParticle,
                                 int pid);

private:
  TFile *file_event;
  TTree *tree_event;
  TTree *tree_mcevent=nullptr;

protected:
  bool is_scattered_electron_detected = true;
  int int_is_scattered_electron_detected = 1; // 0 no 1 yes
  bool is_triggered = true;
  TVector3 c_pv;

  erhic::EventMC *mc_event = nullptr;
  vector<SmearParticleStr> detected_particles;
};

#endif
