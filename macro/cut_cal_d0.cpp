#include "SmearEvent.h"

void cut_cal_d0(string name_input, int pid_system, string name_output,
                string name_mcinput, string is_signal_cal) {
  // pid_system: 0 ECCE_PID, 1 ECCE_TOF_EMC, 2 ECCE_CHERENKOV_EMC, 3
  // EICC_ACCEPTANCE 4 perfect_pid 5 no pid
  bool is_signal = false;
  if (is_signal_cal == "Signal") {
    is_signal = true;
  }

  SmearEvent::CutStr_D0 cutstr_d0;
  TFile *f_mc = new TFile(name_mcinput.c_str());
  SmearEvent smear_event(true, name_input);
  TFile *f_cut = new TFile(name_output.c_str(), "recreate");
  TTree *t_event = (TTree *)f_mc->Get("EICTree");
  cout << t_event << endl;
  smear_event.McTreeInitializing(t_event);

  TTree *t_cut = new TTree("Cut", "Cut");
  t_cut->Branch("inv_mass", &cutstr_d0.inv_mass, "inv_mass/D");
  t_cut->Branch("dca_kaon", &cutstr_d0.dca_kaon, "dca_kaon/D");
  t_cut->Branch("dca_pion", &cutstr_d0.dca_pion, "dca_pion/D");
  t_cut->Branch("cos_theta", &cutstr_d0.cos_theta, "cos_theta/D");
  t_cut->Branch("decay_length_D0", &cutstr_d0.decay_length_D0,
                "decay_length_D0/D");
  t_cut->Branch("dca_12", &cutstr_d0.dca_12, "dca_12/D");
  t_cut->Branch("pT_D0", &cutstr_d0.pT_D0, "pT_D0/D");
  t_cut->Branch("rapidity_D0", &cutstr_d0.rapidity_D0, "rapidity_D0/D");
  t_cut->Branch("p3_mc", cutstr_d0.p3_mc, "p3_mc[3]/D");
  t_cut->Branch("p3_elec", cutstr_d0.p3_elec, "p3_elec[3]/D");
  t_cut->Branch("D0_Sign", &cutstr_d0.D0_Sign, "D0_Sign/I");
  t_cut->Branch("num_final_particles", &cutstr_d0.num_final_particles,
                "num_final_particles/I");
  t_cut->Branch("num_mc_final_charged_particles",
                &cutstr_d0.num_mc_final_charged_particles,
                "num_mc_final_charged_particles/I");

  bool right_signs;
  for (long iEvent = 0; iEvent < smear_event.GetEntries(); iEvent++) {
    smear_event.GetEntry(iEvent);

    if (is_signal_cal != "all") {
      if (is_signal ? !smear_event.isSignal(421) : smear_event.isSignal(421))
        continue;
    }

    int temp_var2 = smear_event.num_final_particles;
    for (int iPion = 0; iPion < smear_event.num_final_particles; iPion++) {
      for (int iKaon = 0; iKaon < smear_event.num_final_particles; iKaon++) {
        if (iPion == iKaon) continue;
        if (pid_system == 0) {
          right_signs =
              smear_event.bool4D0Cut(SmearEvent::ECCE_PID, 2., iPion, iKaon);
        } else if (pid_system == 1) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::ECCE_TOF_EMC, 2.,
                                               iPion, iKaon);
        } else if (pid_system == 2) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::ECCE_CHERENKOV_EMC,
                                               2., iPion, iKaon);
        } else if (pid_system == 3) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::EICC_ACCEPTANCE, 2.,
                                               iPion, iKaon);
        } else if (pid_system == 4) {
          right_signs =
              smear_event.bool4D0Cut(SmearEvent::PERFECT_PID, 2., iPion, iKaon);
        } else if (pid_system == 5) {
          right_signs =
              smear_event.bool4D0Cut(SmearEvent::NO_PID, 2., iPion, iKaon);
        }
        if (right_signs) {
          cutstr_d0 = smear_event.CalCut_D0(iPion, iKaon);
          t_cut->Fill();
        }
      }
    }
  }
  f_cut->cd();
  t_cut->Write();
  f_cut->Close();
}

int main(int argc, char *argv[]) {
  const string name_input = argv[1];
  cut_cal_d0(name_input, atoi(argv[2]), argv[3], argv[4], argv[5]);
  // cut_cal_d0("~/HIUser_sjzhu/work3/code/Smear_Code/bin/d0_decay_test.root"
  // ,3, "cut_d0_test.root"
  // ,"/ustcfs/HICUser/ephy/work/sjzhu/Singularity-master/cvmfs/work/pythia6_generator/3.5_20/generator1/mc.root");
  return 0;
}
