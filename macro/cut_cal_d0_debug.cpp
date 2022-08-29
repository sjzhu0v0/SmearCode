#include "SmearEvent.h"

void cut_cal_d0(string name_input, int pid_system, string name_output) {
  // pid_system: 0 ECCE_PID, 1 ECCE_TOF_EMC, 2 ECCE_CHERENKOV_EMC, 3
  // EICC_ACCEPTANCE 4 perfect_pid 5 no pid
  SmearEvent smear_event(true, name_input);
  SmearEvent::CutStr_D0 cutstr_d0;

  TFile* f_cut = new TFile(name_output.c_str(), "recreate");

  TTree* t_cut = new TTree("Cut", "Cut");
  t_cut->Branch("inv_mass", &cutstr_d0.inv_mass, "inv_mass/D");
  t_cut->Branch("dca_kaon", &cutstr_d0.dca_kaon, "dca_kaon/D");
  t_cut->Branch("dca_pion", &cutstr_d0.dca_pion, "dca_pion/D");
  t_cut->Branch("cos_theta", &cutstr_d0.cos_theta, "cos_theta/D");
  t_cut->Branch("decay_length_D0", &cutstr_d0.decay_length_D0, "decay_length_D0/D");
  t_cut->Branch("dca_D02pv", &cutstr_d0.dca_D02pv, "dca_D02pv/D");
  t_cut->Branch("pT_D0", &cutstr_d0.pT_D0, "pT_D0/D");
  t_cut->Branch("rapidity_D0", &cutstr_d0.rapidity_D0, "rapidity_D0/D");
  t_cut->Branch("p3_mc", cutstr_d0.p3_mc, "p3_mc[3]/D");
  t_cut->Branch("p3_elec", cutstr_d0.p3_elec, "p3_elec[3]/D");

  bool right_signs;
  for (long iEvent = 0; iEvent < smear_event.GetEntries(); iEvent++) {
    smear_event.GetEntry(iEvent);
    int temp_var2 = smear_event.num_final_particles;
    for (int iPion = 0; iPion < smear_event.num_final_particles; iPion++) {
      for (int iKaon = 0; iKaon < smear_event.num_final_particles; iKaon++) {
        if (iPion == iKaon) continue;
        if (pid_system == 0) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::ECCE_PID,3.,iPion,iKaon);
        } else if (pid_system == 1) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::ECCE_TOF_EMC,3.,iPion,iKaon);
        } else if (pid_system == 2) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::ECCE_CHERENKOV_EMC,3.,iPion,iKaon);
        } else if (pid_system == 3) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::EICC_ACCEPTANCE,3.,iPion,iKaon);
        } else if (pid_system == 4) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::PERFECT_PID,3.,iPion,iKaon);
        } else if (pid_system == 5) {
          right_signs = smear_event.bool4D0Cut(SmearEvent::NO_PID,3.,iPion,iKaon);
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

int main(int argc, char* argv[]) {
  // const string name_input = argv[2];
  // cut_cal_d0(name_input, atoi(argv[2]), argv[3]);
  // cut_cal_d0("event_smeared_Det_v1.root", 0, "cut_event_smeared0_Det_v1.root");
  // cut_cal_d0("event_smeared_Det_v1.root", 1, "cut_event_smeared1_Det_v1.root");
  // cut_cal_d0("event_smeared_Det_v1.root", 2, "cut_event_smeared2_Det_v1.root");
  cut_cal_d0("event_smeared_Det_v1.root", 3, "cut_event_smeared3_Det_v1.root");
  // cut_cal_d0("event_smeared_Det_v1.root", 4, "cut_event_smeared4_Det_v1.root");
  // cut_cal_d0("event_smeared_Det_v1.root", 5, "cut_event_smeared5_Det_v1.root");
  return 0;
}
