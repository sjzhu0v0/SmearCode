#include "SmearEvent.h"

void cut_cal_Jpsi(string name_input, int pid_system, string name_output,
                  string name_mcinput) {
  // pid_system: 0 ECCE_PID, 1 ECCE_TOF_EMC, 2 ECCE_CHERENKOV_EMC, 3
  // EICC_ACCEPTANCE 4 perfect_pid 5 no pid
  SmearEvent smear_event(true, name_input);
  SmearEvent::CutStr_Jpsi cutstr_Jpsi;

  TFile *f_cut = new TFile(name_output.c_str(), "recreate");
  TTree *t_event = (TTree *)(new TFile(name_mcinput.c_str()))->Get("EICTree");
  smear_event.McTreeInitializing(t_event);

  TTree *t_cut = new TTree("Cut", "Cut");
  t_cut->Branch("inv_mass", &cutstr_Jpsi.inv_mass, "inv_mass/D");
  t_cut->Branch("dca_electron", &cutstr_Jpsi.dca_electron, "dca_electron/D");
  t_cut->Branch("dca_positron", &cutstr_Jpsi.dca_positron, "dca_positron/D");
  t_cut->Branch("cos_theta", &cutstr_Jpsi.cos_theta, "cos_theta/D");
  t_cut->Branch("decay_length_Jpsi", &cutstr_Jpsi.decay_length_Jpsi,
                "decay_length_Jpsi/D");
  t_cut->Branch("dca_Jpsi2pv", &cutstr_Jpsi.dca_Jpsi2pv, "dca_Jpsi2pv/D");
  t_cut->Branch("pT_Jpsi", &cutstr_Jpsi.pT_Jpsi, "pT_Jpsi/D");
  t_cut->Branch("rapidity_Jpsi", &cutstr_Jpsi.rapidity_Jpsi, "rapidity_Jpsi/D");
  t_cut->Branch("p3_mc", cutstr_Jpsi.p3_mc, "p3_mc[3]/D");
  t_cut->Branch("p3_elec", cutstr_Jpsi.p3_elec, "p3_elec[3]/D");
  t_cut->Branch("num_final_particles", &cutstr_Jpsi.num_final_particles,
                "num_final_particles/I");
  t_cut->Branch("num_mc_final_charged_particles",
                &cutstr_Jpsi.num_mc_final_charged_particles,
                "num_mc_final_charged_particles/I");

  bool right_signs;
  for (long iEvent = 0; iEvent < smear_event.GetEntries(); iEvent++) {
    smear_event.GetEntry(iEvent);
    int temp_var2 = smear_event.num_final_particles;
    for (int iElec = 0; iElec < smear_event.num_final_particles; iElec++) {
      for (int iPosi = 0; iPosi < smear_event.num_final_particles; iPosi++) {
        if (iElec == iPosi)
          continue;
        if (pid_system == 0) {
          right_signs =
              smear_event.bool4JpsiCut(SmearEvent::ECCE_PID, 2., iElec, iPosi);
        } else if (pid_system == 1) {
          right_signs = smear_event.bool4JpsiCut(SmearEvent::ECCE_TOF_EMC, 2.,
                                                 iElec, iPosi);
        } else if (pid_system == 2) {
          right_signs = smear_event.bool4JpsiCut(SmearEvent::ECCE_CHERENKOV_EMC,
                                                 2., iElec, iPosi);
        } else if (pid_system == 3) {
          right_signs = smear_event.bool4JpsiCut(SmearEvent::EICC_ACCEPTANCE,
                                                 2., iElec, iPosi);
        } else if (pid_system == 4) {
          right_signs = smear_event.bool4JpsiCut(SmearEvent::PERFECT_PID, 2.,
                                                 iElec, iPosi);
        } else if (pid_system == 5) {
          right_signs =
              smear_event.bool4JpsiCut(SmearEvent::NO_PID, 2., iElec, iPosi);
        }
        if (right_signs) {
          cutstr_Jpsi = smear_event.CalCut_Jpsi(iElec, iPosi);
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
  cut_cal_Jpsi(name_input, atoi(argv[2]), argv[3], argv[4]);
  // cut_cal_d0("event_smeared_Det_v1.root", 1,
  // "cut_event_smeared_Det_v1.root");
  return 0;
}
