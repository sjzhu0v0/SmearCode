#include "SmearEvent.h"

void cut_cal(string name_input, int pid_system, string name_output,
             string name_mcinput, string is_signal_cal) {
  bool is_signal = false;
  if (is_signal_cal == "Signal") {
    is_signal = true;
  }

  SmearEvent::CutStr_lambda_c cutstr_lambda_c;
  TFile *f_mc = new TFile(name_mcinput.c_str());
  SmearEvent smear_event(true, name_input);
  TFile *f_cut = new TFile(name_output.c_str(), "recreate");
  TTree *t_event = (TTree *)f_mc->Get("EICTree");
  cout << t_event << endl;
  smear_event.McTreeInitializing(t_event);
  TTree *t_cut = new TTree("Cut", "Cut");
  t_cut->Branch("inv_mass", &cutstr_lambda_c.inv_mass, "inv_mass/D");
  t_cut->Branch("dca_proton", &cutstr_lambda_c.dca_proton, "dca_proton/D");
  t_cut->Branch("dca_kaon", &cutstr_lambda_c.dca_kaon, "dca_kaon/D");
  t_cut->Branch("dca_pion", &cutstr_lambda_c.dca_pion, "dca_pion/D");
  t_cut->Branch("cos_theta", &cutstr_lambda_c.cos_theta, "cos_theta/D");
  t_cut->Branch("decay_length_lambda_c", &cutstr_lambda_c.decay_length_lambda_c,
                "decay_length_lambda_c/D");
  t_cut->Branch("maxDca_daughters", &cutstr_lambda_c.maxDca_daughters,
                "maxDca_daughters/D");
  t_cut->Branch("dca_lambda_c2pv", &cutstr_lambda_c.dca_lambda_c2pv,
                "dca_lambda_c2pv/D");
  t_cut->Branch("pT_lambda_c", &cutstr_lambda_c.pT_lambda_c, "pT_lambda_c/D");
  t_cut->Branch("rapidity_lambda_c", &cutstr_lambda_c.rapidity_lambda_c,
                "rapidity_lambda_c/D");
  t_cut->Branch("p3_mc", cutstr_lambda_c.p3_mc, "p3_mc[3]/D");
  t_cut->Branch("p3_elec", cutstr_lambda_c.p3_elec, "p3_elec[3]/D");
  t_cut->Branch("Lambda_c_Sign", &cutstr_lambda_c.Lambda_c_Sign,
                "Lambda_c_Sign/I");
  t_cut->Branch("num_final_particles", &cutstr_lambda_c.num_final_particles,
                "num_final_particles/I");
  t_cut->Branch("num_mc_final_charged_particles",
                &cutstr_lambda_c.num_mc_final_charged_particles,
                "num_mc_final_charged_particles/I");

  bool right_signs;

  for (long iEvent = 0; iEvent < smear_event.GetEntries(); iEvent++) {
    smear_event.GetEntry(iEvent);

    if (is_signal_cal != "all") {
      if (is_signal ? !smear_event.isSignal(421) : smear_event.isSignal(4122))
        continue;
    }

    int temp_var2 = smear_event.num_final_particles;
    for (int iPion = 0; iPion < smear_event.num_final_particles; iPion++) {
      for (int iKaon = 0; iKaon < smear_event.num_final_particles; iKaon++) {
        for (int iProton = 0; iProton < smear_event.num_final_particles;
             iProton++) {
          if (iPion == iKaon || iPion == iProton || iKaon == iProton) continue;
          if (pid_system == 0) {
            right_signs = smear_event.bool4Lambda_cCut(SmearEvent::ECCE_PID, 2.,
                                                       iPion, iKaon, iProton);
          } else if (pid_system == 1) {
            right_signs = smear_event.bool4Lambda_cCut(
                SmearEvent::ECCE_TOF_EMC, 2., iPion, iKaon, iProton);
          } else if (pid_system == 2) {
            right_signs = smear_event.bool4Lambda_cCut(
                SmearEvent::ECCE_CHERENKOV_EMC, 2., iPion, iKaon, iProton);
          } else if (pid_system == 3) {
            right_signs = smear_event.bool4Lambda_cCut(
                SmearEvent::EICC_ACCEPTANCE, 2., iPion, iKaon, iProton);
          } else if (pid_system == 4) {
            right_signs = smear_event.bool4Lambda_cCut(
                SmearEvent::PERFECT_PID, 2., iPion, iKaon, iProton);
          } else if (pid_system == 5) {
            right_signs = smear_event.bool4Lambda_cCut(SmearEvent::NO_PID, 2.,
                                                       iPion, iKaon, iProton);
          }
          if (right_signs) {
            cutstr_lambda_c =
                smear_event.CalCut_Lambda_c(iPion, iKaon, iProton);
            t_cut->Fill();
          }
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
  cut_cal(argv[1], atoi(argv[2]), argv[3], argv[4], argv[5]);
  // cut_cal("~/HIUser_sjzhu/work3/code/Smear_Code/bin/lambda_c_decay_test.root",
  //         3, "cut_lambda_c_test.root",
  //         "/ustcfs/HICUser/ephy/work/sjzhu/Singularity-master/cvmfs/work/"
  //         "pythia6_generator/3.5_20/generator1/mc.root");
  // cut_cal_d0("event_smeared_Det_v1.root", 1,
  // "cut_event_smeared_Det_v1.root");
  return 0;
}
