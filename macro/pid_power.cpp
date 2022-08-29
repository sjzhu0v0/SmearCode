#include "PidMap.h"

void pid_power() {
  gStyle->SetOptStat(0);
  gStyle->SetLabelFont(132, "X");
  gStyle->SetLabelFont(132, "Y");
  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetLabelSize(0.05, "Y");
  gStyle->SetTitleSize(0.05, "X");
  gStyle->SetTitleSize(0.05, "Y");
  ReadMap("Det_v1");

  TH2F* h2_elec2elec[3];
  TH2F* h2_pion2elec[3];

  TH2F* h2_pion2kaon[3];
  TH2F* h2_kaon2kaon[3];

  TH2F* h2_kaon2proton[3];
  TH2F* h2_proton2proton[3];

  TH2F* h2_mc = new TH2F("h2_mc", "h2_mc", 100, 0, 20, 100, -3.5, 3.5);
  TH2F* h2_eff_elec2elec = new TH2F("h2_eff_elec2elec", "h2_eff_elec2elec", 100,
                                    0, 20, 100, -3.5, 3.5);
  TH2F* h2_eff_pion2elec = new TH2F("h2_eff_pion2elec", "h2_eff_pion2elec", 100,
                                    0, 20, 100, -3.5, 3.5);
  TH2F* h2_eff_pion2kaon = new TH2F("h2_eff_pion2kaon", "h2_eff_pion2kaon", 100,
                                    0, 20, 100, -3.5, 3.5);
  TH2F* h2_eff_kaon2kaon = new TH2F("h2_eff_kaon2kaon", "h2_eff_kaon2kaon", 100,
                                    0, 20, 100, -3.5, 3.5);
  TH2F* h2_eff_kaon2proton = new TH2F(
      "h2_eff_kaon2proton", "h2_eff_kaon2proton", 100, 0, 20, 100, -3.5, 3.5);
  TH2F* h2_eff_proton2proton =
      new TH2F("h2_eff_proton2proton", "h2_eff_proton2proton", 100, 0, 20, 100,
               -3.5, 3.5);

  for (int iEta = 0; iEta < 3; iEta++) {
    h2_elec2elec[iEta] =
        new TH2F(Form("h2_elec2elec_%d", iEta), Form("h2_elec2elec_%d", iEta),
                 100, 0, 20, 100, 0, 1);
    h2_pion2elec[iEta] =
        new TH2F(Form("h2_pion2elec_%d", iEta), Form("h2_pion2elec_%d", iEta),
                 100, 0, 20, 100, 0, 1);
    h2_pion2kaon[iEta] =
        new TH2F(Form("h2_pion2kaon_%d", iEta), Form("h2_pion2kaon_%d", iEta),
                 100, 0, 20, 100, 0, 1);
    h2_kaon2kaon[iEta] =
        new TH2F(Form("h2_kaon2kaon_%d", iEta), Form("h2_kaon2kaon_%d", iEta),
                 100, 0, 20, 100, 0, 1);
    h2_kaon2proton[iEta] =
        new TH2F(Form("h2_kaon2proton_%d", iEta),
                 Form("h2_kaon2proton_%d", iEta), 100, 0, 20, 100, 0, 1);
    h2_proton2proton[iEta] =
        new TH2F(Form("h2_proton2proton_%d", iEta),
                 Form("h2_proton2proton_%d", iEta), 100, 0, 20, 100, 0, 1);
  }

  for (int iEtaRegion = 0; iEtaRegion < 3; iEtaRegion++) {
    for (int iEta = 0; iEta < 1000; iEta++) {
      for (int iP = 0; iP < 1000; iP++) {
        double momentum = MRndgen.Uniform(0.01, 25);
        double eta = MRndgen.Uniform(min_cherenkov_eta_region[iEtaRegion],
                                     max_cherenkov_eta_region[iEtaRegion]);
        PidMap pidmap_elec(PidMap::electron, eta, momentum);
        PidMap pidmap_pion(PidMap::pion, eta, momentum);
        PidMap pidmap_kaon(PidMap::kaon, eta, momentum);
        PidMap pidmap_proton(PidMap::proton, eta, momentum);

        h2_mc->Fill(momentum, eta);
        if (pidmap_elec.get_expected_pid_int((PidMap::PID_SYSTEM)2) ==11) {
          h2_eff_elec2elec->Fill(momentum, eta);
        }
        if (pidmap_pion.get_expected_pid_int((PidMap::PID_SYSTEM)2) == 11) {
          h2_eff_pion2elec->Fill(momentum, eta);
        }
        if (pidmap_kaon.get_expected_pid_int((PidMap::PID_SYSTEM)2) == 321) {
          h2_eff_kaon2kaon->Fill(momentum, eta);
        }
        if (pidmap_pion.get_expected_pid_int((PidMap::PID_SYSTEM)2)==321) {
          h2_eff_pion2kaon->Fill(momentum, eta);
        }
        if (pidmap_kaon.get_expected_pid_int((PidMap::PID_SYSTEM)2) ==2212 ) {
          h2_eff_kaon2proton->Fill(momentum, eta);
        }
        if (pidmap_proton.get_expected_pid_int((PidMap::PID_SYSTEM)2) ==2212) {
          h2_eff_proton2proton->Fill(momentum, eta);
        }

        h2_elec2elec[iEtaRegion]->Fill(
            momentum,
            pidmap_elec.get_probability_to_particle(PidMap::electron));
        h2_pion2elec[iEtaRegion]->Fill(
            momentum,
            pidmap_pion.get_probability_to_particle(PidMap::electron));
        h2_pion2kaon[iEtaRegion]->Fill(
            momentum, pidmap_pion.get_probability_to_particle(PidMap::kaon));
        h2_kaon2kaon[iEtaRegion]->Fill(
            momentum, pidmap_kaon.get_probability_to_particle(PidMap::kaon));
        h2_kaon2proton[iEtaRegion]->Fill(
            momentum, pidmap_kaon.get_probability_to_particle(PidMap::proton));
        h2_proton2proton[iEtaRegion]->Fill(
            momentum,
            pidmap_proton.get_probability_to_particle(PidMap::proton));
      }
    }
  }

  h2_eff_elec2elec->Divide(h2_mc);
  h2_eff_pion2elec->Divide(h2_mc);
  h2_eff_kaon2kaon->Divide(h2_mc);
  h2_eff_pion2kaon->Divide(h2_mc);
  h2_eff_proton2proton->Divide(h2_mc);
  h2_eff_kaon2proton->Divide(h2_mc);

  TCanvas* c_pid_efficiency =
      new TCanvas("c_pid_efficiency", "c_pid_efficiency", 1200, 800);
  c_pid_efficiency->Divide(3, 2);
  c_pid_efficiency->cd(1);
  h2_eff_elec2elec->Draw("colz");
  c_pid_efficiency->cd(4);
  h2_eff_pion2elec->Draw("colz");
  c_pid_efficiency->cd(2);
  h2_eff_kaon2kaon->Draw("colz");
  c_pid_efficiency->cd(5);
  h2_eff_pion2kaon->Draw("colz");
  c_pid_efficiency->cd(3);
  h2_eff_proton2proton->Draw("colz");
  c_pid_efficiency->cd(6);
  h2_eff_kaon2proton->Draw("colz");
   c_pid_efficiency->SaveAs("c_pid_efficiency.png");
}

int main() {
  pid_power();
  return 0;
}
