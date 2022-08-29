#include "TFile.h"
#include "TTree.h"
#include "eicsmear/erhic/EventMC.h"
#include "iostream"
#include "string"

using namespace std;

class ParticleSelector {
public:
  ParticleSelector(erhic::EventMC *mcEvent) {
    for (int iTrack = 0; iTrack < mcEvent->GetNTracks(); iTrack++) {
      erhic::ParticleMC *particle = mcEvent->GetTrack(iTrack);
      if (particle->GetStatus() == 21) {
        n_init_particles = iTrack;
      }
    }
  }
  ~ParticleSelector() {}

  bool isParticleDesireddHadron(const erhic::ParticleMC *mcParticle) {
    return isParticleProducedHadron(mcParticle) &&
           !isParticleProducedHadron(mcParticle->GetParent()) &&
           abs(mcParticle->GetParent()->id) != 15;
  }

private:
  int n_init_particles = 0;
  bool isParticleProducedHadron(const erhic::ParticleMC *mcParticle) {
    if (mcParticle->GetIndex() - 1 <= n_init_particles ||
        !isParticleHadron(mcParticle)) {
      return false;
    } else {
      if (mcParticle->id == 2212 &&
          mcParticle->GetParent()->GetIndex() <= n_init_particles) {
        return false;
      } else {
        return true;
      }
    }
  }

  bool isParticleHadron(const erhic::ParticleMC *mcParticle) {
    int idSave = abs(mcParticle->id);
    if (idSave <= 100 || (idSave >= 1000000 && idSave <= 9000000) ||
        idSave >= 9900000)
      return false;
    if (idSave == 130 || idSave == 310)
      return true;
    if (idSave % 10 == 0 || (idSave / 10) % 10 == 0 || (idSave / 100) % 10 == 0)
      return false;
    return true;
  }
};

void Hadron_Mc_info(string InFile, string OutFile) {
  TFile *f_output = new TFile(OutFile.c_str(), "recreate");
  TTree *t_out = new TTree("production_per_event", "prodution_per_event");
  string inputfile = "/data/work/Eic/data/mc.root";
  /*  "/ustcfs/HICUser/ephy/work/sjzhu/Singularity-master/cvmfs/work/"
   "pythia6_generator/3.5_20/generator" +
   InFile + "/mc.root"; */
  TTree *t_event = (TTree *)(new TFile(inputfile.c_str()))->Get("EICTree");

  int num_particles = 0;
  double particles_px[40] = {0.};
  double particles_py[40] = {0.};
  double particles_pz[40] = {0.};
  double particles_vx[40] = {0.};
  double particles_vy[40] = {0.};
  double particles_vz[40] = {0.};
  double particles_z[40] = {0.};
  int particles_pid[40] = {0};
  double p3_scattered_lepton[3] = {0.};
  double Q2 = 0;
  double x = 0;
  double y = 0;

  erhic::EventMC *mcEvent(NULL);

  t_event->SetBranchAddress("event", &mcEvent);

  t_out->Branch("num_particles", &num_particles, "num_particles/I");
  t_out->Branch("particles_px", particles_px, "particles_px[num_particles]/D");
  t_out->Branch("particles_py", particles_py, "particles_py[num_particles]/D");
  t_out->Branch("particles_pz", particles_pz, "particles_pz[num_particles]/D");
  t_out->Branch("particles_vx", particles_vx, "particles_vx[num_particles]/D");
  t_out->Branch("particles_vy", particles_vy, "particles_vy[num_particles]/D");
  t_out->Branch("particles_vz", particles_vz, "particles_vz[num_particles]/D");
  t_out->Branch("particles_z", particles_z, "particles_z[num_particles]/D");
  t_out->Branch("particles_pid", particles_pid,
                "particles_pid[num_particles]/I");
  t_out->Branch("p3_scattered_lepton", p3_scattered_lepton,
                "p3_scattered_lepton[3]/D");
  t_out->Branch("Q2", &Q2, "Q2/D");
  t_out->Branch("x", &x, "x/D");
  t_out->Branch("y", &y, "y/D");

  for (long iEvent = 0; iEvent < t_event->GetEntries(); iEvent++) {
    t_event->GetEntry(iEvent);
    if (iEvent % (t_event->GetEntries() / 10) == 0) {
      cout << iEvent / (t_event->GetEntries() / 10) << "0% completed" << endl;
    }

    Q2 = mcEvent->GetQ2();
    x = mcEvent->GetX();
    y = mcEvent->GetY();
    ParticleSelector particle_selector(mcEvent);
    num_particles = 0;
    for (int iTrack = 0; iTrack < mcEvent->GetNTracks(); iTrack++) {
      const erhic::ParticleMC *mcParticle = mcEvent->GetTrack(iTrack);
      if (particle_selector.isParticleDesireddHadron(mcParticle)) {
        particles_px[num_particles] = mcParticle->Get4Vector().Px();
        particles_py[num_particles] = mcParticle->Get4Vector().Py();
        particles_pz[num_particles] = mcParticle->Get4Vector().Pz();
        particles_vx[num_particles] = mcParticle->GetVertex().x();
        particles_vy[num_particles] = mcParticle->GetVertex().y();
        particles_vz[num_particles] = mcParticle->GetVertex().z();
        particles_z[num_particles] = mcParticle->GetZ();
        particles_pid[num_particles] = mcParticle->id;
        if (abs(mcParticle->id) < 110) {
          cout << mcParticle->GetIndex() << endl;
          mcEvent->Print();
          return;
        }
        num_particles++;
      }
    }

    p3_scattered_lepton[0] = mcEvent->ScatteredLepton()->Get4Vector().Px();
    p3_scattered_lepton[1] = mcEvent->ScatteredLepton()->Get4Vector().Py();
    p3_scattered_lepton[2] = mcEvent->ScatteredLepton()->Get4Vector().Pz();
    t_out->Fill();
  }

  f_output->cd();
  t_out->Write();
  f_output->Close();
}

int main(int argc, char *argv[]) {
  const string InFile = argv[1];
  const string OutFile = argv[2];
  Hadron_Mc_info(InFile, OutFile);
  return 0;
}
