#include "McParticles_byEvent.h"

void ParticleDecay(string ver_det, string name_inputfile, string name_output,
                   int pid, int seed, int n_required_events) {
  //   Pythia pythia;
  // event = pythia.event;
  // settings = pythia.settings;
  // pdt = pythia.particleData;
  ReadMap(ver_det);
  TFile *file_input = new TFile(name_inputfile.c_str());
  TTree *tree = (TTree *)file_input->Get("production_per_event");
  McParticles_byEvent *mcparticle_decay =
      new McParticles_byEvent(name_output.c_str(), tree, seed);
  mcparticle_decay->Run(n_required_events,pid);
}

int main(int argc, char *argv[]) {
  ParticleDecay(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]),
                atoi(argv[6]));
}
