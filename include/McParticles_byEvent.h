#ifndef __McParticles_byEvent_h__
#define __McParticles_byEvent_h__

#include "Pythia8/Pythia.h"
#include "SmearEvent.h"

using namespace Pythia8;

extern const string xmlDB;
extern Pythia pythia;
extern Event &event;
extern Settings &settings;
extern ParticleData &pdt;

class McParticles_byEvent : public SmearEvent {
private:
  TTree *fChain;

  long fCurrent = 0;
  long N_events_gotten = 0;

  void InputMcData(TTree *tree_mcinput);
  void InputMcData(TFile *file_mcinput);
  void fillParticle(int id, float px, float py, float pz, Event &event,
                    ParticleData &pdt, Rndm &rndm, bool atRest = false,
                    bool hasLifetime = true);
  void ParticleDecay_And_Smear(int pid);
  void PythiaInit();

public:
  Int_t num_particles = 0;
  Double_t particles_px[40] = {0.}; //[num_particles]
  Double_t particles_py[40] = {0.}; //[num_particles]
  Double_t particles_pz[40] = {0.}; //[num_particles]
  Int_t particles_pid[40] = {0};    //[num_particles]
  Double_t p3_scattered_lepton[3] = {0.};

  McParticles_byEvent(string fileName, TTree *mctree, int seed);
  ~McParticles_byEvent();

  void Run(long n_required_events, int pid);
};

#endif
