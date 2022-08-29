#include "McParticles_byEvent.h"

const string xmlDB = "/home/sjzhu/work/Software/pythia8/share/Pythia8/xmldoc";
Pythia pythia(xmlDB.c_str());
Event &event = pythia.event;
Settings &settings = pythia.settings;
ParticleData &pdt = pythia.particleData;

McParticles_byEvent::McParticles_byEvent(string fileName, TTree *mctree,
                                         int seed)
    : SmearEvent(false, fileName) {
  gRandom->SetSeed(seed);
  InputMcData(mctree);
}

McParticles_byEvent::~McParticles_byEvent() {}

void McParticles_byEvent::fillParticle(int id, float px, float py, float pz,
                                       Event &event, ParticleData &pdt,
                                       Rndm &rndm, bool atRest,
                                       bool hasLifetime) {
  // Reset event record to allow for new event.
  event.reset();

  // Select particle mass; where relevant according to Breit-Wigner.
  float mm = pdt.mSel(id);

  float ee;
  // Special case when particle is supposed to be at rest.
  if (atRest) {
    ee = mm;
  }
  ee = sqrt(mm * mm + px * px + py * py + pz * pz);
  // Angles as input or uniform in solid angle.

  // Store the particle in the event record.
  int iNew = event.append(id, 1, 0, 0, px, py, pz, ee, mm);

  // Generate lifetime, to give decay away from primary vertex.
  if (hasLifetime)
    event[iNew].tau(event[iNew].tau0() * rndm.exp());
  // cout << event[iNew].tau0() << endl;
}

void McParticles_byEvent::InputMcData(TTree *tree) {
  if (!tree)
    return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("num_particles", &num_particles);
  fChain->SetBranchAddress("particles_px", particles_px);
  fChain->SetBranchAddress("particles_py", particles_py);
  fChain->SetBranchAddress("particles_pz", particles_pz);
  fChain->SetBranchAddress("particles_pid", particles_pid);
  fChain->SetBranchAddress("p3_scattered_lepton", p3_scattered_lepton);
  fCurrent = fChain->GetEntries();
  cout << fCurrent << " event samples used" << endl;
}

void McParticles_byEvent::PythiaInit() {
  pythia.readString(
      "ProcessLevel:all = off"); /*key switchs in examples/main21.cc*/
  pythia.readString("421:onMode = off");
  pythia.readString("421:onIfMatch = 321 -211");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed =  21112311");
  pythia.readString("4122:onMode = off");
  pythia.readString("4122:onIfMatch = 2224 -321");
  pythia.readString("443:onMode = off");
  pythia.readString("443:onIfMatch = 11 -11");

  pythia.init();
}

void McParticles_byEvent::ParticleDecay_And_Smear(int pid) {
  fChain->GetEntry((int)gRandom->Uniform(0, fCurrent));
  for (int iParticle = 0; iParticle < num_particles; iParticle++) {
    if (abs(particles_pid[iParticle]) == pid) {
      fillParticle(particles_pid[iParticle], particles_px[iParticle],
                   particles_py[iParticle], particles_pz[iParticle], event, pdt,
                   pythia.rndm, false, true);
      if (pythia.next()) {
        detected_particles.clear();
        num_final_particles = 0;
        is_triggered = true;
        for (int iTrack = 0; iTrack < event.size(); iTrack++) {
          if (event[iTrack].status() > 0) {
            TVector3 p3particle(event[iTrack].px(), event[iTrack].py(),
                                event[iTrack].pz());
            // cout << p3particle.PseudoRapidity() << "  " << p3particle.Mag()
            // << endl;
            TVector3 v3particle(event[iTrack].xProd(), event[iTrack].yProd(),
                                event[iTrack].zProd());
            SmearParticleStr ParticleStr =
                SmearFunction(p3particle, v3particle, event[iTrack].id());

            is_triggered = is_triggered && ParticleStr.isDetected;
            if (is_triggered) {
              detected_particles.push_back(ParticleStr);
              num_final_particles++;
            } else {
              break;
            }
          }
        }
        if (is_triggered) {
          is_scattered_electron_detected = false;
          smeared_px_scattered_lepton = smeared_py_scattered_lepton =
              smeared_pz_scattered_lepton = 0.;
          true_px_scattered_lepton = true_py_scattered_lepton =
              true_pz_scattered_lepton = 0.;
          TVector3 p3lepton(p3_scattered_lepton);
          double pseudorapidity_scattered_lepton = p3lepton.PseudoRapidity();

          if (pseudorapidity_scattered_lepton > -3.5 &&
              pseudorapidity_scattered_lepton < 3.5) {
            SmearParticleStr scattered_lepton = SmearFunction(p3lepton,
                                                              {
                                                                  0.,
                                                                  0.,
                                                                  0,
                                                              },
                                                              11);
            if (scattered_lepton.isDetected) {
              is_scattered_electron_detected = true;
              int_is_scattered_electron_detected = 1;
              smeared_px_scattered_lepton = scattered_lepton.p_smeared.X();
              smeared_py_scattered_lepton = scattered_lepton.p_smeared.Y();
              smeared_pz_scattered_lepton = scattered_lepton.p_smeared.Z();
            }
          }

          true_px_scattered_lepton = p3_scattered_lepton[0];
          true_py_scattered_lepton = p3_scattered_lepton[1];
          true_pz_scattered_lepton = p3_scattered_lepton[2];

          prime_vx = mvd_dp.GetRandomVx(7);
          prime_vy = mvd_dp.GetRandomVy(7);
          prime_vz = mvd_dp.GetRandomVz(7);

          for (int iDetectedParticle = 0;
               iDetectedParticle < num_final_particles; iDetectedParticle++) {
            true_px[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).p_mc.x();
            true_py[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).p_mc.y();
            true_pz[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).p_mc.z();

            smeared_px[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).p_smeared.x();
            smeared_py[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).p_smeared.y();
            smeared_pz[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).p_smeared.z();

            true_pid[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).pid;
            smeared_dcaxy[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).smeared_dcaxy;
            smeared_dcaz[iDetectedParticle] =
                detected_particles.at(iDetectedParticle).smeared_dcaz;

            for (int iParticleType = 0; iParticleType < n_particle_type;
                 iParticleType++) {
              for (int iDetType = 0; iDetType < n_particle_type; iDetType++) {
                n_sigma_pid_ecce[iDetectedParticle][iParticleType][iDetType] =
                    detected_particles.at(iDetectedParticle)
                        .n_sigma_pid_ecce[iParticleType][iDetType];
                n_sigma_pid_eicc[iDetectedParticle][iParticleType][iDetType] =
                    detected_particles.at(iDetectedParticle)
                        .n_sigma_pid_eicc[iParticleType][iDetType];
              }
            }
          }
          SmearEvent::Fill();
          N_events_gotten++;
        } else {
          continue;
        }
      }
    }
  }
}

void McParticles_byEvent::Run(long n_required_events, int pid) {
  PythiaInit();
  while (N_events_gotten < n_required_events) {
    ParticleDecay_And_Smear(pid);
  }
  SmearEvent::Save();
}
