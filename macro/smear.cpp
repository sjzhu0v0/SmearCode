#include "SmearEvent.h"

void smear(string InFile,string ver_det) {
  ReadMap(ver_det);
  string filePath = "/ustcfs/HICUser/ephy/work/sjzhu/Singularity-master/cvmfs/work/pythia6_generator/3.5_20/generator"+InFile+"/mc.root";
  TTree *t_event = (TTree *)(new TFile(InFile.c_str()))->Get("EICTree");

  erhic::EventMC *mcEvent(NULL);
  t_event->SetBranchAddress("event", &mcEvent);

  SmearEvent *event_smearing = new SmearEvent(false, ("event_smeared_"+ver_det+".root").c_str());

  cout << "Smearing Begin.  " << endl;

  for (long iEvent = 0; iEvent < 10000/* t_event->GetEntries() */; iEvent++)
  {
    if (iEvent % 10000 == 0) {
      cout << iEvent << " events" << endl;
    }
    t_event->GetEntry(iEvent);
    event_smearing->EventInput(mcEvent,iEvent);
    if (event_smearing->IsTriggered()) {
      event_smearing->Fill();
    }
  }
  event_smearing->Save();
}

int main(int argc, char *argv[]) {
  // const string FilePath = argv[1];
  // const string ver_det = argv[2];
  smear("/data/work/Eic/data/mc.root","Det_v3");
  return 0;
}
