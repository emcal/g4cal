#include "Params.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

// calsim: standalone PbWO4 array simulation. FTFP_BERT + EM option4 ("FTFP_BERT_EMZ").
// Writes hits.csv, events.csv, run.json into --out-dir. See --help in Params.hh parse.
int main(int argc, char** argv) {
  Params p = Params::parse(argc, argv);

  G4Random::setTheEngine(new CLHEP::MixMaxRng());
  G4Random::setTheSeed(p.seed);

  auto* rm = new G4RunManager();
  rm->SetVerboseLevel(0);

  auto* det = new DetectorConstruction(p);
  rm->SetUserInitialization(det);

  G4PhysListFactory factory;
  auto* phys = factory.GetReferencePhysList("FTFP_BERT_EMZ");
  phys->SetVerboseLevel(0);
  rm->SetUserInitialization(phys);

  auto* gen = new PrimaryGeneratorAction(p);
  auto* runAct = new RunAction(p);
  rm->SetUserAction(gen);
  rm->SetUserAction(runAct);
  rm->SetUserAction(new EventAction(p, gen, runAct, det));

  auto* ui = G4UImanager::GetUIpointer();
  for (const char* c : {"/control/verbose 0", "/run/verbose 0", "/event/verbose 0",
                        "/tracking/verbose 0", "/process/em/verbose 0",
                        "/process/had/verbose 0"})
    ui->ApplyCommand(c);

  rm->Initialize();
  rm->BeamOn(p.n_events);

  delete rm;
  return 0;
}
