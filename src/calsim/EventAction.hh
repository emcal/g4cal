#ifndef CALSIM_EVENTACTION_HH
#define CALSIM_EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "Params.hh"

class PrimaryGeneratorAction;
class RunAction;
class DetectorConstruction;

// End of event: apply response stages (threshold on attenuated energy -> en_scale ->
// mcsmear-style smearing) and write the hit and truth CSV rows.
class EventAction : public G4UserEventAction {
public:
  EventAction(const Params& p, const PrimaryGeneratorAction* gen, RunAction* run,
              const DetectorConstruction* det)
    : fP(p), fGen(gen), fRun(run), fDet(det) {}

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

private:
  const Params& fP;
  const PrimaryGeneratorAction* fGen;
  RunAction* fRun;
  const DetectorConstruction* fDet;
};

#endif
