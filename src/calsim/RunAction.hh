#ifndef CALSIM_RUNACTION_HH
#define CALSIM_RUNACTION_HH

#include "G4UserRunAction.hh"
#include "Params.hh"
#include <fstream>

// Opens the per-job output CSVs (hits.csv, events.csv) and writes run.json at the end.
class RunAction : public G4UserRunAction {
public:
  explicit RunAction(const Params& p) : fP(p) {}
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  std::ofstream& hitsFile()   { return fHits; }
  std::ofstream& eventsFile() { return fEvents; }

private:
  const Params& fP;
  std::ofstream fHits, fEvents;
};

#endif
