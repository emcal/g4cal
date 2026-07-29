#include "RunAction.hh"
#include "G4Run.hh"
#include <ctime>
#include <iostream>

void RunAction::BeginOfRunAction(const G4Run*) {
  fHits.open(fP.out_dir + "/hits.csv");
  fEvents.open(fP.out_dir + "/events.csv");
  if (!fHits || !fEvents) {
    std::cerr << "calsim: cannot open output CSVs in " << fP.out_dir << "\n";
    std::exit(3);
  }
  fHits << "event,ix,iy,edep_true,e_vis";
  if (fP.timing_on) fHits << ",t_vis";
  fHits << "\n";
  fEvents << "event,pdg,ekin,x,y,dx,dy,dz\n";
}

void RunAction::EndOfRunAction(const G4Run* run) {
  fHits.close();
  fEvents.close();

  std::ofstream j(fP.out_dir + "/run.json");
  std::string s = fP.json();
  // splice processed-event count and timestamp into the params JSON
  const auto pos = s.rfind("\n}");
  std::time_t t = std::time(nullptr);
  char ts[32];
  std::strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%S", std::gmtime(&t));
  s.insert(pos, ",\n  \"n_processed\": " + std::to_string(run->GetNumberOfEvent())
                + ",\n  \"timestamp\": \"" + ts + "\"");
  j << s;
}
