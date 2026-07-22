#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "CrystalSD.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <cstdio>

void EventAction::BeginOfEventAction(const G4Event*) {
  fDet->sd()->BeginEvent();
}

void EventAction::EndOfEventAction(const G4Event* evt) {
  const int id = evt->GetEventID();
  const int nx = fP.nx;

  char line[256];
  std::snprintf(line, sizeof line, "%d,%d,%.6g,%.4f,%.4f,0,0,1\n",
                id, fGen->pdg(), fGen->ekinGeV(), fGen->xmm(), fGen->ymm());
  fRun->eventsFile() << line;

  auto& hits = fRun->hitsFile();
  for (const auto& [cell, h] : fDet->sd()->hits()) {
    const double edepGeV = h.edep / GeV;
    const double eattGeV = h.eatt / GeV;

    // response chain order per halld: attenuate (done in SD) -> 5 MeV threshold on the
    // attenuated unscaled energy (hitECAL.c:275) -> scale -> smear (ECALSmearer.cc:80-97)
    double evis = 0.;
    if (eattGeV > fP.thresh_att_gev) {
      evis = eattGeV * (fP.smear_on ? fP.en_scale : 1.0);
      if (fP.smear_on && evis > 0.) {
        const double E = evis;
        const double var = fP.meas_a * fP.meas_a / E + fP.meas_b * fP.meas_b / (E * E)
                         + fP.meas_c * fP.meas_c
                         - fP.intr_a * fP.intr_a / E - fP.intr_b * fP.intr_b / (E * E);
        if (var > 0.) evis *= 1.0 + G4RandGauss::shoot(0., std::sqrt(var));
        if (evis < 0.) evis = 0.;
      }
    }

    if (edepGeV < fP.store_min_gev && evis <= 0.) continue;
    const int ix = cell % nx, iy = cell / nx;
    if (fP.time_on) {
      double tvis = (h.eatt > 0.) ? h.twsum / h.eatt / ns : 0.;
      tvis += G4RandGauss::shoot(0., fP.t_sigma_ns);
      std::snprintf(line, sizeof line, "%d,%d,%d,%.6g,%.6g,%.4f\n",
                    id, ix, iy, edepGeV, evis, tvis);
    } else {
      std::snprintf(line, sizeof line, "%d,%d,%d,%.6g,%.6g\n", id, ix, iy, edepGeV, evis);
    }
    hits << line;
  }

  if ((id + 1) % 1000 == 0 || id + 1 == fP.n_events)
    std::printf("progress %d/%d\n", id + 1, fP.n_events);
}
