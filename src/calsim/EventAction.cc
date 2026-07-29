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
  const int nx = fP.crystals_nx;

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
    if (eattGeV > fP.hit_threshold_gev) {
      evis = eattGeV * (fP.smearing_on ? fP.energy_scale : 1.0);
      if (fP.smearing_on && evis > 0.) {
        // extra smearing = quadrature difference between the measured resolution model
        // and the intrinsic Geant4 part: sigma^2/E^2 = stochastic^2/E + noise^2/E^2 + constant^2
        const double E = evis;
        const double var = fP.measured_res_stochastic * fP.measured_res_stochastic / E
                         + fP.measured_res_noise * fP.measured_res_noise / (E * E)
                         + fP.measured_res_constant * fP.measured_res_constant
                         - fP.intrinsic_res_stochastic * fP.intrinsic_res_stochastic / E
                         - fP.intrinsic_res_noise * fP.intrinsic_res_noise / (E * E);
        if (var > 0.) evis *= 1.0 + G4RandGauss::shoot(0., std::sqrt(var));
        if (evis < 0.) evis = 0.;
      }
    }

    if (edepGeV < fP.store_min_gev && evis <= 0.) continue;
    const int ix = cell % nx, iy = cell / nx;
    if (fP.timing_on) {
      double tvis = (h.eatt > 0.) ? h.twsum / h.eatt / ns : 0.;
      tvis += G4RandGauss::shoot(0., fP.time_sigma_ns);
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
