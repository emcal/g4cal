#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cstdlib>
#include <iostream>

PrimaryGeneratorAction::PrimaryGeneratorAction(const Params& p) : fP(p) {
  fGun = std::make_unique<G4ParticleGun>(1);
  auto* def = G4ParticleTable::GetParticleTable()->FindParticle(fP.particle);
  if (!def) {
    std::cerr << "calsim: unknown particle " << fP.particle << "\n";
    std::exit(2);
  }
  fGun->SetParticleDefinition(def);
  fGun->SetParticleMomentumDirection({0., 0., 1.});
  fPdg = def->GetPDGEncoding();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
  const double pitch = fP.pitch_mm();   // mm

  fEkin = (fP.energy_min_gev == fP.energy_max_gev)
            ? fP.energy_min_gev
            : fP.energy_min_gev + (fP.energy_max_gev - fP.energy_min_gev) * G4UniformRand();

  // center of the "central" crystal (for even nx this is cell nx/2, off-center by pitch/2)
  const double cx0 = (fP.crystals_nx / 2 - 0.5 * (fP.crystals_nx - 1)) * pitch;
  const double cy0 = (fP.crystals_ny / 2 - 0.5 * (fP.crystals_ny - 1)) * pitch;

  if (fP.gun_mode == "face") {
    fX = (G4UniformRand() - 0.5) * fP.crystals_nx * pitch;
    fY = (G4UniformRand() - 0.5) * fP.crystals_ny * pitch;
  } else if (fP.gun_mode == "central") {
    fX = cx0 + (G4UniformRand() - 0.5) * pitch;
    fY = cy0 + (G4UniformRand() - 0.5) * pitch;
  } else if (fP.gun_mode == "point") {
    fX = fP.gun_x_mm;
    fY = fP.gun_y_mm;
  } else if (fP.gun_mode == "grid") {
    const int n = fP.gun_grid_n;
    const long i = fCount % (long(n) * n);
    const int gx = int(i % n), gy = int(i / n);
    fX = cx0 + ((gx + 0.5) / n - 0.5) * pitch;
    fY = cy0 + ((gy + 0.5) / n - 0.5) * pitch;
  } else {
    std::cerr << "calsim: unknown gun mode " << fP.gun_mode << "\n";
    std::exit(2);
  }
  ++fCount;

  const double z0 = -0.5 * fP.crystal_length_mm - 200.0;   // 200 mm upstream of front face
  fGun->SetParticleEnergy(fEkin * GeV);               // kinetic energy
  fGun->SetParticlePosition({fX * mm, fY * mm, z0 * mm});
  fGun->GeneratePrimaryVertex(evt);
}
