#ifndef CALSIM_PRIMARYGENERATORACTION_HH
#define CALSIM_PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "Params.hh"
#include <memory>

// Particle gun, normal incidence (+z), kinetic energy in GeV.
// Impact modes: face (uniform over full array face), central (uniform over the
// central crystal, cell (nx/2, ny/2)), grid (fixed grid_n x grid_n scan over the
// central crystal, cycled event by event).
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  explicit PrimaryGeneratorAction(const Params& p);
  void GeneratePrimaries(G4Event* evt) override;

  // truth of the last generated event
  int    pdg()  const { return fPdg; }
  double ekinGeV() const { return fEkin; }
  double xmm()  const { return fX; }
  double ymm()  const { return fY; }

private:
  const Params& fP;
  std::unique_ptr<G4ParticleGun> fGun;
  int fPdg = 0;
  double fEkin = 0, fX = 0, fY = 0;
  long fCount = 0;
};

#endif
