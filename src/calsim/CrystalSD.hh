#ifndef CALSIM_CRYSTALSD_HH
#define CALSIM_CRYSTALSD_HH

#include "G4VSensitiveDetector.hh"
#include "Params.hh"
#include <map>

// Accumulates per-crystal energy per event. Attenuation (response stage b) is applied
// per step because it depends on the step's z inside the crystal: light travels to the
// readout at the BACK (+z) face, weight = exp(-(z_back - z_step)/lambda) (hitECAL.c:172).
class CrystalSD : public G4VSensitiveDetector {
public:
  struct CellAcc {
    double edep = 0;    // raw deposit (G4 units)
    double eatt = 0;    // attenuated deposit
    double twsum = 0;   // sum of t * attenuated deposit, for energy-weighted time
  };

  CrystalSD(const G4String& name, const Params& p, int copyDepth)
    : G4VSensitiveDetector(name), fP(p), fDepth(copyDepth) {}

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

  void BeginEvent() { fHits.clear(); }
  const std::map<int, CellAcc>& hits() const { return fHits; }

private:
  const Params& fP;
  int fDepth;   // touchable depth of the volume carrying cell_id as copy number
  std::map<int, CellAcc> fHits;
};

#endif
