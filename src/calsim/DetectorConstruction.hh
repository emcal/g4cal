#ifndef CALSIM_DETECTORCONSTRUCTION_HH
#define CALSIM_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "Params.hh"

class G4LogicalVolume;
class CrystalSD;

// n_x by n_y grid of PbWO4 crystals wrapped in Tedlar, centered at the origin,
// front face at z = -crystal_z/2, beam along +z. Copy number of the module
// placement = cell_id = iy*nx + ix (0-based).
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit DetectorConstruction(const Params& p) : fP(p) {}
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  CrystalSD* sd() const { return fSD; }

private:
  const Params& fP;
  G4LogicalVolume* fCrystalLV = nullptr;
  int fCopyDepth = 1;   // crystal sits inside the wrap module; module carries cell_id
  mutable CrystalSD* fSD = nullptr;
};

#endif
