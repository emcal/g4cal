#include "DetectorConstruction.hh"
#include "CrystalSD.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();
  auto* air = nist->FindOrBuildMaterial("G4_AIR");
  auto* pwo = nist->FindOrBuildMaterial("G4_PbWO4");   // 8.28 g/cm3, matches hdds recipe

  // Tedlar C2H3F, 1.49 g/cm3 (hdds Material_HDDS.xml:697-702)
  auto* tedlar = new G4Material("Tedlar", 1.49 * g / cm3, 3);
  tedlar->AddElement(nist->FindOrBuildElement("C"), 2);
  tedlar->AddElement(nist->FindOrBuildElement("H"), 3);
  tedlar->AddElement(nist->FindOrBuildElement("F"), 1);

  const double pitch = fP.pitch_mm() * mm;
  const double cxy   = fP.crystal_side_mm * mm;
  const double cz    = fP.crystal_length_mm * mm;

  const double worldHx = 0.5 * fP.crystals_nx * pitch + 100. * mm;
  const double worldHy = 0.5 * fP.crystals_ny * pitch + 100. * mm;
  const double worldHz = 0.5 * cz + 400. * mm;

  auto* worldS  = new G4Box("world", worldHx, worldHy, worldHz);
  auto* worldLV = new G4LogicalVolume(worldS, air, "world");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0);

  const bool haveWrap = fP.wrap_thickness_mm > 0.;
  G4LogicalVolume* moduleLV = nullptr;
  if (haveWrap) {
    auto* moduleS = new G4Box("module", 0.5 * pitch, 0.5 * pitch, 0.5 * cz);
    moduleLV = new G4LogicalVolume(moduleS, tedlar, "module");
  }
  auto* crystalS = new G4Box("crystal", 0.5 * cxy, 0.5 * cxy, 0.5 * cz);
  fCrystalLV = new G4LogicalVolume(crystalS, pwo, "crystal");
  if (haveWrap) {
    new G4PVPlacement(nullptr, {}, fCrystalLV, "crystal", moduleLV, false, 0, fP.check_overlaps);
    fCopyDepth = 1;
  } else {
    moduleLV = fCrystalLV;   // place crystals directly; cell_id on the crystal itself
    fCopyDepth = 0;
  }

  for (int iy = 0; iy < fP.crystals_ny; ++iy) {
    for (int ix = 0; ix < fP.crystals_nx; ++ix) {
      const double x = (ix - 0.5 * (fP.crystals_nx - 1)) * pitch;
      const double y = (iy - 0.5 * (fP.crystals_ny - 1)) * pitch;
      const int cell = iy * fP.crystals_nx + ix;
      new G4PVPlacement(nullptr, {x, y, 0.}, moduleLV, "module", worldLV, false, cell,
                        fP.check_overlaps);
    }
  }
  return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
  fSD = new CrystalSD("crystalSD", fP, fCopyDepth);
  G4SDManager::GetSDMpointer()->AddNewDetector(fSD);
  SetSensitiveDetector(fCrystalLV, fSD);
}
