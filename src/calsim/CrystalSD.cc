#include "CrystalSD.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

G4bool CrystalSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  const double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return false;

  const auto* pre  = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto* touch = pre->GetTouchable();
  const int cell = touch->GetCopyNumber(fDepth);

  const G4ThreeVector mid = 0.5 * (pre->GetPosition() + post->GetPosition());
  // local frame of the crystal volume: z in [-crystal_z/2, +crystal_z/2]
  const G4ThreeVector local = touch->GetHistory()->GetTopTransform().TransformPoint(mid);
  const double dist = 0.5 * fP.crystal_z_mm * mm - local.z();   // to back (readout) face

  double w = 1.0;
  if (fP.att_on) w = std::exp(-dist / (fP.lambda_mm * mm));

  auto& h = fHits[cell];
  h.edep += edep;
  h.eatt += edep * w;
  if (fP.time_on) {
    const double t = 0.5 * (pre->GetGlobalTime() + post->GetGlobalTime())
                     + dist / (fP.veff_mm_ns * mm / ns);
    h.twsum += t * edep * w;
  }
  return true;
}
