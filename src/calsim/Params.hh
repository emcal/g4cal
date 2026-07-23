#ifndef CALSIM_PARAMS_HH
#define CALSIM_PARAMS_HH

#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

// All lengths in mm, energies in GeV, times in ns (converted to Geant4 units at use site).
// Response-chain constants: harvested from halld_sim, see REFERENCES.md.
struct Params {
  // geometry
  int nx = 3, ny = 3;
  double crystal_xy_mm = 20.55;   // halld ECAL crystal (hdds CrystalECAL_HDDS.xml:140)
  double crystal_z_mm  = 200.0;
  double wrap_mm       = 0.175;   // Tedlar per lateral face (hdds CrystalECAL_HDDS.xml:136-137)

  // gun
  std::string particle = "e-";
  double e_min_gev = 5.0, e_max_gev = 5.0;  // kinetic energy (G4 /gun/energy convention)
  std::string gun_mode = "central";         // face | central | grid | point
  int grid_n = 5;
  double gun_x_mm = 0.0, gun_y_mm = 0.0;    // fixed impact for gun_mode=point
  int n_events = 100;
  long seed = 1;

  // output
  std::string out_dir = ".";
  std::string run_id  = "test";

  // response chain (see notes/DECISIONS.md; order: attenuate -> threshold -> scale -> smear)
  bool   att_on     = true;
  double lambda_mm  = 2000.0;   // hitECAL.c:27 ATTEN_LENGTH 200 cm
  double veff_mm_ns = 130.0;    // hitECAL.c:28 C_EFFECTIVE 13 cm/ns
  bool   smear_on   = true;
  double en_scale   = 1.0962;   // ECALSmearer.cc:13, applied before smearing
  double meas_a = 0.0308, meas_b = 0.010, meas_c = 0.007;   // ECALSmearer.cc:16-18
  double intr_a = 0.0171216, intr_b = 0.0155070;            // ECALSmearer.cc:21-22
  double thresh_att_gev = 0.005;   // hitECAL.c:33, on attenuated unscaled energy
  double store_min_gev  = 0.0001;  // store hit if edep_true above this (or e_vis > 0)
  bool   time_on    = false;
  double t_sigma_ns = 0.4;         // ECALSmearer.cc:26

  bool check_overlaps = false;

  double pitch_mm() const { return crystal_xy_mm + 2.0 * wrap_mm; }

  static bool parseBool(const std::string& v) {
    return v == "1" || v == "on" || v == "true" || v == "yes";
  }

  static Params parse(int argc, char** argv) {
    Params p;
    for (int i = 1; i + 1 < argc; i += 2) {
      std::string k = argv[i], v = argv[i + 1];
      if      (k == "--nx") p.nx = std::atoi(v.c_str());
      else if (k == "--ny") p.ny = std::atoi(v.c_str());
      else if (k == "--crystal-xy") p.crystal_xy_mm = std::atof(v.c_str());
      else if (k == "--crystal-z")  p.crystal_z_mm  = std::atof(v.c_str());
      else if (k == "--wrap")       p.wrap_mm       = std::atof(v.c_str());
      else if (k == "--particle")   p.particle = v;
      else if (k == "--e-min")      p.e_min_gev = std::atof(v.c_str());
      else if (k == "--e-max")      p.e_max_gev = std::atof(v.c_str());
      else if (k == "--energy")     { p.e_min_gev = p.e_max_gev = std::atof(v.c_str()); }
      else if (k == "--gun-mode")   p.gun_mode = v;
      else if (k == "--grid-n")     p.grid_n = std::atoi(v.c_str());
      else if (k == "--gun-x")      p.gun_x_mm = std::atof(v.c_str());
      else if (k == "--gun-y")      p.gun_y_mm = std::atof(v.c_str());
      else if (k == "--events")     p.n_events = std::atoi(v.c_str());
      else if (k == "--seed")       p.seed = std::atol(v.c_str());
      else if (k == "--out-dir")    p.out_dir = v;
      else if (k == "--run-id")     p.run_id = v;
      else if (k == "--att")        p.att_on = parseBool(v);
      else if (k == "--lambda")     p.lambda_mm = std::atof(v.c_str());
      else if (k == "--veff")       p.veff_mm_ns = std::atof(v.c_str());
      else if (k == "--smear")      p.smear_on = parseBool(v);
      else if (k == "--en-scale")   p.en_scale = std::atof(v.c_str());
      else if (k == "--meas-a")     p.meas_a = std::atof(v.c_str());
      else if (k == "--meas-b")     p.meas_b = std::atof(v.c_str());
      else if (k == "--meas-c")     p.meas_c = std::atof(v.c_str());
      else if (k == "--intr-a")     p.intr_a = std::atof(v.c_str());
      else if (k == "--intr-b")     p.intr_b = std::atof(v.c_str());
      else if (k == "--thresh-att") p.thresh_att_gev = std::atof(v.c_str());
      else if (k == "--store-min")  p.store_min_gev = std::atof(v.c_str());
      else if (k == "--time")       p.time_on = parseBool(v);
      else if (k == "--t-sigma")    p.t_sigma_ns = std::atof(v.c_str());
      else if (k == "--check-overlaps") p.check_overlaps = parseBool(v);
      else { std::cerr << "calsim: unknown option " << k << "\n"; std::exit(2); }
    }
    return p;
  }

  std::string json() const {
    std::ostringstream o;
    o << "{\n"
      << "  \"run_id\": \"" << run_id << "\",\n"
      << "  \"nx\": " << nx << ", \"ny\": " << ny << ",\n"
      << "  \"crystal_xy_mm\": " << crystal_xy_mm << ", \"crystal_z_mm\": " << crystal_z_mm
      << ", \"wrap_mm\": " << wrap_mm << ", \"pitch_mm\": " << pitch_mm() << ",\n"
      << "  \"particle\": \"" << particle << "\", \"e_min_gev\": " << e_min_gev
      << ", \"e_max_gev\": " << e_max_gev << ",\n"
      << "  \"gun_mode\": \"" << gun_mode << "\", \"grid_n\": " << grid_n
      << ", \"gun_x_mm\": " << gun_x_mm << ", \"gun_y_mm\": " << gun_y_mm
      << ", \"direction\": [0, 0, 1],\n"
      << "  \"n_events\": " << n_events << ", \"seed\": " << seed << ",\n"
      << "  \"att_on\": " << (att_on ? "true" : "false")
      << ", \"lambda_mm\": " << lambda_mm << ", \"veff_mm_ns\": " << veff_mm_ns << ",\n"
      << "  \"smear_on\": " << (smear_on ? "true" : "false")
      << ", \"en_scale\": " << en_scale << ",\n"
      << "  \"meas\": [" << meas_a << ", " << meas_b << ", " << meas_c << "]"
      << ", \"intrinsic\": [" << intr_a << ", " << intr_b << "],\n"
      << "  \"thresh_att_gev\": " << thresh_att_gev
      << ", \"store_min_gev\": " << store_min_gev << ",\n"
      << "  \"time_on\": " << (time_on ? "true" : "false")
      << ", \"t_sigma_ns\": " << t_sigma_ns << "\n"
      << "}\n";
    return o.str();
  }
};

#endif
