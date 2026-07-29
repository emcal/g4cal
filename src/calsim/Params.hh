#ifndef CALSIM_PARAMS_HH
#define CALSIM_PARAMS_HH

#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

// All lengths in mm, energies in GeV, times in ns (converted to Geant4 units at use site).
// Response-chain constants: harvested from halld_sim, see REFERENCES.md.
//
// Resolution terms follow the standard calorimeter model
//   sigma_E/E = sqrt(stochastic^2/E + noise^2/E^2 + constant^2)
// "measured_res_*" is the target (measured halld ECAL) resolution; "intrinsic_res_*" is
// the intrinsic Geant4 part; smearing adds the quadrature difference between the two.
struct Params {
  // geometry
  int crystals_nx = 3, crystals_ny = 3;   // grid size in crystals
  double crystal_side_mm   = 20.55;   // square side; halld ECAL crystal (hdds CrystalECAL_HDDS.xml:140)
  double crystal_length_mm = 200.0;
  double wrap_thickness_mm = 0.175;   // Tedlar per lateral face (hdds CrystalECAL_HDDS.xml:136-137)

  // gun
  std::string particle = "e-";
  double energy_min_gev = 5.0, energy_max_gev = 5.0;  // kinetic energy (G4 /gun/energy convention)
  std::string gun_mode = "central";                   // face | central | grid | point
  int gun_grid_n = 5;                                 // grid gun mode: gun_grid_n x gun_grid_n scan points
  double gun_x_mm = 0.0, gun_y_mm = 0.0;              // fixed impact for gun_mode=point
  int n_events = 100;
  long seed = 1;

  // physics: --hadronic off disables gamma/electro/muon-nuclear, the only
  // EM->hadronic doorway for gamma/e- primaries (EM-only showers, profile tables)
  bool hadronic_on = true;

  // output
  std::string out_dir = ".";
  std::string run_id  = "test";

  // response chain (see notes/DECISIONS.md; order: attenuate -> threshold -> scale -> smear)
  bool   attenuation_on    = true;
  double atten_length_mm   = 2000.0;  // hitECAL.c:27 ATTEN_LENGTH 200 cm
  double light_speed_mm_ns = 130.0;   // effective light speed for timing; hitECAL.c:28 C_EFFECTIVE 13 cm/ns
  bool   smearing_on       = true;
  double energy_scale      = 1.0962;  // ECALSmearer.cc:13, applied before smearing
  double measured_res_stochastic = 0.0308, measured_res_noise = 0.010, measured_res_constant = 0.007;  // ECALSmearer.cc:16-18
  double intrinsic_res_stochastic = 0.0171216, intrinsic_res_noise = 0.0155070;                        // ECALSmearer.cc:21-22
  double hit_threshold_gev = 0.005;   // hitECAL.c:33, on attenuated unscaled energy
  double store_min_gev     = 0.0001;  // store hit if edep_true above this (or e_vis > 0)
  bool   timing_on         = false;
  double time_sigma_ns     = 0.4;     // ECALSmearer.cc:26

  bool check_overlaps = false;

  double pitch_mm() const { return crystal_side_mm + 2.0 * wrap_thickness_mm; }

  static bool parseBool(const std::string& v) {
    return v == "1" || v == "on" || v == "true" || v == "yes";
  }

  static Params parse(int argc, char** argv) {
    Params p;
    for (int i = 1; i + 1 < argc; i += 2) {
      std::string k = argv[i], v = argv[i + 1];
      if      (k == "--crystals-nx")        p.crystals_nx = std::atoi(v.c_str());
      else if (k == "--crystals-ny")        p.crystals_ny = std::atoi(v.c_str());
      else if (k == "--crystal-side-mm")    p.crystal_side_mm   = std::atof(v.c_str());
      else if (k == "--crystal-length-mm")  p.crystal_length_mm = std::atof(v.c_str());
      else if (k == "--wrap-thickness-mm")  p.wrap_thickness_mm = std::atof(v.c_str());
      else if (k == "--particle")           p.particle = v;
      else if (k == "--energy-min-gev")     p.energy_min_gev = std::atof(v.c_str());
      else if (k == "--energy-max-gev")     p.energy_max_gev = std::atof(v.c_str());
      else if (k == "--energy-gev")         { p.energy_min_gev = p.energy_max_gev = std::atof(v.c_str()); }
      else if (k == "--gun-mode")           p.gun_mode = v;
      else if (k == "--gun-grid-n")         p.gun_grid_n = std::atoi(v.c_str());
      else if (k == "--gun-x-mm")           p.gun_x_mm = std::atof(v.c_str());
      else if (k == "--gun-y-mm")           p.gun_y_mm = std::atof(v.c_str());
      else if (k == "--events")             p.n_events = std::atoi(v.c_str());
      else if (k == "--seed")               p.seed = std::atol(v.c_str());
      else if (k == "--hadronic")           p.hadronic_on = parseBool(v);
      else if (k == "--out-dir")            p.out_dir = v;
      else if (k == "--run-id")             p.run_id = v;
      else if (k == "--attenuation")        p.attenuation_on = parseBool(v);
      else if (k == "--atten-length-mm")    p.atten_length_mm = std::atof(v.c_str());
      else if (k == "--light-speed-mm-ns")  p.light_speed_mm_ns = std::atof(v.c_str());
      else if (k == "--smearing")           p.smearing_on = parseBool(v);
      else if (k == "--energy-scale")       p.energy_scale = std::atof(v.c_str());
      else if (k == "--measured-res-stochastic")  p.measured_res_stochastic = std::atof(v.c_str());
      else if (k == "--measured-res-noise")       p.measured_res_noise = std::atof(v.c_str());
      else if (k == "--measured-res-constant")    p.measured_res_constant = std::atof(v.c_str());
      else if (k == "--intrinsic-res-stochastic") p.intrinsic_res_stochastic = std::atof(v.c_str());
      else if (k == "--intrinsic-res-noise")      p.intrinsic_res_noise = std::atof(v.c_str());
      else if (k == "--hit-threshold-gev")  p.hit_threshold_gev = std::atof(v.c_str());
      else if (k == "--store-min-gev")      p.store_min_gev = std::atof(v.c_str());
      else if (k == "--timing")             p.timing_on = parseBool(v);
      else if (k == "--time-sigma-ns")      p.time_sigma_ns = std::atof(v.c_str());
      else if (k == "--check-overlaps")     p.check_overlaps = parseBool(v);
      else { std::cerr << "calsim: unknown option " << k << "\n"; std::exit(2); }
    }
    return p;
  }

  std::string json() const {
    std::ostringstream o;
    o << "{\n"
      << "  \"run_id\": \"" << run_id << "\",\n"
      << "  \"crystals_nx\": " << crystals_nx << ", \"crystals_ny\": " << crystals_ny << ",\n"
      << "  \"crystal_side_mm\": " << crystal_side_mm
      << ", \"crystal_length_mm\": " << crystal_length_mm
      << ", \"wrap_thickness_mm\": " << wrap_thickness_mm
      << ", \"pitch_mm\": " << pitch_mm() << ",\n"
      << "  \"particle\": \"" << particle << "\", \"energy_min_gev\": " << energy_min_gev
      << ", \"energy_max_gev\": " << energy_max_gev << ",\n"
      << "  \"gun_mode\": \"" << gun_mode << "\", \"gun_grid_n\": " << gun_grid_n
      << ", \"gun_x_mm\": " << gun_x_mm << ", \"gun_y_mm\": " << gun_y_mm
      << ", \"direction\": [0, 0, 1],\n"
      << "  \"n_events\": " << n_events << ", \"seed\": " << seed << ",\n"
      << "  \"hadronic_on\": " << (hadronic_on ? "true" : "false") << ",\n"
      << "  \"attenuation_on\": " << (attenuation_on ? "true" : "false")
      << ", \"atten_length_mm\": " << atten_length_mm
      << ", \"light_speed_mm_ns\": " << light_speed_mm_ns << ",\n"
      << "  \"smearing_on\": " << (smearing_on ? "true" : "false")
      << ", \"energy_scale\": " << energy_scale << ",\n"
      << "  \"measured_res\": [" << measured_res_stochastic << ", " << measured_res_noise
      << ", " << measured_res_constant << "]"
      << ", \"intrinsic_res\": [" << intrinsic_res_stochastic << ", " << intrinsic_res_noise << "],\n"
      << "  \"hit_threshold_gev\": " << hit_threshold_gev
      << ", \"store_min_gev\": " << store_min_gev << ",\n"
      << "  \"timing_on\": " << (timing_on ? "true" : "false")
      << ", \"time_sigma_ns\": " << time_sigma_ns << "\n"
      << "}\n";
    return o.str();
  }
};

#endif
