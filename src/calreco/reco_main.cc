// calreco: run ilreco island reconstruction over a calsim hits.csv, write reco.csv.
// The array size is runtime configuration (--nx/--ny); the calreco_3x3 /
// calreco_10x10 / calreco_100x100 binaries are identical and kept only so the
// production scripts' paths stay valid.

#include "adapter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

namespace {

struct Args {
  std::string hits, out, profile;
  int nx = 0, ny = 0;
  double pitch_mm = 20.35;
  double cal = 1.0;        // global calibration constant (DECISIONS.md, stage d)
  double min_e_gev = 0.0;  // extra hit threshold on calibrated e_vis (block-threshold knob)
  std::string energy_col = "e_vis";   // or edep_true, for chain-OFF comparisons
  bool per_cluster = false;  // one row per cluster (multi-particle studies) vs leading only
};

Args parseArgs(int argc, char** argv) {
  Args a;
  for (int i = 1; i + 1 < argc; i += 2) {
    std::string k = argv[i], v = argv[i + 1];
    if      (k == "--hits") a.hits = v;
    else if (k == "--out") a.out = v;
    else if (k == "--profile") a.profile = v;
    else if (k == "--nx") a.nx = std::atoi(v.c_str());
    else if (k == "--ny") a.ny = std::atoi(v.c_str());
    else if (k == "--pitch") a.pitch_mm = std::atof(v.c_str());
    else if (k == "--cal") a.cal = std::atof(v.c_str());
    else if (k == "--min-e") a.min_e_gev = std::atof(v.c_str());
    else if (k == "--energy-col") a.energy_col = v;
    else if (k == "--per-cluster") a.per_cluster = (v == "1" || v == "on");
    else { std::fprintf(stderr, "calreco: unknown option %s\n", k.c_str()); std::exit(2); }
  }
  if (a.hits.empty() || a.out.empty() || a.profile.empty() || a.nx <= 0 || a.ny <= 0) {
    std::fprintf(stderr,
      "usage: calreco --hits hits.csv --out reco.csv --profile prof_pwo.dat "
      "--nx N --ny N [--pitch mm] [--cal c] [--min-e GeV] [--energy-col e_vis|edep_true]\n");
    std::exit(2);
  }
  return a;
}

struct Hit { int ix, iy; double e; };

void processEvent(long evt, std::vector<Hit>& hits, const Args& a, FILE* out) {
  static std::vector<int> vix, viy;
  static std::vector<double> ve;
  vix.clear(); viy.clear(); ve.clear();
  for (const auto& h : hits) {
    const double e = h.e * a.cal;
    if (e > a.min_e_gev && e > 0.) { vix.push_back(h.ix); viy.push_back(h.iy); ve.push_back(e); }
  }

  CalCluster cl[64];
  int ncl = 0;
  if (!ve.empty())
    ncl = cal_reco_event((int)ve.size(), vix.data(), viy.data(), ve.data(), 64, cl);
  if (ncl < 0) { std::fprintf(stderr, "calreco: bad input in event %ld\n", evt); std::exit(3); }

  const int nkeep = std::min(ncl, 64);
  if (a.per_cluster) {
    // one row per cluster, energy-descending icl index; zero-cluster events get icl=-1
    static std::vector<int> order;
    order.clear();
    for (int i = 0; i < nkeep; ++i) order.push_back(i);
    std::sort(order.begin(), order.end(),
              [&](int i, int j) { return cl[i].e > cl[j].e; });
    if (nkeep == 0) {
      std::fprintf(out, "%ld,-1,0,0,nan,nan,nan,0\n", evt);
    } else {
      for (size_t r = 0; r < order.size(); ++r) {
        const CalCluster& c = cl[order[r]];
        const double x = (c.x - 0.5 * (a.nx - 1)) * a.pitch_mm;
        const double y = (c.y - 0.5 * (a.ny - 1)) * a.pitch_mm;
        std::fprintf(out, "%ld,%zu,%d,%.6g,%.4f,%.4f,%.4g,%d\n",
                     evt, r, ncl, c.e, x, y, c.chi2, c.size);
      }
    }
    hits.clear();
    return;
  }
  const CalCluster* lead = nullptr;
  for (int i = 0; i < nkeep; ++i)
    if (!lead || cl[i].e > lead->e) lead = &cl[i];

  if (lead) {
    const double x = (lead->x - 0.5 * (a.nx - 1)) * a.pitch_mm;
    const double y = (lead->y - 0.5 * (a.ny - 1)) * a.pitch_mm;
    std::fprintf(out, "%ld,%d,%.6g,%.4f,%.4f,%.4g,%d\n",
                 evt, ncl, lead->e, x, y, lead->chi2, lead->size);
  } else {
    std::fprintf(out, "%ld,0,0,nan,nan,nan,0\n", evt);
  }
  hits.clear();
}

}  // namespace

int main(int argc, char** argv) {
  Args a = parseArgs(argc, argv);
  if (cal_reco_init(a.profile.c_str(), a.nx, a.ny) != 0) return 2;

  FILE* in = std::fopen(a.hits.c_str(), "r");
  if (!in) { std::fprintf(stderr, "calreco: cannot open %s\n", a.hits.c_str()); return 3; }
  FILE* out = std::fopen(a.out.c_str(), "w");
  if (!out) { std::fprintf(stderr, "calreco: cannot open %s\n", a.out.c_str()); return 3; }
  if (a.per_cluster) std::fprintf(out, "event,icl,n_clusters,e,x,y,chi2,size\n");
  else               std::fprintf(out, "event,n_clusters,e,x,y,chi2,size\n");

  char line[512];
  if (!std::fgets(line, sizeof line, in)) { std::fclose(in); std::fclose(out); return 0; }
  // header -> column index of the chosen energy column
  int ecol = -1;
  {
    int col = 0;
    for (char* tok = std::strtok(line, ",\n"); tok; tok = std::strtok(nullptr, ",\n"), ++col)
      if (a.energy_col == tok) ecol = col;
    if (ecol < 0) {
      std::fprintf(stderr, "calreco: column %s not in %s\n", a.energy_col.c_str(), a.hits.c_str());
      return 3;
    }
  }

  std::vector<Hit> hits;
  long cur = -1;
  long nev = 0;
  while (std::fgets(line, sizeof line, in)) {
    long evt; int ix, iy;
    double vals[8];
    int col = 0;
    char* tok = std::strtok(line, ",\n");
    evt = std::atol(tok);
    tok = std::strtok(nullptr, ",\n"); ix = std::atoi(tok);
    tok = std::strtok(nullptr, ",\n"); iy = std::atoi(tok);
    for (col = 3; tok && col < 8; ++col) {
      tok = std::strtok(nullptr, ",\n");
      if (tok) vals[col] = std::atof(tok);
    }
    if (evt != cur) {
      if (cur >= 0) { processEvent(cur, hits, a, out); ++nev; }
      cur = evt;
    }
    hits.push_back({ix, iy, vals[ecol]});
  }
  if (cur >= 0) { processEvent(cur, hits, a, out); ++nev; }

  std::fclose(in);
  std::fclose(out);
  std::fprintf(stderr, "calreco: %ld events with hits reconstructed\n", nev);
  return 0;
}
