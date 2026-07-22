// Round-trip unit test for the ilreco adapter: an off-by-one in the cell
// mapping would silently shift every position by one crystal, so a single
// crystal fired at known (ix, iy) must reconstruct at exactly that cell.
// The grid size comes from TEST_GRID_N (one test binary per production preset).
// Usage: test_adapter <path-to-prof_pwo.dat>

#include "adapter.h"

#include <cmath>
#include <cstdio>

#ifndef TEST_GRID_N
#define TEST_GRID_N 10
#endif

static int n_failures = 0;

static void check(bool ok, const char* what) {
  std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
  if (!ok) ++n_failures;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: test_adapter prof_pwo.dat\n");
    return 2;
  }
  const int NX = TEST_GRID_N, NY = TEST_GRID_N;
  if (cal_reco_init(argv[1], NX, NY) != 0) return 2;
  std::printf("grid %dx%d\n", cal_reco_ncol(), cal_reco_nrow());

  // single-crystal round trips: corners, center
  const int probe_cells[][2] = {{0, 0}, {NX - 1, NY - 1}, {NX / 2, NY / 2},
                                {0, NY - 1}};
  for (const auto& cell : probe_cells) {
    const int ix = cell[0], iy = cell[1];
    const double energy = 1.0;
    CalCluster clusters[8];
    const int n_found = cal_reco_event(1, &ix, &iy, &energy, 8, clusters);
    char message[128];
    std::snprintf(message, sizeof message, "single hit at (%d,%d): 1 cluster", ix, iy);
    check(n_found == 1, message);
    if (n_found == 1) {
      std::snprintf(message, sizeof message,
                    "single hit at (%d,%d): position rounds to cell (got %.3f,%.3f)",
                    ix, iy, clusters[0].x, clusters[0].y);
      check(std::lround(clusters[0].x) == ix && std::lround(clusters[0].y) == iy,
            message);
      std::snprintf(message, sizeof message,
                    "single hit at (%d,%d): energy sane (got %.3f)", ix, iy,
                    clusters[0].e);
      check(clusters[0].e > 0.9 && clusters[0].e < 1.3, message);
    }
  }

  // two adjacent crystals, asymmetric split: one cluster between them,
  // nearer the hot one
  {
    const int center_x = NX / 2 - 1, center_y = NY / 2;
    const int ix[2] = {center_x, center_x + 1}, iy[2] = {center_y, center_y};
    const double energies[2] = {0.7, 0.3};
    CalCluster clusters[8];
    const int n_found = cal_reco_event(2, ix, iy, energies, 8, clusters);
    check(n_found >= 1, "two-hit split: cluster found");
    if (n_found >= 1) {
      check(clusters[0].x > center_x && clusters[0].x < center_x + 0.5,
            "two-hit split: x between cells, nearer the 70% one");
      check(std::lround(clusters[0].y) == center_y, "two-hit split: y at the row");
      check(clusters[0].e > 0.9 && clusters[0].e < 1.35, "two-hit split: energy sane");
    }
  }

  // empty event
  {
    CalCluster clusters[8];
    check(cal_reco_event(0, nullptr, nullptr, nullptr, 8, clusters) == 0,
          "empty event: 0 clusters");
  }

  // out-of-range cell rejected
  {
    const int ix = NX, iy = 0;
    const double energy = 1.0;
    CalCluster clusters[8];
    check(cal_reco_event(1, &ix, &iy, &energy, 8, clusters) == -1,
          "out-of-range cell rejected");
  }

  std::printf("%s (%d failures)\n", n_failures ? "TEST FAILED" : "ALL PASSED",
              n_failures);
  return n_failures ? 1 : 0;
}
