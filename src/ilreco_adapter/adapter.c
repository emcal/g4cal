#include "adapter.h"

#include <ilreco.h>
#include <stdio.h>

/* one hidden context: the calreco tools are single-threaded */
static ilreco_config* config = NULL;
static ilreco_workspace* workspace = NULL;
static int grid_cols = 0;
static int grid_rows = 0;

int cal_reco_init(const char* profile_path, int n_cols, int n_rows) {
  char error[256] = {0};
  if (config) {
    ilreco_workspace_destroy(workspace);
    ilreco_config_destroy(config);
    config = NULL;
    workspace = NULL;
    grid_cols = grid_rows = 0;
  }
  config = ilreco_config_create(n_cols, n_rows, profile_path, error, sizeof error);
  if (!config) {
    fprintf(stderr, "cal_reco_init: %s\n", error);
    return -1;
  }
  workspace = ilreco_workspace_create(config);
  if (!workspace) {
    fprintf(stderr, "cal_reco_init: workspace allocation failed\n");
    ilreco_config_destroy(config);
    config = NULL;
    return -1;
  }
  grid_cols = n_cols;
  grid_rows = n_rows;
  return 0;
}

int cal_reco_ncol(void) { return grid_cols; }
int cal_reco_nrow(void) { return grid_rows; }

int cal_reco_event(int nhits, const int* ix, const int* iy, const double* e_gev,
                   int max_out, CalCluster* out) {
  enum { MAX_HITS = 16384, MAX_CLUSTERS = 256 };
  static ilreco_hit hits[MAX_HITS];          /* off the stack; single-threaded */
  static ilreco_cluster clusters[MAX_CLUSTERS];

  if (!config || nhits < 0 || nhits > MAX_HITS || max_out < 0) return -1;

  for (int i = 0; i < nhits; ++i) {
    hits[i].col = ix[i];
    hits[i].row = iy[i];
    hits[i].e = e_gev[i];
  }

  /* ilreco validates the cells itself (out-of-array or masked -> -1) */
  const int n_found = ilreco_reconstruct(workspace, hits, nhits,
                                         clusters,
                                         max_out < MAX_CLUSTERS ? max_out
                                                                : MAX_CLUSTERS);
  if (n_found < 0) return -1;

  const int n_out = n_found < max_out ? n_found : max_out;
  for (int i = 0; i < n_out && i < MAX_CLUSTERS; ++i) {
    out[i].e = clusters[i].e;
    out[i].x = clusters[i].x;
    out[i].y = clusters[i].y;
    out[i].chi2 = clusters[i].chi2;
    out[i].size = clusters[i].size;
  }
  return n_found;
}
