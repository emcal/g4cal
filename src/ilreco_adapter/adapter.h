#ifndef CAL_ILRECO_ADAPTER_H
#define CAL_ILRECO_ADAPTER_H

/*
 * Thin wrapper around ilreco.s context API for the g4cal tools.
 *
 * Cells are 0-based (ix, iy), as everywhere in g4cal; energies in GeV.
 * Cluster positions come back in 0-based cell units (x == ix means the
 * center of column ix). Conversion to mm is the caller's job; for an array
 * centered on the origin:
 *   x_mm = (x - 0.5 * (nx - 1)) * pitch_mm
 *
 * Single-threaded: one hidden config + workspace, which matches how the
 * calreco tools run. For multithreaded use call ilreco directly.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double e;      /* cluster energy, GeV (ilreco containment-corrected) */
  double x;      /* 0-based column coordinate, cell units */
  double y;      /* 0-based row coordinate, cell units */
  double chi2;   /* profile-fit chi2 / ndof */
  int    size;   /* number of hits in the cluster */
} CalCluster;

/* Create the reconstruction context for an n_cols x n_rows array using the
 * given shower-profile file. Call once (calling again replaces the context).
 * Returns 0 on success, -1 on failure (message on stderr). */
int cal_reco_init(const char* profile_path, int n_cols, int n_rows);

/* Array dimensions of the current context (0 before cal_reco_init). */
int cal_reco_ncol(void);
int cal_reco_nrow(void);

/* Run island reconstruction on one event. Inputs are parallel arrays of
 * nhits 0-based cells with energies in GeV. Clusters are written to
 * out[0 .. min(ret, max_out) - 1], energy-descending; returns the total
 * number found, or -1 on invalid input (cell outside the array, ...). */
int cal_reco_event(int nhits, const int* ix, const int* iy, const double* e_gev,
                   int max_out, CalCluster* out);

#ifdef __cplusplus
}
#endif

#endif
