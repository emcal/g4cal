"""Shared fitting helpers: iterative Gaussian peak fit and the resolution-curve fit
sigma_E/E = a/sqrt(E) (+) b/E (+) c (quadrature sum).
"""
import numpy as np
from scipy.optimize import curve_fit


def gauss(x, n, mu, sig):
    return n * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit_gauss_iterative(values, nsigma=2.0, iterations=4, bins=80):
    """Iteratively clipped Gaussian fit. Returns (mu, sig, mu_err, sig_err) or None."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 100:
        return None
    mu, sig = np.median(v), 1.4826 * np.median(np.abs(v - np.median(v)))  # robust start
    if sig <= 0:
        return None
    for _ in range(iterations):
        w = v[(v > mu - nsigma * sig) & (v < mu + nsigma * sig)]
        if len(w) < 100:
            return None
        counts, edges = np.histogram(w, bins=bins)
        centers = 0.5 * (edges[1:] + edges[:-1])
        try:
            popt, pcov = curve_fit(gauss, centers, counts,
                                   p0=[counts.max(), w.mean(), w.std()],
                                   sigma=np.sqrt(np.maximum(counts, 1)),
                                   absolute_sigma=True, maxfev=5000)
        except RuntimeError:
            return None
        mu, sig = popt[1], abs(popt[2])
    perr = np.sqrt(np.diag(pcov))
    return mu, sig, perr[1], perr[2]


def resolution_model(E, stochastic, noise, constant):
    """sigma_E/E as a function of E [GeV]: standard calorimeter three-term model."""
    return np.sqrt(stochastic**2 / E + noise**2 / E**2 + constant**2)


def res_terms(res_cfg):
    """[stochastic, noise, constant] for resolution_model from a config measured_res /
    intrinsic_res block; blocks without a constant term (intrinsic_res) get 0."""
    return [float(res_cfg.stochastic), float(res_cfg.noise), float(res_cfg.get("constant", 0.0))]


def fit_resolution(E, res, res_err):
    """Fit sigma/E vs E (quadrature model). Returns (popt, perr, chi2, ndf)."""
    E, res, res_err = map(np.asarray, (E, res, res_err))
    popt, pcov = curve_fit(resolution_model, E, res, p0=[0.03, 0.01, 0.007],
                           sigma=res_err, absolute_sigma=True,
                           bounds=([0, 0, 0], [1, 1, 1]), maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    chi2 = float(np.sum(((resolution_model(E, *popt) - res) / res_err) ** 2))
    return popt, perr, chi2, len(E) - 3
