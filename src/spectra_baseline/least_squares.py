"""Implementations of least squares smoothing algorithims."""


from __future__ import annotations

import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import splu

__all__ = [
    "als",
]


def _als(
    spec: np.ndarray,
    DDT: sparse.csc_matrix,
    p: float,
    ratio: float,
    maxiter: int,
    w: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    L = spec.shape[0]
    # must be csc format for the SuperLU library used by scipy.linalg.splu
    W = sparse.diags(w, 0, shape=(L, L), format="csc")
    for _ in range(maxiter):
        Z = W + DDT
        LU = splu(Z)
        z = LU.solve(w * spec)
        w_last = w
        w = p * (spec > z) + (1 - p) * (spec < z)
        # use this direct setting instead of `setdiag` as this is much faster
        # ok in this specific case because we are garunteed that W is a diagonal matrix.
        W.data = w[None]
        if np.linalg.norm(w - w_last) < ratio:
            return z, w
    return z, w


def als(
    spectra: np.ndarray,
    lam: float = 1e5,
    p: float = 0.01,
    ratio: float = 0.1,
    maxiter: int = 20,
    w_init: np.ndarray = None,
) -> np.ndarray:
    """
    Run the als baselining algorithim on a single spectra.

    Parameters
    ----------
    spectra : array-like ([N], wns)
        The spectra to calculate the baseline of.
    lam : float, default 1e5
        TODO
    p : float, default 0.01
        TODO
    ratio : float, default 0.1
        The cutoff in the change of the norm of the baseline for convergence.
    maxiter: int, default 20
        The maximum number of iterations
    w_init : array-like (wns,)
        The initial w values to use. If None and *spectra* is multidimensional
        then w0 will be bootstrapped from the initial values.
        TODO: allow passing an array of w_init

    Notes
    -----
    Code is an optimized version of: https://stackoverflow.com/a/29185844/835607
    """
    spectra = np.atleast_2d(spectra)

    Z = np.zeros_like(spectra)
    # pre-compute DDT as it doesn't change when doing multiple specta
    L = spectra.shape[1]
    D = sparse.diags([1, -2, 1], offsets=[-1, 0, 1], shape=(L, L), format="csc")
    DDT = lam * D.dot(D.transpose())
    update_w_init = False
    if w_init is None:
        update_w_init = True
        w_init = np.ones(L)

    for i, spec in enumerate(spectra):
        Z[i], w = _als(spec, DDT, p, ratio, maxiter, w_init)
        if update_w_init:
            w_init = w
    # squeeze to return to original shape in case of only one spectra
    return (spectra - Z).squeeze()
