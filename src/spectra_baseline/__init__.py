"""A collection of methodbaselining spectra."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("spectra-baseline")
except PackageNotFoundError:
    __version__ = "uninstalled"

__author__ = "Ian Hunt-Isaak"
__email__ = "ianhuntisaak@gmail.com"


# def baseline(spectra, lam=1e5, p=0.01, ratio: float = 0.1, maxiter: int = 20, w_init: np.ndarray = None):
#     """
#     Baseline spectra using the als method.

#     Parameters
#     ----------
#     spectra : ([N], wns)
#         1d or 2d array of spectra to be baselined.
#     w_init : array-like (wns,), optional 
#         The initial values for w0 for baselining. If not
#         provided and processing multiple spectra then w_init will be bootstrapped
#         using the initial computations.

#     returns
#     -------
#     baselined : array-like (N, wns)
#     """
#     spectra = np.atleast_2d(spectra)
#     Z = np.zeros_like(spectra)
#     for i, spec in enumerate(spectra):
#         Z[i], w_init = als(spec, lam, p, ratio, maxiter, w_init)
#     return spectra - Z
