"""A collection of methodbaselining spectra."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("spectra-baseline")
except PackageNotFoundError:
    __version__ = "uninstalled"

__author__ = "Ian Hunt-Isaak"
__email__ = "ianhuntisaak@gmail.com"
