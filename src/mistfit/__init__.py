# src/mist_stellar_fitter/__init__.py
from .core import fit_stars_with_minimint

__all__ = ["mistfit"]

try:
    from importlib.metadata import version
    __version__ = version("mist-stellar-fitter")
except Exception:
    __version__ = "0+unknown"
