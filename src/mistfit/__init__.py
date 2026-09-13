"""mistfit -- MIST isochrone fitting with minimint + dynesty."""
from .core import fit_stars_with_minimint

__all__ = ["fit_stars_with_minimint"]

try:
    from importlib.metadata import version
    __version__ = version("mistfit")
except Exception:  # pragma: no cover - package not installed
    __version__ = "0+unknown"
