"""radiation_unit_conversion: convert between astronomical radiation units.

A Python implementation of the STScI radiation unit conversions
(https://www.stsci.edu/~strolger/docs/UNITS.txt).
"""
from importlib.metadata import PackageNotFoundError, version

from .units import (
    conversion,
    flambda2fnu,
    fnu2flambda,
    frequency2wavelength,
    wavelength2frequency,
)

try:
    __version__ = version("radiation_unit_conversion")
except PackageNotFoundError:  # package is not installed (e.g. running from source)
    __version__ = "0.0.0"

__all__ = [
    "conversion",
    "flambda2fnu",
    "fnu2flambda",
    "frequency2wavelength",
    "wavelength2frequency",
    "__version__",
]
