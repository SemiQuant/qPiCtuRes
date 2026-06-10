"""qPiCtuRes - interactive qPCR analysis and plotting."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("qpictures")
except PackageNotFoundError:  # package is not installed
    __version__ = "0.0.0+local"

from .analysis import (
    AdjustName,
    PairedSummary,
    ReferenceMatch,
    StandardCurve,
    TestName,
    absolute_quantification,
    fit_standard_curve,
    paired_relative_quantification,
    relative_quantification,
)
from .io import read_thermo_well_results
from .plot import (
    absolute_quantity_barplot,
    paired_relative_quantity_barplot,
    relative_quantity_barplot,
    standard_curve_plot,
    standard_curve_with_quantities_figure,
)

__all__ = [
    "__version__",
    "read_thermo_well_results",
    "relative_quantification",
    "paired_relative_quantification",
    "fit_standard_curve",
    "absolute_quantification",
    "relative_quantity_barplot",
    "paired_relative_quantity_barplot",
    "standard_curve_plot",
    "absolute_quantity_barplot",
    "standard_curve_with_quantities_figure",
    "ReferenceMatch",
    "PairedSummary",
    "StandardCurve",
    "TestName",
    "AdjustName",
]
