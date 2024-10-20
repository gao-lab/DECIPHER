r"""
`Spider` package for spatial omics data independent embedding
"""

import warnings
from importlib.metadata import version

ignore_warnings = [
    "The 'nopython' keyword",
    "is_categorical_dtype is deprecated",
    "Setting `dl_pin_memory_gpu_training`",
    "`use_gpu` is deprecated in v1.0",
    "The AnnData.concatenate",
    "SparseDataset is deprecated and will be removed",
    "UserWarning: No data for colormapping provided via",
    "FutureWarning: The default value of 'ignore' for the `na_action`",
    "Matplotlib created a",
    "The feature generate_power_seq is currently marked under review.",
    "The feature FeatureMapContrastiveTask is currently marked under review.",
    "The feature AmdimNCELoss is currently marked under review.",
]
for ignore_warning in ignore_warnings:
    warnings.filterwarnings("ignore", message=f".*{ignore_warning}.*")

from .cls import Spider  # noqa
from .utils import CFG  # noqa

name = "Spider"
__version__ = version(name)
__author__ = "Chen-Rui Xia"
