r"""
Shared data
"""
import scanpy as sc
from pytest import fixture


@fixture
def adata() -> sc.AnnData:
    r"""
    Generate omics data
    """
    return sc.datasets.pbmc3k()
