"""
Input datafiles for tests
"""

__all__ = [
        "PSF_FRAME_ZERO_PKA",
        "PSF_FRAME_NINETY_PKA",
        "INDEXERR_FRAME14_GRO",
        "INDEXERR_FRAME14_XTC",
        "ATTERR_FRAME1_PDB",
        "ATTERR_FRAME1_XTC",

]

from importlib.resources import files

_DATA = files(__package__).joinpath("data")

PSF_FRAME_ZERO_PKA = str(_DATA.joinpath("psf_ref_frame0.pka"))
PSF_FRAME_NINETY_PKA = str(_DATA.joinpath("psf_ref_frame90.pka"))
# Issue #10
INDEXERR_FRAME14_GRO = str(_DATA.joinpath("indexerr_frame14.gro"))
INDEXERR_FRAME14_XTC = str(_DATA.joinpath("indexerr_frame14.xtc"))
ATTERR_FRAME1_PDB = str(_DATA.joinpath("atterr_frame1.pdb"))
ATTERR_FRAME1_XTC = str(_DATA.joinpath("atterr_frame1.xtc"))

del _DATA
