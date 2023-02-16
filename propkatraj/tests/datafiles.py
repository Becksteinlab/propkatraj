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

from pkg_resources import resource_filename

PSF_FRAME_ZERO_PKA = resource_filename(__name__,
                                       './data/psf_ref_frame0.pka')
PSF_FRAME_NINETY_PKA = resource_filename(__name__,
                                         './data/psf_ref_frame90.pka')
# Issue #10
INDEXERR_FRAME14_GRO = resource_filename(__name__,
                                         './data/indexerr_frame14.gro')
INDEXERR_FRAME14_XTC = resource_filename(__name__,
                                         './data/indexerr_frame14.xtc')
ATTERR_FRAME1_PDB = resource_filename(__name__,
                                      './data/atterr_frame1.pdb')
ATTERR_FRAME1_XTC = resource_filename(__name__,
                                      './data/atterr_frame1.xtc')
