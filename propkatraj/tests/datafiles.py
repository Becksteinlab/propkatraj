"""
Input datafiles for tests
"""

__all__ = [
        "PSF_FRAME_ZERO_PKA",
        "PSF_FRAME_NINETY_PKA",
]

from pkg_resources import resource_filename

PSF_FRAME_ZERO_PKA = resource_filename(__name__,
                                       './datafiles/psf_ref_frame0.pka')
PSF_FRAME_NINETY_PKA = resource_filename(__name__,
                                         './datafiles/psf_ref_frame90.pka')
