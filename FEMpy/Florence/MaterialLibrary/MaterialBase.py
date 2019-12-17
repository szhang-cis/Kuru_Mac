#from __future__ import print_function
#import numpy as np
#from Florence.Utils import insensitive
from warnings import warn

# BASE CLASS FOR ALL MATERIAL MODELS - SHOULD NOT BE USED DIRECTLY
class Material(object):
    """Base class for all material models"""

    def __init__(self, mtype, ndim, energy_type="internal_energy", density=None,
        is_compressible=True, is_incompressible=False, is_nearly_incompressible=False,
        is_nonisotropic=True,is_anisotropic=False,is_transversely_isotropic=False, anisotropic_orientations=None,
        **kwargs):


        # SAFETY CHECKS
        if not isinstance(mtype, str):
            raise TypeError("Type of material model should be given as a string")
        if not isinstance(energy_type, str):
            raise TypeError("Material energy can either be 'internal_energy' or 'enthalpy'")

        self.energy_type = energy_type

        # MATERIAL CONSTANTS
        self.rho = density

        # SET ALL THE OPTIONAL KEYWORD ARGUMENTS
        for i in kwargs.keys():
            if "__" not in i:
                setattr(self,i,kwargs[i])

        self.mtype = mtype
        self.ndim = ndim
        self.nvar = self.ndim

        self.H_Voigt = None

        if self.H_Voigt is not None:
            self.H_VoigtSize = self.H_Voigt.shape[0]

        self.is_compressible = is_compressible
        self.is_nearly_incompressible = is_nearly_incompressible
        self.is_incompressible = is_incompressible

        self.is_anisotropic = is_anisotropic
        self.is_transversely_isotropic = is_transversely_isotropic
        self.is_nonisotropic = is_nonisotropic
        self.anisotropic_orientations = anisotropic_orientations

        self.has_low_level_dispatcher = False

