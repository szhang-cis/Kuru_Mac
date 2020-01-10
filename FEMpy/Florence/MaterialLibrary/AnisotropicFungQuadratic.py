from __future__ import division
import numpy as np
from numpy import einsum
from .MaterialBase import Material
from Florence.Tensor import trace, Voigt, makezero


class AnisotropicFungQuadratic(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = u/2*J**(-2/3)*(C:I-3) + \sum k1/(2*k2)*(exp(k2*(FN.FN-1)**2)-1)
        U(J) = k/2*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fibre-like part, where could be possible
        more than one fiber family.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(AnisotropicFungQuadratic, self).__init__(mtype, ndim, **kwargs)
        self.nvar = self.ndim
        self.is_transversely_isotropic = True
        self.energy_type = "internal_energy"
        self.nature = "nonlinear"
        self.fields = "mechanics"

        if self.ndim==3:
            self.H_VoigtSize = 6
        else:
            self.H_VoigtSize = 3

        # LOW LEVEL DISPATCHER
        #self.has_low_level_dispatcher = True
        self.has_low_level_dispatcher = False

#    def KineticMeasures(self,F,ElectricFieldx=0, elem=0):
#        from Florence.MaterialLibrary.LLDispatch._AnisotropicFungQuadratic_ import KineticMeasures
#        return KineticMeasures(self, F, self.anisotropic_orientations[elem][:,None])


    def Hessian(self,StrainTensors,growth_remodeling=None,elem=0,gcounter=0):

        kappa = self.kappa
        mu = self.mu
        k1 = self.k1
        k2 = self.k2

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        b = StrainTensors['b'][gcounter]
        F = StrainTensors['F'][gcounter]

        if self.ndim == 3:
            trb = trace(b)
        elif self.ndim == 2:
            trb = trace(b) + 1

        H_Voigt = 2.*mu*J**(-5./3.) * (1./9.*trb*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b) - 1./3.*einsum('ij,kl',b,I) + \
                1./6.*trb*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) ) + \
                kappa*((2.*J-1.)*einsum('ij,kl',I,I) - (J-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))

        # Anisotropic contibution
        for i_fibre in range(2):
            N = self.anisotropic_orientations[i_fibre][elem][:,None]
            FN = np.dot(F,N)[:,0]
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            expo = np.exp(k2*(innerFN-1.)**2)
            H_Voigt += 4.*k1/J*(1.+2.*k2*(innerFN-1.)**2)*expo*einsum('ij,kl',outerFN,outerFN)

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,growth_remodeling=None,elem=0,gcounter=0):

        kappa = self.kappa
        mu = self.mu
        k1 = self.k1
        k2 = self.k2

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        b = StrainTensors['b'][gcounter]
        F = StrainTensors['F'][gcounter]

        if self.ndim == 3:
            trb = trace(b)
        elif self.ndim == 2:
            trb = trace(b) + 1

        stress = mu*J**(-5./3.)*(b - 1./3.*trb*I) + kappa*(J-1.)*I

        # Anisotropic contibution
        for i_fibre in range(2):
            N = self.anisotropic_orientations[i_fibre][elem][:,None]
            FN = np.dot(F,N)[:,0]
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            expo = np.exp(k2*(innerFN-1.)**2)
            stress += 2.*k1/J*(innerFN-1.)*expo*outerFN

        return stress
