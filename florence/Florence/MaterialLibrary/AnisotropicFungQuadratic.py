from __future__ import division
import numpy as np
from numpy import einsum
from .MaterialBase import Material
from Florence.Tensor import trace, Voigt, makezero


class AnisotropicFungQuadratic(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = u/2*J**(-2/3)*(C:I) + 1/4*c1/c2*(exp(c2*(FN.FN-1)**2)-1)
        U(J) = k*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where could be possible
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


    def Hessian(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

        mu = self.mu
        c1 = self.c1
        c2 = self.c2
        kappa = self.kappa

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        b = StrainTensors['b'][gcounter]
        F = StrainTensors['F'][gcounter]
        N1 = self.anisotropic_orientations['N1'][elem][:,None]
        FN1 = np.dot(F,N1)[:,0]
        innerFN1 = einsum('i,i',FN1,FN1)
        outerFN1 = einsum('i,j',FN1,FN1)
        N2 = self.anisotropic_orientations['N2'][elem][:,None]
        FN2 = np.dot(F,N2)[:,0]
        innerFN2 = einsum('i,i',FN2,FN2)
        outerFN2 = einsum('i,j',FN2,FN2)

        if self.ndim == 3:
            trb = trace(b)
        elif self.ndim == 2:
            trb = trace(b) + 1

        expo_1 = np.exp(c2*(innerFN1-1.)**2.)
        expo_2 = np.exp(c2*(innerFN2-1.)**2.)

        H_Voigt = 2.*mu*J**(-5./3.) * (1./9.*trb*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b) - 1./3.*einsum('ij,kl',b,I) + \
                1./6.*trb*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) ) + \
                2.*c1/J * (1.+2.*c2*(innerFN1-1.)**2.)*expo_1*einsum('ij,kl',outerFN1,outerFN1) + \
                2.*c1/J * (1.+2.*c2*(innerFN2-1.)**2.)*expo_2*einsum('ij,kl',outerFN2,outerFN2) + \
                2.*kappa*(2.*J-1.)*einsum('ij,kl',I,I) - 1./2.*kappa*(J-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I))

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        b = StrainTensors['b'][gcounter]
        F = StrainTensors['F'][gcounter]
        N1 = self.anisotropic_orientations['N1'][elem][:,None]
        FN1 = np.dot(F,N1)[:,0]
        innerFN1 = einsum('i,i',FN1,FN1)
        outerFN1 = einsum('i,j',FN1,FN1)
        N2 = self.anisotropic_orientations['N2'][elem][:,None]
        FN2 = np.dot(F,N2)[:,0]
        innerFN2 = einsum('i,i',FN2,FN2)
        outerFN2 = einsum('i,j',FN2,FN2)

        mu = self.mu
        c1 = self.c1
        c2 = self.c2
        kappa = self.kappa

        if self.ndim == 3:
            trb = trace(b)
        elif self.ndim == 2:
            trb = trace(b) + 1

        expo_1 = np.exp(c2*(innerFN1-1.)**2.)
        expo_2 = np.exp(c2*(innerFN2-1.)**2.)

        stress = mu*J**(-5./3.)*(b - 1./3.*trb*I) + \
            c1/J*(innerFN1-1.)*expo_1*outerFN1 + \
            c1/J*(innerFN2-1.)*expo_2*outerFN2 + \
            2.*kappa*(J-1.)*I

        return stress
