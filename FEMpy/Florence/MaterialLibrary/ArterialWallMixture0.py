from __future__ import division
import numpy as np
from sys import exit
from numpy import einsum
from .MaterialBase import Material
from Florence.Tensor import trace, Voigt, makezero


class ArterialWallMixture0(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/4*c1m/c2m*(exp(c2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/4*c1c/c2c*(exp(c2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where it could be possible
        more than one fiber family and it includes the deposition stretches.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(ArterialWallMixture0, self).__init__(mtype, ndim, **kwargs)
        self.nvar = self.ndim
        self.is_transversely_isotropic = True
        self.energy_type = "internal_energy"
        self.nature = "nonlinear"
        self.fields = "mechanics"
        self.has_deposition_stretch = True

        if self.ndim==3:
            self.H_VoigtSize = 6
        else:
            self.H_VoigtSize = 3

        # LOW LEVEL DISPATCHER
        #self.has_low_level_dispatcher = True
        self.has_low_level_dispatcher = False
        # GROWTH-REMODELING AND DEPOSITION STRETCH FOR TRY SENTENCE IN CHECKDATA
        #self.growth_remodeling = np.ones((8,12),dtype=np.float64)

    def Hessian(self,StrainTensors,growth_remodeling=None,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        if J<0.0:
            print(F)
            print(J)
            exit('Deformation Gradient is negative at element: '+str(elem)+' gauss: '+str(gcounter))
        # Directional vector for element
        Axial = self.anisotropic_orientations[5][elem][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Tangential = self.anisotropic_orientations[1][elem][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Normal = self.anisotropic_orientations[0][elem][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        #ELASTIN
        kappa = self.kappa*self.mixture_density[0]
        mu = self.mu*self.mixture_density[0]
        Gh_ela = self.deposition_stretch['Elastin']
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        J_ela = np.linalg.det(F_ela)
        b_ela = np.dot(F_ela,F_ela.T)

        if self.ndim == 3:
            trb_ela = trace(b_ela)
        elif self.ndim == 2:
            trb_ela = trace(b_ela) + 1.

        H_Voigt = 2.*mu*J_ela**(-2./3.)*(1./9.*trb_ela*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b_ela) - 1./3.*einsum('ij,kl',b_ela,I) + \
                1./6.*trb_ela*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) )/J + \
                kappa*J_ela*((2.*J_ela-1.)*einsum('ij,kl',I,I) - \
                (J_ela-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))/J

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[fibre_i][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.deposition_stretch['Muscle'],FN)
                k1 = self.k1m*self.mixture_density[1]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.deposition_stretch['Collagen'],FN)
                k1 = self.k1c*self.mixture_density[fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Anisotropic Stiffness for this key
            expo = np.exp(k2*(innerFN-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            H_Voigt += 4.*k1/J*(1.+2.*k2*(innerFN-1.)**2)*expo*einsum('ij,kl',outerFN,outerFN)
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.mixture_density[1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                H_Voigt += -2.*(s_act/(den0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*\
                        einsum('ij,kl',outerFN,outerFN)/J

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,growth_remodeling=None,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        if J<0.0:
            print(F)
            print(J)
            exit('Deformation Gradient is negative at element: '+str(elem)+' gauss: '+str(gcounter))
        # Directional vector for element
        Axial = self.anisotropic_orientations[5][elem][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Tangential = self.anisotropic_orientations[1][elem][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Normal = self.anisotropic_orientations[0][elem][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        #ELASTIN
        kappa = self.kappa*self.mixture_density[0]
        mu = self.mu*self.mixture_density[0]
        Gh_ela = self.deposition_stretch['Elastin']
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        J_ela = np.linalg.det(F_ela)
        b_ela = np.dot(F_ela,F_ela.T)

        if self.ndim == 3:
            trb_ela = trace(b_ela)
        elif self.ndim == 2:
            trb_ela = trace(b_ela) + 1.

        stress = mu*J_ela**(-2./3.)*(b_ela - (1./3.)*trb_ela*I)/J + kappa*J_ela*(J_ela-1.)*I/J

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[fibre_i][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.deposition_stretch['Muscle'],FN)
                k1 = self.k1m*self.mixture_density[1]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.deposition_stretch['Collagen'],FN)
                k1 = self.k1c*self.mixture_density[fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Anisotropic Stress for this fibre
            expo = np.exp(k2*(innerFN-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            stress += 2.*k1/J*(innerFN-1.)*expo*outerFN
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.mixture_density[1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                stress += (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN/J

        return stress

