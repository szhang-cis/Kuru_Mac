from __future__ import division
import numpy as np
from numpy import einsum
from .MaterialBase import Material
from Florence.Tensor import trace, Voigt, makezero


class ArterialWallMixture(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/4*c1m/c2m*(exp(c2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/4*c1c/c2c*(exp(c2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where could be possible
        more than one fiber family.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(ArterialWallMixture, self).__init__(mtype, ndim, **kwargs)
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

#    def KineticMeasures(self,F,ElectricFieldx=0, elem=0):
#        from Florence.MaterialLibrary.LLDispatch._AnisotropicFungQuadratic_ import KineticMeasures
#        return KineticMeasures(self, F, self.anisotropic_orientations[elem][:,None])

    def Hessian(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        # Directional vector for element
        Axial = self.anisotropic_orientations['co4'][elem][gcounter][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Tangential = self.anisotropic_orientations['co1'][elem][gcounter][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Normal = self.anisotropic_orientations['thick'][elem][gcounter][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Directional tensors
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal

        #ELASTIN
        kappa = self.kappa*self.mix_density['ela'][elem][gcounter]
        mu2D = self.mu2D[elem][gcounter]*self.mix_density['ela'][elem][gcounter]
        mu3D = self.mu3D*self.mix_density['ela'][elem][gcounter]
        Gh_ela = self.deposition_stretch['ela'][elem][gcounter]
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        J_ela = np.linalg.det(F_ela)
        b_ela = np.dot(F_ela,F_ela.T)
        b_tan = np.dot(F_ela,np.dot(outerTangential,F_ela.T))
        C_ela = np.dot(F_ela.T,F_ela)
        det_2D = np.linalg.det(np.dot(outerTangential,np.dot(C_ela,outerTangential))+outerNormal)

        if self.ndim == 3:
            trb_ela = trace(b_ela)
        elif self.ndim == 2:
            trb_ela = trace(b_ela) + 1.

        #SMC AND COLLAGEN FIBRES
        total_anisotropic = 0.
        for key in ['co1','co2','co3','co4','smc']:
            anisotropic_term = 0.
            # Fibre direction
            N = self.anisotropic_orientations[key][elem][gcounter][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if key is 'smc':
                FN = np.dot(self.deposition_stretch['smc'],FN)
                c1 = self.c1m*self.mix_density[key][elem][gcounter]
                c2 = self.c2m
            elif key is not 'smc':
                FN = np.dot(self.deposition_stretch['col'],FN)
                c1 = self.c1c*self.mix_density[key][elem][gcounter]
                c2 = self.c2c

            # Tensor of reference direction
            outerN = einsum('i,j',N,N)
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Anisotropic Stiffness for this key
            expo = np.exp(c2*(innerFN-1.)**2.)
            anisotropic_term = 2.*c1*(1.+2.*c2*(innerFN-1.)**2.)*expo*\
                    einsum('ij,kl',outerFN,outerFN)
            # Active stress for SMC
            if key is 'smc':
                den0 = 1050.0e-9
                s_act = 54.0e-3*self.mix_density[key][elem][gcounter]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                active_elasticity = -2.*(s_act/(den0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*\
                        einsum('ij,kl',outerFN,outerFN)
                anisotropic_term += active_elasticity

            # Addition to total anisotropic Stress
            total_anisotropic += anisotropic_term


        H_Voigt = 2.*mu3D*J_ela**(-2./3.)*(1./9.*trb_ela*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b_ela) - 1./3.*einsum('ij,kl',b_ela,I) + \
                1./6.*trb_ela*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) )/J + \
                2.*mu2D/det_2D*(einsum('ij,kl',outerTangential,outerTangential) + \
                1./2.*(einsum('il,jk',outerTangential,outerTangential) + \
                einsum('ik,jl',outerTangential,outerTangential)))/J + \
                total_anisotropic/J + \
                2.*kappa*J_ela*((2.*J_ela-1.)*einsum('ij,kl',I,I) - \
                (J_ela-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))/J

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        # Directional vector for element
        Axial = self.anisotropic_orientations['co4'][elem][gcounter][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Tangential = self.anisotropic_orientations['co1'][elem][gcounter][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Normal = self.anisotropic_orientations['thick'][elem][gcounter][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Total growth gradient deformation
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal

        #ELASTIN
        kappa = self.kappa*self.mix_density['ela'][elem][gcounter]
        mu2D = self.mu2D[elem][gcounter]*self.mix_density['ela'][elem][gcounter]
        mu3D = self.mu3D*self.mix_density['ela'][elem][gcounter]
        Gh_ela = self.deposition_stretch['ela'][elem][gcounter]
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        J_ela = np.linalg.det(F_ela)
        b_ela = np.dot(F_ela,F_ela.T)
        b_tan = np.dot(F_ela,np.dot(outerTangential,F_ela.T))
        C_ela = np.dot(F_ela.T,F_ela)
        det_2D = np.linalg.det(np.dot(outerTangential,np.dot(C_ela,outerTangential))+outerNormal)

        if self.ndim == 3:
            trb_ela = trace(b_ela)
        elif self.ndim == 2:
            trb_ela = trace(b_ela) + 1

        #SMC AND COLLAGEN FIBRES
        total_anisotropic = 0.
        for key in ['co1','co2','co3','co4','smc']:
            anisotropic_term = 0.
            # Fibre direction
            N = self.anisotropic_orientations[key][elem][gcounter][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if key is 'smc':
                FN = np.dot(self.deposition_stretch['smc'],FN)
                c1 = self.c1m*self.mix_density[key][elem][gcounter]
                c2 = self.c2m
            elif key is not 'smc':
                FN = np.dot(self.deposition_stretch['col'],FN)
                c1 = self.c1c*self.mix_density[key][elem][gcounter]
                c2 = self.c2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Anisotropic Stress for this fibre
            expo = np.exp(c2*(innerFN-1.)**2.)
            anisotropic_term = c1*(innerFN-1.)*expo*outerFN
            # Active stress for SMC
            if key is 'smc':
                den0 = 1050.0e-9
                s_act = 54.0e-3*self.mix_density[key][elem][gcounter]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                active_stress = (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*outerFN
                anisotropic_term += active_stress

            # Addition to total anisotropic Stress
            total_anisotropic += anisotropic_term


        stress = mu3D*J_ela**(-2./3.)*(b_ela - 1./3.*trb_ela*I)/J + \
            mu2D*(b_tan - outerTangential/det_2D)/J + \
            total_anisotropic/J + \
            2.*kappa*J_ela*(J_ela-1.)*I/J

        return stress
