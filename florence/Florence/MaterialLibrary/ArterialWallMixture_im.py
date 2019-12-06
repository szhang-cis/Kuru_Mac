from __future__ import division
import numpy as np
from numpy import einsum
from .MaterialBase import Material
from Florence.Tensor import trace, Voigt, makezero


class ArterialWallMixture_im(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/4*c1m/c2m*(exp(c2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/4*c1c/c2c*(exp(c2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where could be possible
        more than one fiber family.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(ArterialWallMixture_im, self).__init__(mtype, ndim, **kwargs)
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
        Axial = self.anisotropic_orientations[5][elem][gcounter][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Tangential = self.anisotropic_orientations[2][elem][gcounter][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Normal = self.anisotropic_orientations[0][elem][gcounter][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Total growth gradient deformation
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        # Growth deformation gradient
        F_g = self.Growth[elem][gcounter]*outerNormal + outerTangential
        # Total inelastical deformation gradient (inverted)
        F_gr_inv = np.linalg.inv(F_g)

        #ELASTIN
        kappa = self.kappa*self.GrowthRemodeling[0][elem][gcounter]
        mu3D = self.mu3D*self.GrowthRemodeling[0][elem][gcounter]
        Gh_ela = self.deposition_stretch['ela'][elem][gcounter]
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        F_ela_e = np.dot(F_ela,F_gr_inv)
        J_ela_e = np.linalg.det(F_ela_e)
        b_ela_e = np.dot(F_ela_e,F_ela_e.T)
        F_push = np.dot(F_ela_e,F_gr_inv.T)
        Normal_push = np.dot(F_push,Normal)
        outerNormal_push = einsum('i,j',Normal_push,Normal_push)
        F_grNormal = np.dot(F_g,Normal)
        innerF_grNormal = einsum('i,i',F_grNormal,F_grNormal)

        if self.ndim == 3:
            trb_ela_e = trace(b_ela_e)
        elif self.ndim == 2:
            trb_ela_e = trace(b_ela_e) + 1.

        #SMC AND COLLAGEN FIBRES
        total_anisotropic = 0.
        for idx in [1,2,3,4,5]:
            anisotropic_term = 0.0
            # Fibre direction
            N = self.anisotropic_orientations[idx][elem][gcounter][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if idx is 1:
                FN = np.dot(self.deposition_stretch['smc'],FN)
                c1 = self.c1m*self.GrowthRemodeling[1][elem][gcounter]
                c2 = self.c2m
            elif idx is not 1:
                FN = np.dot(self.deposition_stretch['col'],FN)
                c1 = self.c1c*self.GrowthRemodeling[idx][elem][gcounter]
                c2 = self.c2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = self.GrowthRemodeling[idx+5][elem][gcounter]
            # Elastic deformation
            innerN_e = innerFN/lambda_r**2
            outerN_e = outerFN/lambda_r**2
            # Anisotropic Stiffness for this key
            expo = np.exp(c2*(innerN_e-1.)**2.)
            anisotropic_term = 2.*c1*(1.+2.*c2*(innerN_e-1.)**2.)*expo*einsum('ij,kl',outerN_e,outerN_e) #- \
                             #2.*c1*(1.-1./innerN_e)*expo*einsum('ij,kl',outerN_e,outerN_e)
            # Active stress for SMC
            if idx is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.GrowthRemodeling[1][elem][gcounter]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                active_stress = -2.*(s_act/(den0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*\
                        einsum('ij,kl',outerFN,outerFN)
                anisotropic_term += active_stress

            # Addition to total anisotropic Stress
            total_anisotropic += anisotropic_term


        H_Voigt = 2.*mu3D*J_ela_e**(-2./3.)*(1./9.*trb_ela_e*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b_ela_e) - 1./3.*einsum('ij,kl',b_ela_e,I) + \
                1./6.*trb_ela_e*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) )/J + \
                total_anisotropic/J + \
                2.*kappa*J_ela_e*((2.*J_ela_e-1.)*einsum('ij,kl',I,I) - \
                (J_ela_e-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))/J #- \
                #2.*mu3D*J_ela_e**(-2./3.)*innerF_grNormal*einsum('ij,kl',outerNormal_push,I)/J

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        self.TotalDeformation['F'][elem][gcounter][:,:] = F
        self.TotalDeformation['J'][elem][gcounter] = J
        # Directional vector for element
        Axial = self.anisotropic_orientations[5][elem][gcounter][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Tangential = self.anisotropic_orientations[2][elem][gcounter][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Normal = self.anisotropic_orientations[0][elem][gcounter][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Total growth gradient deformation
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.Growth[elem][gcounter]*outerNormal + outerTangential
        # Total inelastical deformation gradient (inverted)
        F_gr_inv = np.linalg.inv(F_g)

        #ELASTIN
        kappa = self.kappa*self.GrowthRemodeling[0][elem][gcounter]
        mu3D = self.mu3D*self.GrowthRemodeling[0][elem][gcounter]
        Gh_ela = self.deposition_stretch['ela'][elem][gcounter]
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        F_ela_e = np.dot(F_ela,F_gr_inv)
        J_ela_e = np.linalg.det(F_ela_e)
        b_ela_e = np.dot(F_ela_e,F_ela_e.T)

        if self.ndim == 3:
            trb_ela_e = trace(b_ela_e)
        elif self.ndim == 2:
            trb_ela_e = trace(b_ela_e) + 1.

        #SMC AND COLLAGEN FIBRES
        total_anisotropic = 0.0
        for idx in [1,2,3,4,5]:
            anisotropic_term = 0.0
            # Fibre direction
            N = self.anisotropic_orientations[idx][elem][gcounter][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if idx is 1:
                FN = np.dot(self.deposition_stretch['smc'],FN)
                c1 = self.c1m*self.GrowthRemodeling[1][elem][gcounter]
                c2 = self.c2m
            elif idx is not 1:
                FN = np.dot(self.deposition_stretch['col'],FN)
                c1 = self.c1c*self.GrowthRemodeling[idx][elem][gcounter]
                c2 = self.c2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = self.GrowthRemodeling[idx+5][elem][gcounter]
            # Elastic deformation
            innerN_e = innerFN/lambda_r**2
            outerN_e = outerFN/lambda_r**2
            # Anisotropic Stress for this fibre
            expo = np.exp(c2*(innerN_e-1.)**2.)
            anisotropic_term = c1*(innerN_e-1.)*expo*outerN_e
            # Active stress for SMC
            if idx is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.GrowthRemodeling[1][elem][gcounter]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                active_stress = (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*outerFN
                anisotropic_term += active_stress

            # Addition to total anisotropic Stress
            total_anisotropic += anisotropic_term


        stress = mu3D*J_ela_e**(-2./3.)*(b_ela_e - (1./3.)*trb_ela_e*I)/J + \
            total_anisotropic/J + \
            2.*kappa*J_ela_e*(J_ela_e-1.)*I/J

        return stress
