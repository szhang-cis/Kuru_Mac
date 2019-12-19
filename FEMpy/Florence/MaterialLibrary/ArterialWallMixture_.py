from __future__ import division
import numpy as np
from sys import exit
from numpy import einsum
from .MaterialBase import Material
from Florence.Tensor import trace, Voigt, makezero


class ArterialWallMixture_(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/4*c1m/c2m*(exp(c2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/4*c1c/c2c*(exp(c2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where could be possible
        more than one fiber family.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(ArterialWallMixture_, self).__init__(mtype, ndim, **kwargs)
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
        self.growth_remodeling = np.ones((8,12),dtype=np.float64)

    def Hessian(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

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

        # Total growth gradient deformation
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.growth_remodeling[gcounter][11]*outerNormal + outerTangential
        # Total inelastical deformation gradient (inverted)
        F_gr_inv = np.linalg.inv(F_g)

        #ELASTIN
        kappa = self.kappa*self.growth_remodeling[gcounter][0]
        mu3D = self.mu3D*self.growth_remodeling[gcounter][0]
        Gh_ela = self.Deposition['Matrix']
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
        total_anisotropic = 0.
        for key in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[key][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            FN = np.dot(self.Deposition['Fibre'][key-1],FN)
            if key is 1:
                c1 = self.c1m*self.growth_remodeling[gcounter][1]
                c2 = self.c2m
            elif key is not 1:
                c1 = self.c1c*self.growth_remodeling[gcounter][key]
                c2 = self.c2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = self.growth_remodeling[gcounter][5+key]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # Elastic deformation
            innerN_e = innerFN/lambda_r**2
            outerN_e = outerFN/lambda_r**2
            # Anisotropic Stiffness for this key
            expo = np.exp(c2*(innerN_e-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            anisotropic_term = 2.*c1*(1.+2.*c2*(innerN_e-1.)**2.)*expo*einsum('ij,kl',outerN_e,outerN_e)
            # Active stress for SMC
            if key is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.growth_remodeling[gcounter][1]
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
                (J_ela_e-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))/J

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,ElectricFieldx=0,elem=0,gcounter=0):

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

        # Total growth gradient deformation
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.growth_remodeling[gcounter][11]*outerNormal + outerTangential
        # Total inelastical deformation gradient (inverted)
        F_gr_inv = np.linalg.inv(F_g)

        #ELASTIN
        kappa = self.kappa*self.growth_remodeling[gcounter][0]
        mu3D = self.mu3D*self.growth_remodeling[gcounter][0]
        Gh_ela = self.Deposition['Matrix']
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
        total_anisotropic = 0.
        for key in [1,2,3,4,5]:
            anisotropic_term = 0.
            # Fibre direction
            N = self.anisotropic_orientations[key][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            FN = np.dot(self.Deposition['Fibre'][key-1],FN)
            if key is 1:
                c1 = self.c1m*self.growth_remodeling[gcounter][1]
                c2 = self.c2m
            elif key is not 1:
                c1 = self.c1c*self.growth_remodeling[gcounter][key]
                c2 = self.c2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = self.growth_remodeling[gcounter][5+key]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # Elastic deformation
            innerN_e = innerFN/lambda_r**2
            outerN_e = outerFN/lambda_r**2
            # Anisotropic Stress for this fibre
            expo = np.exp(c2*(innerN_e-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            anisotropic_term = c1*(innerN_e-1.)*expo*outerN_e
            # Active stress for SMC
            if key is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.growth_remodeling[gcounter][1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                active_stress = (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN
                anisotropic_term += active_stress

            # Addition to total anisotropic Stress
            total_anisotropic += anisotropic_term


        stress = mu3D*J_ela_e**(-2./3.)*(b_ela_e - (1./3.)*trb_ela_e*I)/J + \
            total_anisotropic/J + \
            2.*kappa*J_ela_e*(J_ela_e-1.)*I/J

        return stress

    def ConstituentStress(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]

        #SMC AND COLLAGEN FIBRES
        fibre_stress = np.zeros((5),dtype=np.float64)
        softness = np.zeros((5),dtype=np.float64)
        for key in [1,2,3,4,5]:
            anisotropic_term = 0.
            # Fibre direction
            N = self.anisotropic_orientations[key][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            FN = np.dot(self.Deposition['Fibre'][key-1],FN)
            if key is 1:
                c1 = self.c1m*self.growth_remodeling[gcounter][1]
                c2 = self.c2m
            elif key is not 1:
                c1 = self.c1c*self.growth_remodeling[gcounter][key]
                c2 = self.c2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = self.growth_remodeling[gcounter][5+key]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # Elastic deformation
            innerN_e = innerFN/lambda_r**2
            outerN_e = outerFN/lambda_r**2
            # Anisotropic Stress for this fibre
            anisotropic_term = c1*(innerN_e-1.)*np.exp(c2*(innerN_e-1.)**2.)*outerN_e
            # Active stress for SMC
            if key is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.growth_remodeling[gcounter][1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                active_stress = (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*outerFN
                anisotropic_term += active_stress
            # Stress in the direction of the component
            outerN = einsum('i,j',N,N)
            fibre_stress[key-1] = np.tensordot(anisotropic_term,outerN)
            # Fibre softness for remodeling
            stiffness = (4.*c2*innerN_e*(innerN_e-1.)**2 + 4.*innerN_e-2.)*c1*\
                np.sqrt(innerN_e)*np.exp(c2*(innerN_e-1.)**2.)
            softness[key-1]=np.sqrt(innerFN)/(innerN_e*stiffness)

        return fibre_stress,softness
