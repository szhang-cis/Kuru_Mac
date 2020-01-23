from __future__ import division
import numpy as np
from sys import exit
from numpy import einsum
from .MaterialBase import Material
from Kuru.Tensor import trace, Voigt, makezero


class ArterialWallMixtureGR(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/4*c1m/c2m*(exp(c2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/4*c1c/c2c*(exp(c2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa*(J-1)**2

        This is a Nearly Incompressible NeoHookean and Fiber-like part, where it could be possible
        more than one fiber family and it includes the deposition stretches.

    """

    def __init__(self, ndim, **kwargs):
        mtype = type(self).__name__
        super(ArterialWallMixtureGR, self).__init__(mtype, ndim, **kwargs)
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

    def Hessian(self,StrainTensors,growth_remodeling=None,radial_stretch=None,mu2D=0.,elem=0,gcounter=0):

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

        # Growth Gradient Deformation
        if growth_remodeling is None:
            growth_remodeling = np.zeros((12),dtype=np.float64)
            growth_remodeling[0] = 241.5
            growth_remodeling[1] = 157.5
            growth_remodeling[2] = 65.1
            growth_remodeling[3] = 260.4
            growth_remodeling[4] = 260.4
            growth_remodeling[5] = 65.1
            growth_remodeling[6:12] = 1.0
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = growth_remodeling[11]*outerNormal + outerTangential
        F_g_inv = np.linalg.inv(F_g)

        #ELASTIN
        if radial_stretch is None:
            radial_stretch = 1./(self.deposition_stretch['Tangential']*self.deposition_stretch['Axial'])
        kappa = self.kappa*growth_remodeling[0]
        mu3D = self.mu3D*growth_remodeling[0]
        Gh_ela = np.eye(3,3,dtype=np.float64)
        Gh_ela[0,0] = radial_stretch
        Gh_ela[1,1] = self.deposition_stretch['Tangential']
        Gh_ela[2,2] = self.deposition_stretch['Axial']
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        F_ela_e = np.dot(F_ela,F_g_inv)
        J_ela_e = np.linalg.det(F_ela_e)
        b_ela_e = np.dot(F_ela_e,F_ela_e.T)
        Normal_gr = np.dot(F_g,Normal)
        Normal_gr = Normal_gr/np.linalg.norm(Normal_gr)
        outerNormal_gr = einsum('i,j',Normal_gr,Normal_gr)
        outerTangential_gr = I - outerNormal_gr
        b_ela_tan = np.dot(F_ela_e,np.dot(outerTangential_gr,F_ela_e.T))
        C_ela_e = np.dot(F_ela_e.T,F_ela_e)
        C_ela_tan = np.dot(outerTangential_gr,np.dot(C_ela_e,outerTangential_gr))
        det_2D = np.linalg.det(C_ela_tan + outerNormal_gr)

        if self.ndim == 3:
            trb_ela_e = trace(b_ela_e)
        elif self.ndim == 2:
            trb_ela_e = trace(b_ela_e) + 1.

        H_Voigt = 2.*mu3D*J_ela_e**(-2./3.)*(1./9.*trb_ela_e*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b_ela_e) - 1./3.*einsum('ij,kl',b_ela_e,I) + \
                1./6.*trb_ela_e*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) )/J + \
                mu2D*growth_remodeling[0]*(2.*einsum('ij,kl',outerTangential_gr,outerTangential_gr) + \
                einsum('il,jk',outerTangential_gr,outerTangential_gr) + \
                einsum('ik,jl',outerTangential_gr,outerTangential_gr))/(det_2D*J) +\
                kappa*J_ela_e*((2.*J_ela_e-1.)*einsum('ij,kl',I,I) - \
                (J_ela_e-1.)*(einsum('ik,jl',I,I) + einsum('il,jk',I,I)))/J

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[fibre_i][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.deposition_stretch['Muscle'],FN)
                k1 = self.k1m*growth_remodeling[1]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.deposition_stretch['Collagen'],FN)
                k1 = self.k1c*growth_remodeling[fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling stretch in fibre direction
            lambda_r = growth_remodeling[fibre_i+5]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            # Anisotropic Stiffness for this key
            expo = np.exp(k2*(innerFN_e-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            H_Voigt += 4.*k1/J*(1.+2.*k2*(innerFN_e-1.)**2)*expo*einsum('ij,kl',outerFN_e,outerFN_e)
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*growth_remodeling[1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                H_Voigt += -2.*(s_act/(den0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*\
                        einsum('ij,kl',outerFN,outerFN)/J

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,growth_remodeling=None,radial_stretch=None,mu2D=0.,elem=0,gcounter=0):

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

        # Growth Gradient Deformation
        if growth_remodeling is None:
            growth_remodeling = np.zeros((12),dtype=np.float64)
            growth_remodeling[0] = 241.5
            growth_remodeling[1] = 157.5
            growth_remodeling[2] = 65.1
            growth_remodeling[3] = 260.4
            growth_remodeling[4] = 260.4
            growth_remodeling[5] = 62.1
            growth_remodeling[6:12] = 1.0
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = growth_remodeling[11]*outerNormal + outerTangential
        F_g_inv = np.linalg.inv(F_g)

        #ELASTIN
        if radial_stretch is None:
            radial_stretch = 1./(self.deposition_stretch['Tangential']*self.deposition_stretch['Axial'])
        kappa = self.kappa*growth_remodeling[0]
        mu3D = self.mu3D*growth_remodeling[0]
        Gh_ela = np.eye(3,3,dtype=np.float64)
        Gh_ela[0,0] = radial_stretch
        Gh_ela[1,1] = self.deposition_stretch['Tangential']
        Gh_ela[2,2] = self.deposition_stretch['Axial']
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        F_ela_e = np.dot(F_ela,F_g_inv)
        J_ela_e = np.linalg.det(F_ela_e)
        b_ela_e = np.dot(F_ela_e,F_ela_e.T)
        Normal_gr = np.dot(F_g,Normal)
        Normal_gr = Normal_gr/np.linalg.norm(Normal_gr)
        outerNormal_gr = einsum('i,j',Normal_gr,Normal_gr)
        outerTangential_gr = I - outerNormal_gr
        b_ela_tan = np.dot(F_ela_e,np.dot(outerTangential_gr,F_ela_e.T))
        C_ela_e = np.dot(F_ela_e.T,F_ela_e)
        C_ela_tan = np.dot(outerTangential_gr,np.dot(C_ela_e,outerTangential_gr))
        det_2D = np.linalg.det(C_ela_tan + outerNormal_gr)

        if self.ndim == 3:
            trb_ela_e = trace(b_ela_e)
        elif self.ndim == 2:
            trb_ela_e = trace(b_ela_e) + 1.

        stress = mu3D*J_ela_e**(-2./3.)*(b_ela_e - (1./3.)*trb_ela_e*I)/J +\
                mu2D*growth_remodeling[0]*(b_ela_tan - outerTangential_gr/det_2D)/J +\
                kappa*J_ela_e*(J_ela_e-1.)*I/J

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[fibre_i][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.deposition_stretch['Muscle'],FN)
                k1 = self.k1m*growth_remodeling[1]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.deposition_stretch['Collagen'],FN)
                k1 = self.k1c*growth_remodeling[fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling stretch in fibre direction
            lambda_r = growth_remodeling[fibre_i+5]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            # Anisotropic Stress for this fibre
            expo = np.exp(k2*(innerFN_e-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            stress += 2.*k1/J*(innerFN_e-1.)*expo*outerFN_e
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*growth_remodeling[1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                stress += (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN/J

        return stress

    def ConstituentStress(self,StrainTensors,growth_remodeling=None,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]

        # Default Growth and Remodeling
        if growth_remodeling is None:
            growth_remodeling = np.zeros((12),dtype=np.float64)
            growth_remodeling[0] = 241.5
            growth_remodeling[1] = 157.5
            growth_remodeling[2] = 65.1
            growth_remodeling[3] = 260.4
            growth_remodeling[4] = 260.4
            growth_remodeling[5] = 62.1
            growth_remodeling[6:12] = 1.0
        #SMC AND COLLAGEN FIBRES
        fibre_stress = np.zeros((5),dtype=np.float64)
        softness = np.zeros((5),dtype=np.float64)
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[fibre_i][elem][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.deposition_stretch['Muscle'],FN)
                k1 = self.k1m*growth_remodeling[1]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.deposition_stretch['Collagen'],FN)
                k1 = self.k1c*growth_remodeling[fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = growth_remodeling[5+fibre_i]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # Elastic deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            # Anisotropic Stress for this fibre
            anisotropic_term = 2.*k1/J*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*growth_remodeling[1]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                anisotropic_term += (s_act/(den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN/J
            # Stress in the direction of the component
            outerN = einsum('i,j',N,N)
            fibre_stress[fibre_i-1] = np.tensordot(anisotropic_term,outerN)
            # Fibre softness for remodeling
            stiffness = 2.*(4.*k2*innerFN_e*(innerFN_e-1.)**2 + 4.*innerFN_e-2.)*k1*\
                np.sqrt(innerFN_e)*np.exp(k2*(innerFN_e-1.)**2)
            softness[fibre_i-1]=np.sqrt(innerFN)/(innerFN_e*stiffness)

        return fibre_stress,softness
