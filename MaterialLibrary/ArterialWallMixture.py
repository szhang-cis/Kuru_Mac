#from __future__ import division
import numpy as np
from sys import exit
from numpy import einsum
from .MaterialBase import Material
from Kuru.Tensor import trace, Voigt #, makezero


class ArterialWallMixture(Material):
    """A incompressible anisotropic Fung Quadratic model with the energy given by:

        W(C) = rho_e*mu/2*J**(-2/3)*(C:I) + rho_m*1/2*k1m/k2m*(exp(k2m*(FN.FN-1)**2)-1)
                                          + rho_c*1/2*k1c/k2c*(exp(k2c*(FN.FN-1)**2)-1)
        U(J) = rho_e*kappa/2*(J-1)**2

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

        if self.ndim==3:
            self.H_VoigtSize = 6
        else:
            self.H_VoigtSize = 3

        # LOW LEVEL DISPATCHER
        self.has_low_level_dispatcher = True

        # FIELD VARIABLES AS GROWTH_&_REMODELING AND/OR DEPOSITION STRETCHES, ETC
        self.has_growth_remodeling = True
        self.has_state_variables = True

    def KineticMeasures(self,F, elem=0):
        # first three direction are use for rotation matrix (Normal, Tangential and Axial)
        anisotropic_orientations = self.anisotropic_orientations[elem,:,:]
        # State variables are filled first with remodeling then with densities and growth
        # remodeling (0-8 is elastin, 9 muscle, 10-13 collagen), densities (14-19) growth (20) 
        state_variables = self.StateVariables
        from Kuru.MaterialLibrary.LLDispatch._ArterialWallMixture_ import KineticMeasures
        return KineticMeasures(self,np.ascontiguousarray(F),
                np.ascontiguousarray(anisotropic_orientations),
                np.ascontiguousarray(state_variables))


    def Hessian(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        # Directional vector for element
        Normal = self.anisotropic_orientations[elem][0][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Tangential = self.anisotropic_orientations[elem][1][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Axial = self.anisotropic_orientations[elem][2][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Growth tensor (same for all constituent)
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.StateVariables[gcounter][20]*outerNormal + outerTangential
        # Remodeling tensor for elastin
        F_r = self.StateVariables[gcounter][0:9].reshape(3,3)
        F_r = np.dot(Rotation.T,np.dot(F_r,Rotation))
        # Elastin Growth and Remodeling tensor
        F_gr = np.dot(F_r,F_g)
        F_gr_inv = np.linalg.inv(F_gr)

        # ELASTIN
        mu = self.mu*self.StateVariables[gcounter][14]
        F_e = np.dot(F,F_gr_inv)
        J_e = np.linalg.det(F_e)
        b_e = np.dot(F_e,F_e.T)

        if self.ndim == 3:
            trb = trace(b_e)
        elif self.ndim == 2:
            trb = trace(b_e) + 1.

        H_Voigt = 2.*mu*(J_e**(-2./3.)/J)*(1./9.*trb*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b_e) - 1./3.*einsum('ij,kl',b_e,I) + \
                1./6.*trb*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) )

        if self.is_nearly_incompressible:
            H_Voigt += self.StateVariables[gcounter][14]*self.pressure*(J_e/J)*\
                (einsum('ij,kl',I,I)-(einsum('ik,jl',I,I)+einsum('il,jk',I,I)))
        else:
            kappa = self.kappa*self.StateVariables[gcounter][14]
            H_Voigt += kappa*(J_e/J)*((2.*J_e-1.)*einsum('ij,kl',I,I)-(J_e-1.)*(einsum('ik,jl',I,I)+einsum('il,jk',I,I)))

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            # Remodeling stretch in fibre direction
            lambda_r = self.StateVariables[gcounter][fibre_i+8]
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            if fibre_i is 1:
            	if innerFN_e >= 1.:
                	k1 = self.k1tm*self.StateVariables[gcounter][15]
                	k2 = self.k2tm
                else:
                	k1 = self.k1cm*self.StateVariables[gcounter][15]
                	k2 = self.k2cm
	            # Anisotropic Stiffness for this key
	        	H_Voigt += 4.*k1/J*(1.+2.*k2*(innerFN_e-1.)**2)*np.exp(k2*(innerFN_e-1.)**2)*einsum('ij,kl',outerFN_e,outerFN_e)
                # Active stress for SMC
                den0 = self.rho
                s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]
                stretch_m = self.maxi_active_stretch
                stretch_0 = self.zero_active_stretch
                stretch_a = self.active_stretch
            	H_Voigt += -2.*(s_act/(J*den0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*\
                        einsum('ij,kl',outerFN,outerFN)
            elif fibre_i is not 1:
            	if innerFN_e >=1.:
                	k1 = self.k1tc*self.StateVariables[gcounter][14+fibre_i]
                	k2 = self.k2tc
                else:
                	k1 = self.k1cc*self.StateVariables[gcounter][14+fibre_i]
                	k2 = self.k2cc	
	            # Anisotropic Stiffness for this key
	        	H_Voigt += 4.*k1/J*(1.+2.*k2*(innerFN_e-1.)**2)*np.exp(k2*(innerFN_e-1.)**2)*einsum('ij,kl',outerFN_e,outerFN_e)

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        # Directional vector for element
        Normal = self.anisotropic_orientations[elem][0][:,None]
        Normal = np.dot(I,Normal)[:,0]
        Tangential = self.anisotropic_orientations[elem][1][:,None]
        Tangential = np.dot(I,Tangential)[:,0]
        Axial = self.anisotropic_orientations[elem][2][:,None]
        Axial = np.dot(I,Axial)[:,0]
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]

        # Growth tensor (same for all constituent)
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.StateVariables[gcounter][20]*outerNormal + outerTangential
        # Remodeling tensor for elastin
        F_r = np.eye(3,3,dtype=np.float64)
        F_r[0,0] = self.StateVariables[gcounter][0]
        F_r[1,1] = self.StateVariables[gcounter][4]
        F_r[2,2] = self.StateVariables[gcounter][8]
        F_r = np.dot(Rotation.T,np.dot(F_r,Rotation))
        # Elastin Growth and Remodeling tensor
        F_gr = np.dot(F_r,F_g)
        F_gr_inv = np.linalg.inv(F_gr)

        # ELASTIN
        mu = self.mu*self.StateVariables[gcounter][14]
        F_e = np.dot(F,F_gr_inv)
        J_e = np.linalg.det(F_e)
        b_e = np.dot(F_e,F_e.T)

        if self.ndim == 3:
            trb = trace(b_e)
        elif self.ndim == 2:
            trb = trace(b_e) + 1.

        stress = mu*(J_e**(-2./3.)/J)*(b_e-1./3.*trb*I)

        if self.is_nearly_incompressible:
            stress += self.StateVariables[gcounter][14]*self.pressure*(J_e/J)*I
        else:
            kappa = self.kappa*self.StateVariables[gcounter][14]
            stress += kappa*(J_e-1.)*(J_e/J)*I

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            # Remodeling stretch in fibre direction
            lambda_r = self.StateVariables[gcounter][fibre_i+8]
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            if fibre_i is 1:
            	if innerFN_e >= 1.0:
                    k1 = self.k1tm*self.StateVariables[gcounter][15]
                    k2 = self.k2tm
            	else:
                	k1 = self.k1cm*self.StateVariables[gcounter][15]
                	k2 = self.k2cm
	            # Anisotropic Stiffness for this key
	        	stress += 2.*k1/J*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e
	            # Active stress for SMC
                den0 = self.rho
                s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]
                stretch_m = self.maxi_active_stretch
                stretch_0 = self.zero_active_stretch
                stretch_a = self.active_stretch
                stress += (s_act/(den0*innerFN*J))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN
            elif fibre_i is not 1:
            	if innerFN_e >= 1.0:
                    k1 = self.k1tc*self.StateVariables[gcounter][14+fibre_i]
                    k2 = self.k2tc
            	else:
                	k1 = self.k1cc*self.StateVariables[gcounter][14+fibre_i]
                	k2 = self.k2cc
	            # Anisotropic Stiffness for this key
	        	stress += 2.*k1/J*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e

        return stress

    def ConstituentStress(self,StrainTensors,elem=0,gcounter=0):

        density0 = self.rho
        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]

        #SMC AND COLLAGEN FIBRES
        fibre_stress = np.zeros((5),dtype=np.float64)
        softness = np.zeros((5),dtype=np.float64)
        for fibre_i in [1,2,3,4,5]:
            fraction = self.StateVariables[gcounter][14+fibre_i]/(J*density0)
            if fraction == 0.:
                continue
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            # Remodeling along the fibre
            lambda_r = self.StateVariables[gcounter][8+fibre_i]
            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            # Elastic deformation
            innerFN_e = innerFN/lambda_r**2
            if fibre_i is 1:
                if innerFN_e >= 1.0:
                    k1 = self.k1tm*self.StateVariables[gcounter][15]
                    k2 = self.k2tm
                else:
                    k1 = self.k1cm*self.StateVariables[gcounter][15]
                    k2 = self.k2cm
	            # Anisotropic Stress for this fibre
	            fibre_stress[0] = 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*innerFN_e/(J*fraction)
	            # Active stress for SMC
                s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]
                stretch_m = self.maxi_active_stretch
                stretch_0 = self.zero_active_stretch
                stretch_a = self.active_stretch
                fibre_stress[0] += (s_act/density0)*(1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)/(J*fraction)
            elif fibre_i is not 1:
                if innerFN_e >= 1.0:
                    k1 = self.k1tc*self.StateVariables[gcounter][14+fibre_i]
                    k2 = self.k2tc
                else:
                    k1 = self.k1cc*self.StateVariables[gcounter][14+fibre_i]
                    k2 = self.k2cc
	            # Anisotropic Stress for this fibre
	            fibre_stress[fibre_i-1] = 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*innerFN_e/(J*fraction)

            # Fibre softness for remodeling
            stiffness = (8.*k2*innerFN_e*(innerFN_e-1.)**2 + 8.*innerFN_e-4.)*k1*\
                np.sqrt(innerFN_e)*np.exp(k2*(innerFN_e-1.)**2)/(J*fraction)
            softness[fibre_i-1] = np.sqrt(innerFN)/(innerFN_e*stiffness)

        return fibre_stress,softness


    def LLConstituentStress(self,F,elem=0):

        density0 = self.rho
        gpoints = F.shape[0]
        fibre_stress = np.zeros((gpoints,5),dtype=np.float64)
        softness = np.zeros((gpoints,5),dtype=np.float64)
        I = np.eye(self.ndim,self.ndim,dtype=np.float64)
        for gcounter in range(gpoints):
            Fp = F[gcounter,:,:]
            J = np.linalg.det(Fp)

            #SMC AND COLLAGEN FIBRES
            for fibre_i in [1,2,3,4,5]:
                fraction = self.StateVariables[gcounter][14+fibre_i]/(J*density0)
                if fraction == 0.:
                    continue
                anisotropic_term = 0.
                # Fibre direction
                N = self.anisotropic_orientations[elem][fibre_i][:,None]
                N = np.dot(I,N)[:,0]
                FN = np.dot(Fp,N)
                # Remodeling along the fibre
                lambda_r = self.StateVariables[gcounter][8+fibre_i]
                # TOTAL deformation
                innerFN = einsum('i,i',FN,FN)
                # Elastic deformation
                innerFN_e = innerFN/lambda_r**2
                if fibre_i is 1:
                    if innerFN_e >= 1.0:
                        k1 = self.k1tm*self.StateVariables[gcounter][15]
                        k2 = self.k2tm
                    else:
                        k1 = self.k1cm*self.StateVariables[gcounter][15]
                        k2 = self.k2cm
	                # Anisotropic Stress for this fibre
	                fibre_stress[gcounter][0] = 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*innerFN_e/(J*fraction)
	                # Active stress for SMC
                    s_act = self.maxi_active_stress*self.StateVariables[gcounter][15]
                    stretch_m = self.maxi_active_stretch
                    stretch_0 = self.zero_active_stretch
                    stretch_a = self.active_stretch
                    fibre_stress[gcounter][0] += (s_act/density0)*(1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)/(J*fraction)
                elif fibre_i is not 1:
                    if innerFN_e >= 1.0:
                        k1 = self.k1tc*self.StateVariables[gcounter][14+fibre_i]
                        k2 = self.k2tc
                    else:
                        k1 = self.k1cc*self.StateVariables[gcounter][14+fibre_i]
                        k2 = self.k2cc
                    # Anisotropic Stress for this fibre
                    fibre_stress[gcounter][fibre_i-1] = 2.*k1*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*innerFN_e/(J*fraction)

                # Fibre softness for remodeling
                stiffness = (8.*k2*innerFN_e*(innerFN_e-1.)**2 + 8.*innerFN_e-4.)*k1*\
                    np.sqrt(innerFN_e)*np.exp(k2*(innerFN_e-1.)**2)/(J*fraction)
                softness[gcounter][fibre_i-1] = np.sqrt(innerFN)/(innerFN_e*stiffness)

        return fibre_stress,softness
