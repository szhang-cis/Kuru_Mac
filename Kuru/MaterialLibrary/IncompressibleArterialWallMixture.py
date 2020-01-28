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
        self.is_incompressible = True
        self.is_transversely_isotropic = True
        self.energy_type = "internal_energy"
        self.nature = "nonlinear"
        #self.fields = "mechanics"
        self.fields = "mechanics"

        if self.ndim==3:
            self.H_VoigtSize = 6
        else:
            self.H_VoigtSize = 3

        # LOW LEVEL DISPATCHER
        #self.has_low_level_dispatcher = True
        self.has_low_level_dispatcher = False

        # FIELD VARIABLES AS GROWTH_&_REMODELING AND/OR DEPOSITION STRETCHES, ETC
        self.has_field_variables = True
        #self.has_field_variables = False

    def MappingFieldVariables(self,mesh,function_space,elem=0):
        """
        Function to map the field variables to the Gauss points in the element.
        """
        Bases = function_space.Bases

        # Field variables at element nodes
        ElemFieldVariables = self.field_variables[mesh.elements[elem,:],:]

        # Field variables at gauss points
        self.FieldVariables = np.einsum('ij,ik->jk',Bases,ElemFieldVariables)

    def KineticMeasures(self,F, elem=0):
        # first three direction are use for rotation matrix (Normal, Tangential and Axial)
        anisotropic_orientations = self.anisotropic_orientations[elem,:,:]
        # Field variables are filled first with deposition_stretch then with densities, remodeling and growth
        # deposition (0-8 is elastin, 9 muscle, 10 collagen) then densities (11-16) remodeling (17-21) 
        # and growth (22)
        field_variables = self.FieldVariables
        from Kuru.MaterialLibrary.LLDispatch._ArterialWallMixture_ import KineticMeasures
        return KineticMeasures(self,np.ascontiguousarray(F),np.ascontiguousarray(anisotropic_orientations),
                np.ascontiguousarray(field_variables))


    def Hessian(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        if J<0.0:
            print(F)
            print(J)
            exit('Deformation Gradient is negative at element: '+str(elem)+' gauss: '+str(gcounter))
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

        # Growth tensor definition
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.FieldVariables[gcounter][22]*outerNormal + outerTangential
        F_g_inv = np.linalg.inv(F_g)

        #ELASTIN
        kappa = self.kappa*self.FieldVariables[gcounter][11]
        mu = self.mu*self.FieldVariables[gcounter][11]
        Gh_ela = np.eye(3,3,dtype=np.float64)
        Gh_ela[0,0] = self.FieldVariables[gcounter][0]
        Gh_ela[1,1] = self.FieldVariables[gcounter][4]
        Gh_ela[2,2] = self.FieldVariables[gcounter][8]
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        F_ela_e = np.dot(F_ela,F_g_inv)
        J_ela = np.linalg.det(F_ela_e)
        b_ela = np.dot(F_ela_e,F_ela_e.T)

        if self.ndim == 3:
            trb = trace(b_ela)
        elif self.ndim == 2:
            trb = trace(b_ela) + 1.

        H_Voigt = 2.*mu*(J_ela**(-2./3.)/J)*(1./9.*trb*einsum('ij,kl',I,I) - \
                1./3.*einsum('ij,kl',I,b_ela) - 1./3.*einsum('ij,kl',b_ela,I) + \
                1./6.*trb*(einsum('il,jk',I,I) + einsum('ik,jl',I,I)) ) + \
                kappa*(J_ela/J)*((2.*J-1.)*einsum('ij,kl',I,I)-(J-1.)*(einsum('ik,jl',I,I)+einsum('il,jk',I,I)))

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.FieldVariables[gcounter][9],FN)
                k1 = self.k1m*self.FieldVariables[gcounter][12]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.FieldVariables[gcounter][10],FN)
                k1 = self.k1c*self.FieldVariables[gcounter][11+fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling stretch in fibre direction
            lambda_r = self.FieldVariables[gcounter][fibre_i+16]
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
                s_act = 54.0e3*self.FieldVariables[gcounter][12]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                H_Voigt += -2.*(s_act/(J*den0*innerFN**2))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*\
                        einsum('ij,kl',outerFN,outerFN)

        H_Voigt = Voigt(H_Voigt ,1)

        return H_Voigt

    def CauchyStress(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]
        if J<0.0:
            print(F)
            print(J)
            exit('Deformation Gradient is negative at element: '+str(elem)+' gauss: '+str(gcounter))
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

        # Growth tensor definition
        outerNormal = einsum('i,j',Normal,Normal)
        outerTangential = I - outerNormal
        F_g = self.FieldVariables[gcounter][22]*outerNormal + outerTangential
        F_g_inv = np.linalg.inv(F_g)

        #ELASTIN
        kappa = self.kappa*self.FieldVariables[gcounter][11]
        mu = self.mu*self.FieldVariables[gcounter][11]
        Gh_ela = np.eye(3,3,dtype=np.float64)
        Gh_ela[0,0] = self.FieldVariables[gcounter][0]
        Gh_ela[1,1] = self.FieldVariables[gcounter][4]
        Gh_ela[2,2] = self.FieldVariables[gcounter][8]
        Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
        F_ela = np.dot(F,Gh_ela)
        F_ela_e = np.dot(F_ela,F_g_inv)
        J_ela = np.linalg.det(F_ela_e)
        b_ela = np.dot(F_ela_e,F_ela_e.T)

        if self.ndim == 3:
            trb = trace(b_ela)
        elif self.ndim == 2:
            trb = trace(b_ela) + 1.

        stress = mu*(J_ela**(-2./3.)/J)*(b_ela-1./3.*trb*I) + kappa*(J-1.)*(J_ela/J)*I

        #SMC AND COLLAGEN FIBRES
        for fibre_i in [1,2,3,4,5]:
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.FieldVariables[gcounter][9],FN)
                k1 = self.k1m*self.FieldVariables[gcounter][12]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.FieldVariables[gcounter][10],FN)
                k1 = self.k1c*self.FieldVariables[gcounter][11+fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling stretch in fibre direction
            lambda_r = self.FieldVariables[gcounter][fibre_i+16]
            # ELASTIC deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            # Anisotropic Stiffness for this key
            expo = np.exp(k2*(innerFN_e-1.)**2)
            if np.isnan(expo) or np.isinf(expo): exit('Fibre model is NaN or Inf: '+str(expo))
            stress += 2.*k1/J*(innerFN_e-1.)*expo*outerFN_e
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.FieldVariables[gcounter][12]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                stress += (s_act/(den0*innerFN*J))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN

        return stress

    def ConstituentStress(self,StrainTensors,elem=0,gcounter=0):

        I = StrainTensors['I']
        J = StrainTensors['J'][gcounter]
        F = StrainTensors['F'][gcounter]

        #SMC AND COLLAGEN FIBRES
        fibre_stress = np.zeros((5),dtype=np.float64)
        softness = np.zeros((5),dtype=np.float64)
        for fibre_i in [1,2,3,4,5]:
            anisotropic_term = 0.
            # Fibre direction
            N = self.anisotropic_orientations[elem][fibre_i][:,None]
            N = np.dot(I,N)[:,0]
            FN = np.dot(F,N)
            if fibre_i is 1:
                FN = np.dot(self.FieldVariables[gcounter][9],FN)
                k1 = self.k1m*self.FieldVariables[gcounter][12]
                k2 = self.k2m
            elif fibre_i is not 1:
                FN = np.dot(self.FieldVariables[gcounter][10],FN)
                k1 = self.k1c*self.FieldVariables[gcounter][11+fibre_i]
                k2 = self.k2c

            # TOTAL deformation
            innerFN = einsum('i,i',FN,FN)
            outerFN = einsum('i,j',FN,FN)
            # Remodeling along the fibre
            lambda_r = self.FieldVariables[gcounter][16+fibre_i]
            if np.isclose(lambda_r,0.0): exit('Remodeling in fibre is zero at element: '+str(elem))
            # Elastic deformation
            innerFN_e = innerFN/lambda_r**2
            outerFN_e = outerFN/lambda_r**2
            # Anisotropic Stress for this fibre
            anisotropic_term = 2.*k1/J*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e
            # Active stress for SMC
            if fibre_i is 1:
                den0 = 1050.0
                s_act = 54.0e3*self.FieldVariables[gcounter][12]
                stretch_m = 1.4
                stretch_a = 1.0
                stretch_0 = 0.8
                anisotropic_term += (s_act/(J*den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN
            # Stress in the direction of the component
            outerN = einsum('i,j',N,N)
            fibre_stress[fibre_i-1] = np.tensordot(anisotropic_term,outerN)
            # Fibre softness for remodeling
            stiffness = 2.*(4.*k2*innerFN_e*(innerFN_e-1.)**2 + 4.*innerFN_e-2.)*k1*\
                np.sqrt(innerFN_e)*np.exp(k2*(innerFN_e-1.)**2)
            softness[fibre_i-1]=np.sqrt(innerFN)/(innerFN_e*stiffness)

        return fibre_stress,softness


    def LLConstituentStress(self,F,elem=0):

        gpoints = F.shape[0]
        fibre_stress = np.zeros((gpoints,5),dtype=np.float64)
        softness = np.zeros((gpoints,5),dtype=np.float64)
        I = np.eye(self.ndim,self.ndim,dtype=np.float64)
        for gcounter in range(gpoints):
            Fp = F[gcounter,:,:]
            J = np.linalg.det(Fp)

            #SMC AND COLLAGEN FIBRES
            for fibre_i in [1,2,3,4,5]:
                anisotropic_term = 0.
                # Fibre direction
                N = self.anisotropic_orientations[elem][fibre_i][:,None]
                N = np.dot(I,N)[:,0]
                FN = np.dot(Fp,N)
                if fibre_i is 1:
                    FN = np.dot(self.FieldVariables[gcounter][9],FN)
                    k1 = self.k1m*self.FieldVariables[gcounter][12]
                    k2 = self.k2m
                elif fibre_i is not 1:
                    FN = np.dot(self.FieldVariables[gcounter][10],FN)
                    k1 = self.k1c*self.FieldVariables[gcounter][11+fibre_i]
                    k2 = self.k2c

                # TOTAL deformation
                innerFN = einsum('i,i',FN,FN)
                outerFN = einsum('i,j',FN,FN)
                # Remodeling along the fibre
                lambda_r = self.FieldVariables[gcounter][16+fibre_i]
                # Elastic deformation
                innerFN_e = innerFN/lambda_r**2
                outerFN_e = outerFN/lambda_r**2
                # Anisotropic Stress for this fibre
                anisotropic_term = 2.*k1/J*(innerFN_e-1.)*np.exp(k2*(innerFN_e-1.)**2)*outerFN_e
                # Active stress for SMC
                if fibre_i is 1:
                    den0 = 1050.0
                    s_act = 54.0e3*self.FieldVariables[gcounter][12]
                    stretch_m = 1.4
                    stretch_a = 1.0
                    stretch_0 = 0.8
                    anisotropic_term += (s_act/(J*den0*innerFN))*\
                        (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2)*outerFN
                # Stress in the direction of the component
                outerN = einsum('i,j',N,N)
                fibre_stress[gcounter][fibre_i-1] = np.tensordot(anisotropic_term,outerN)
                # Fibre softness for remodeling
                stiffness = 2.*(4.*k2*innerFN_e*(innerFN_e-1.)**2 + 4.*innerFN_e-2.)*k1*\
                    np.sqrt(innerFN_e)*np.exp(k2*(innerFN_e-1.)**2)
                softness[gcounter][fibre_i-1]=np.sqrt(innerFN)/(innerFN_e*stiffness)

        return fibre_stress,softness
