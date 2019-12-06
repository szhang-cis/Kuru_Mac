# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Florence
sys.path.append(os.path.join(os.path.expanduser("~"),"florence"))
#import Florence
from Florence import *

def Directions(mesh,LagrangeGaussCoords,fibre_direction):
    ngauss = LagrangeGaussCoords.shape[1]
    # Geometric definitions per element
    tangential = np.zeros((mesh.nelem,ngauss,ndim),dtype=np.float64)
    directrix = [0.,1.,0.]
    # Loop throught the element in the mesh
    for elem in range(mesh.nelem):
        for gcounter in range(ngauss):
            # Geometric definitions per element
            tangential[elem,gcounter,:] = np.cross(directrix,LagrangeGaussCoords[elem,gcounter,:])
            tangential[elem,gcounter,:] = tangential[elem,gcounter,:]/\
                    np.linalg.norm(tangential[elem,gcounter,:])
            fibre_direction['thick'][elem,gcounter,:] = np.cross(tangential[elem,gcounter,:],directrix)
            fibre_direction['thick'][elem,gcounter,:] = fibre_direction['thick'][elem,gcounter,:]/\
                    np.linalg.norm(fibre_direction['thick'][elem,gcounter,:])
            # Define the anisotropic orientations
            fibre_direction['smc'][elem,gcounter,:] = np.multiply(tangential[elem,gcounter,:],np.sin(np.pi/2))
            fibre_direction['co1'][elem,gcounter,:] = np.multiply(tangential[elem,gcounter,:],np.sin(np.pi/2))
            fibre_direction['co2'][elem,gcounter,:] = np.multiply(directrix,np.cos(np.pi/4)) + \
                    np.multiply(tangential[elem,gcounter,:],np.sin(np.pi/4))
            fibre_direction['co3'][elem,gcounter,:] = np.multiply(directrix,np.cos(np.pi/4)) - \
                    np.multiply(tangential[elem,gcounter,:],np.sin(np.pi/4))
            fibre_direction['co4'][elem,gcounter,:] = np.multiply(directrix,np.cos(0.))

def ComputeStress(material,elem,gcounter,stress,Stress_i,Softness):

    parameter2D = material.mu2D
    mu3D = material.mu3D
    c1m = material.c1m
    c2m = material.c2m
    c1c = material.c1c
    c2c = material.c2c
    kappa = material.kappa
    G_h = material.deposition_stretch
    density = material.mix_density
    Growth = material.Growth
    Remodeling = material.Remodeling
    TotalDeformation = material.TotalDeformation

    I = np.eye(3,3,dtype=np.float64)
    F = TotalDeformation['F'][elem][gcounter]
    J = TotalDeformation['J'][elem][gcounter]
    # Directional vector for element
    Axial = material.anisotropic_orientations['co4'][elem][gcounter][:,None]
    Axial = np.dot(I,Axial)[:,0]
    Tangential = material.anisotropic_orientations['co1'][elem][gcounter][:,None]
    Tangential = np.dot(I,Tangential)[:,0]
    Normal = material.anisotropic_orientations['thick'][elem][gcounter][:,None]
    Normal = np.dot(I,Normal)[:,0]
    Rotation = np.eye(3,3,dtype=np.float64)
    for i in range(3):
        Rotation[0,i] = Normal[i]
        Rotation[1,i] = Tangential[i]
        Rotation[2,i] = Axial[i]

    # Total growth gradient deformation
    outerNormal = einsum('i,j',Normal,Normal)
    F_g = Growth[elem][gcounter]*outerNormal + (I - outerNormal)
    # Total inelastical deformation gradient
    F_gr = F_g
    F_gr_inv = np.linalg.inv(F_gr)
    # Structural inelastic Tensors
    #Normal_gr = np.dot(F_gr,Normal)
    #innerNormal_gr = einsum('i,i',Normal_gr,Normal_gr)
    #Normal_gr = Normal_gr/np.sqrt(innerNormal_gr)
    #outerNormal_gr = einsum('i,j',Normal_gr,Normal_gr)
    outerTangential = I - outerNormal

    #ELASTIN
    Gh_ela = G_h['ela'][elem][gcounter]
    Gh_ela = np.dot(Rotation.T,np.dot(Gh_ela,Rotation))
    mu2D = parameter2D[elem][gcounter]*density['ela'][elem][gcounter]
    mu3D = mu3D*density['ela'][elem][gcounter]
    kappa = kappa*density['ela'][elem][gcounter]
    F_ela = np.dot(F,Gh_ela)
    F_ela_e = np.dot(F_ela,F_gr_inv)
    J_ela_e = np.linalg.det(F_ela_e)
    b_ela_e = np.dot(F_ela_e,F_ela_e.T)
    b_tan_e = np.dot(F_ela_e,np.dot(outerTangential,F_ela_e.T))
    C_ela_e = np.dot(F_ela_e.T,F_ela_e)
    det_2D = np.linalg.det(np.dot(outerTangential,np.dot(C_ela_e,outerTangential))+outerNormal)

    if material.ndim == 3:
        trb_ela_e = trace(b_ela_e)
    elif material.ndim == 2:
        trb_ela_e = trace(b_ela_e) + 1

    stress['ela'][elem][gcounter][:,:]=2.*kappa*J_ela_e*(J_ela_e-1.)*I/J+\
        mu3D*J_ela_e**(-2./3.)*(b_ela_e - 1./3.*trb_ela_e*I)/J + \
        mu2D*(b_tan_e - outerTangential/det_2D)/J

    #SMC AND COLLAGEN FIBRES
    for key in ['co1','co2','co3','co4','smc']:
        # Fibre direction
        N = material.anisotropic_orientations[key][elem][gcounter][:,None]
        N = np.dot(I,N)[:,0]
        FN = np.dot(F,N)
        if key is 'smc':
            FN = np.dot(G_h['smc'],FN)
            c1 = c1m*density[key][elem][gcounter]
            c2 = c2m
        elif key is not 'smc':
            FN = np.dot(G_h['col'],FN)
            c1 = c1c*density[key][elem][gcounter]
            c2 = c2c

        # Tensor of reference direction
        outerN = einsum('i,j',N,N)
        # TOTAL deformation
        innerFN = einsum('i,i',FN,FN)
        outerFN = einsum('i,j',FN,FN)
        # Remodeling tensor
        F_r = Remodeling[key][elem][gcounter]*outerN + \
                1./np.sqrt(Remodeling[key][elem][gcounter])*(I-outerN)
        # Inelastic deformation
        F_gr = np.dot(F_r,F_g)
        N_gr = np.dot(F_gr,N)
        innerN_gr = einsum('i,i',N_gr,N_gr)
        # Elastic deformation
        innerN_e = innerFN/innerN_gr
        outerN_e = outerFN/innerN_gr
        # Passive Stress for this fibre
        expo = np.exp(c2*(innerN_e-1.)**2.)
        fibre_stress = c1*(innerN_e-1.)*expo*outerN_e/J
        # Active stress for SMC
        if key is 'smc':
            den0_tot = 1050.0e-9
            s_act = 54.0e-3*density[key][elem][gcounter]
            stretch_m = 1.4
            stretch_a = 1.0
            stretch_0 = 0.8
            passive_stress = (s_act/(den0_tot*innerFN))*\
                    (1.-((stretch_m-stretch_a)/(stretch_m-stretch_0))**2.)*\
                    outerFN/J
            fibre_stress += passive_stress

        stress[key][elem][gcounter][:,:] = fibre_stress
        Stress_i[key][elem][gcounter] = np.tensordot(fibre_stress,outerN)
        # Fibre softness for remodeling
        Stiffness = innerN_e*(innerN_e-1.)**2.
        Stiffness = (4.*c2*Stiffness+4.*innerN_e-2.)*c1
        Stiffness = Stiffness*np.sqrt(innerN_e)
        Stiffness = Stiffness*expo
        Softness[key][elem][gcounter]=np.sqrt(innerFN)/(innerN_e*Stiffness)


def ElastinDepositionStretch(material,R_gauss,nelem,NoGauss):

    NewtonRaphson_failed = False
    pressure = 13.332e-3
    R = 10.0
    H = 1.41
    I = np.eye(3,3,dtype=np.float64)
    
    stress = {}
    Stress_i = {}
    Softness = {}
    stress['ela'] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
    for key in ['co1','co2','co3','co4','smc']:
        stress[key] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
        Stress_i[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
        Softness[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
    
    Residual = np.ones(2,dtype=np.float64)
    Stretch = np.ones(2,dtype=np.float64)
    print('Looking for Radial Deposition Stretch at Gauss points')
    for elem in range(nelem):
        if NewtonRaphson_failed:
            print('ElastinDepositionStretch stop at elem'+str(elem))
            break

        for gcounter in range(NoGauss):
            if NewtonRaphson_failed:
                print('ElastinDepositionStretch stop at elem'+str(elem)+' gcounter '+str(gcounter))
                break

            # Normal direction vector by element at Gauss point
            Normal = material.anisotropic_orientations['thick'][elem][gcounter][:,None]
            Normal = np.dot(I,Normal)[:,0]
            outerNormal = einsum('i,j',Normal,Normal)
            # Cuasi Newton-Raphson method
            Iter = 0
            Residual[0] = 1.0
            Residual[1] = 1.1
            Stretch[0] = 0.60
            Stretch[1] = 0.61
            norm_residual = np.absolute(Residual[0]/(pressure*(1.0-(R_gauss[elem][gcounter]-R)/H)))
            error = np.absolute(Stretch[0]-Stretch[1])/Stretch[0]
            while norm_residual>1.0e-5 or Iter==0:
                # values at trials
                for i in range(2):
                    material.deposition_stretch['ela'][elem][gcounter][0,0] = Stretch[i]
                    ComputeStress(material=material, elem=elem, gcounter=gcounter,
                            stress=stress, Stress_i=Stress_i, Softness=Softness)
                    # Total radial Cauchy stress
                    MixtureStress = 0.0
                    for key in ['ela','co1','co2','co3','co4','smc']:
                        MixtureStress += stress[key][elem][gcounter][:,:]

                    RadialStress = np.tensordot(MixtureStress,outerNormal)
                    Residual[i] = pressure*(1.0-(R_gauss[elem][gcounter]-R)/H)+RadialStress

                # secant
                Derivative = (Residual[1]-Residual[0])/(Stretch[1]-Stretch[0])
                Stretch[1] = Stretch[0]
                Stretch[0] = Stretch[0] - Residual[0]/Derivative
                norm_residual = np.absolute(Residual[0]/(pressure*(1.0-(R_gauss[elem][gcounter]-R)/H)))
                error = np.absolute(Stretch[0]-Stretch[1])/Stretch[0]
                Iter += 1
                material.deposition_stretch['ela'][elem][gcounter][0,0] = Stretch[0]

                # Check solution
                if error<1.0e-5:
                    material.deposition_stretch['ela'][elem][gcounter][0,0] = Stretch[0]
                    break

                # Check Residual and iterations
                if np.isnan(Residual[0]) or Residual[0]>1e06 or Iter>50:
                    NewtonRaphson_failed = True
                    print('ElastinDepositionStretch did not converge!')
                    break

    print('Radial Deposition Stretch calculated at Gauss points')

def CircumferentialStress(material,nelem,NoGauss):

    NewtonRaphson_failed = False
    pressure = 13.332e-3
    R = 10.0
    H = 1.41
    I = np.eye(3,3,dtype=np.float64)
    
    stress = {}
    Stress_i = {}
    Softness = {}
    stress['ela'] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
    for key in ['co1','co2','co3','co4','smc']:
        stress[key] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
        Stress_i[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
        Softness[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
    
    Residual = np.ones(2,dtype=np.float64)
    parameter = np.ones(2,dtype=np.float64)
    print('Looking for neo-Hookean parameter 2D at Gauss points')
    for elem in range(nelem):
        if NewtonRaphson_failed:
            print('CircumferentialStress stop at elem'+str(elem))
            break

        for gcounter in range(NoGauss):
            if NewtonRaphson_failed:
                print('CircumferentialStress stop at elem'+str(elem)+' gcounter '+str(gcounter))
                break

            # Normal direction vector by element at Gauss point
            Circ = material.anisotropic_orientations['co1'][elem][gcounter][:,None]
            Circ = np.dot(I,Circ)[:,0]
            outerCirc = einsum('i,j',Circ,Circ)
            # Cuasi Newton-Raphson method
            Iter = 0
            Residual[0] = 1.0
            Residual[1] = 1.1
            parameter[0] = 35.0e3
            parameter[1] = 15.0e3
            norm_residual = np.absolute(Residual[0]*H/(pressure*R))
            error = np.absolute(parameter[0]-parameter[1])/parameter[0]
            while norm_residual>1.0e-5 or Iter==0:
                # values at trials
                for i in range(2):
                    material.mu2D[elem][gcounter] = parameter[i]
                    ComputeStress(material=material, elem=elem, gcounter=gcounter,
                            stress=stress, Stress_i=Stress_i, Softness=Softness)
                    # Total radial Cauchy stress
                    MixtureStress = 0.0
                    for key in ['ela','co1','co2','co3','co4','smc']:
                        MixtureStress += stress[key][elem][gcounter][:,:]

                    CircStress = np.tensordot(MixtureStress,outerCirc)
                    Residual[i] = pressure*R/H - CircStress

                # secant
                Derivative = (Residual[1]-Residual[0])/(parameter[1]-parameter[0])
                parameter[1] = parameter[0]
                parameter[0] = parameter[0] - Residual[0]/Derivative
                norm_residual = np.absolute(Residual[0]*H/(pressure*R))
                error = np.absolute(parameter[0]-parameter[1])/parameter[0]
                Iter += 1
                material.mu2D[elem][gcounter] = parameter[0]

                # Check solution
                if error<1.0e-5:
                    material.mu2D[elem][gcounter] = parameter[0]
                    break

                # Check Residual and iterations
                if np.isnan(Residual[0]) or Residual[0]>1e06 or Iter>50:
                    NewtonRaphson_failed = True
                    print('CircumferentialStress did not converge!')
                    break

    print('neo-Hookean parameter 2D calculated at Gauss points')

def HomeostaticStress(material,Stress_H,nelem,NoGauss):
   
    stress = {}
    Stress_i = {}
    Softness = {}
    stress['ela'] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
    for key in ['co1','co2','co3','co4','smc']:
        stress[key] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
        Stress_i[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
        Softness[key] = np.zeros((nelem,NoGauss),dtype=np.float64)

    for elem in range(nelem):
        for gcounter in range(NoGauss):
            ComputeStress(material=material, elem=elem, gcounter=gcounter,
                    stress=stress, Stress_i=Stress_i, Softness=Softness)

    #SMC AND COLLAGEN FIBRES
    for key in ['co1','co2','co3','co4','smc']:
        # Stress for this key or fibre
        Stress_H[key] = Stress_i[key]


def GrowthRemodeling(time,material,Stress_H,density_rate,Remodeling_rate,
        Y_square,den0_tot,den0_e,nelem,NoGauss):
    #**** compute rates of G&R at t_n **** \dot{F}_{gr}(t_n)
    density = material.mix_density
    Remodeling = material.Remodeling
    
    stress = {}
    Stress_i = {}
    Softness = {}
    density_rate = {}
    Remodeling_rate = {}
    stress['ela'] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
    for key in ['co1','co2','co3','co4','smc']:
        stress[key] = np.zeros((nelem,NoGauss,3,3),dtype=np.float64)
        Stress_i[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
        Softness[key] = np.zeros((nelem,NoGauss),dtype=np.float64)
        density_rate[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
        Remodeling_rate[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)

    for elem in range(nelem):
        for gcounter in range(NoGauss):
            ComputeStress(material=material, elem=elem, gcounter=gcounter,
                    stress=stress, Stress_i=Stress_i, Softness=Softness)
    
    # Elastin degradation [days]
    L_dam = 10.0
    t_dam = 40.0
    D_max = 0.5
    T_ela = 101.0*365.25
    # elastin density rate
    exponential = np.exp(-0.5*Y_square/L_dam**2 - time/t_dam)
    density_rate['ela']= -density['ela']/T_ela - (D_max/t_dam)*np.multiply(exponential,den0_e)
    # Fribres turnover and remodeling
    turnover = 101.0
    gain = 0.05/turnover
    #SMC AND COLLAGEN FIBRES
    for key in ['co1','co2','co3','co4','smc']:
        # fibres density rates
        DeltaStress = Stress_i[key] - Stress_H[key]
        density_rate[key] = gain*np.multiply(density[key],np.divide(DeltaStress,Stress_H[key]))
        # fibres remodeling rates
        lambda_r_dot = np.divide(density_rate[key],density[key]) + 1.0/turnover
        lambda_r_dot = np.multiply(lambda_r_dot,DeltaStress)
        Remodeling_rate[key] = np.multiply(lambda_r_dot,Softness[key])

    den_tot = np.zeros((nelem,NoGauss),dtype=np.float64)
    density['ela'] = density['ela'] + Delta_t*density_rate['ela']
    den_tot += density['ela']
    # SMC and Collagen turnover [months]
    for key in ['co1','co2','co3','co4','smc']:
    # SMC and Collagen turnover [months]
        density[key] = density[key] + Delta_t*density_rate[key]
        den_tot += density[key]
    # SMC and Collagen remodeling [months]
        Remodeling[key] = Remodeling[key] + Delta_t*Remodeling_rate[key]

    # Total Growth magnitud
    Growth = np.divide(den_tot,den0_tot)

    # Update
    material.mix_density = density
    material.Remodeling = Remodeling
    material.Growth = Growth

def homogenized_CMT():
    """A hyperelastic explicit dynamics example using Mooney Rivlin model
        of a cylinder column under compression with cubic (p=3) hexahedral elements
    """
    ProblemPath = PWD(__file__)
    mesh_file = ProblemPath + '/Quarter_Cylinder.msh'

#===============  MESH PROCESING  ==========================
    # Build mesh with Florence tools from GMSH mesh
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex",read_surface_info=True)
    ndim = mesh.InferSpatialDimension()
    #Boolean arrays for boundary condition in Dirichlet
    BottomSurface = np.zeros(mesh.nnode,dtype=bool)
    TopSurface = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_Z = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_X = np.zeros(mesh.nnode,dtype=bool)
    #Boolean array for boundary condition in Neumann
    InnerSurface = np.zeros(mesh.faces.shape[0],dtype=bool)
    InnerFaces = []
    for idface in range(mesh.faces.shape[0]):
        if mesh.face_to_surface[idface] == 13:
            for i in range(mesh.faces.shape[1]):
                BottomSurface[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 21:
            for i in range(mesh.faces.shape[1]):
                TopSurface[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 4:
            for i in range(mesh.faces.shape[1]):
                Symmetry_Z[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 26:
            for i in range(mesh.faces.shape[1]):
                Symmetry_X[mesh.faces[idface,i]] = True
        elif mesh.face_to_surface[idface] == 25:
            InnerFaces.append([int(i) for i in mesh.faces[idface,:]])
            InnerSurface[idface] = True
    
    DirichletBoundary = {}
    DirichletBoundary['Bottom'] = BottomSurface
    DirichletBoundary['Top'] = TopSurface
    DirichletBoundary['SymmetryZ'] = Symmetry_Z
    DirichletBoundary['SymmetryX'] = Symmetry_X

    InnerFaces = np.array(InnerFaces,copy=True)
    NeumannBoundary = {}
    NeumannBoundary['InnerLogic'] = InnerSurface
    NeumannBoundary['InnerFaces'] = InnerFaces

#==================  FORMULATION  =========================
    formulation = DisplacementFormulation(mesh)
    NoGauss = formulation.function_spaces[0].AllGauss.shape[0]
    Bases = formulation.function_spaces[0].Bases

    LagrangeGaussCoords = np.zeros((mesh.nelem,NoGauss,3),dtype=np.float64)
    for elem in range(mesh.nelem):
        LagrangeCoords = mesh.points[mesh.elements[elem,:],:]
        ElemLagrangeGaussCoords = np.einsum('ij,ik->jk',Bases,LagrangeCoords)
        for i in range(3):
            LagrangeGaussCoords[elem,:,i] = ElemLagrangeGaussCoords[:,i]

    R_gauss = np.square(LagrangeGaussCoords[:,:,0]) + \
            np.square(LagrangeGaussCoords[:,:,2])
    R_gauss = np.sqrt(R_gauss)

#===============  MATERIAL DEFINITION  ====================
    # Kinematic Definiton
    TotalDeformation = {}
    TotalDeformation['F'] = np.zeros((mesh.nelem,NoGauss,3,3),dtype=np.float64)
    TotalDeformation['F'][:,:,0,0] = 1.0
    TotalDeformation['F'][:,:,1,1] = 1.0
    TotalDeformation['F'][:,:,2,2] = 1.0
    TotalDeformation['J'] = np.ones((mesh.nelem,NoGauss),dtype=np.float64)
    # Dictionary for fibres orientations
    fibre_direction = {}
    for key in ['co1','co2','co3','co4','smc','thick']:
        fibre_direction[key] = np.zeros((mesh.nelem,NoGauss,ndim),dtype=np.float64)

    # Dictionary for Density mixture and Inelastic gradient
    total_density = 1050.0e-9
    density = {}
    for key in ['ela','co1','co2','co3','co4','smc']:
        density[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)

    # Dictionary for Directional Stresses and Remodeling directional magnitude
    Stress_H = {}
    Remodeling = {}
    for key in ['co1','co2','co3','co4','smc']:
        Stress_H[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
        Remodeling[key] = np.ones((mesh.nelem,NoGauss),dtype=np.float64)

    # Densities per constituent
    density['ela'][:,:] = 0.23*total_density
    density['smc'][:,:] = 0.15*total_density
    density['co1'][:,:] = 0.062*total_density
    density['co2'][:,:] = 0.248*total_density
    density['co3'][:,:] = 0.248*total_density
    density['co4'][:,:] = 0.062*total_density

    # Growth magnitude
    Growth = np.ones((mesh.nelem,NoGauss),dtype=np.float64)
    
    # Deposition Stretch
    DepositionStretch = {}
    DepositionStretch['ela'] = np.zeros((mesh.nelem,NoGauss,3,3),dtype=np.float64)
    DepositionStretch['ela'][:,:,2,2] = 1.25
    DepositionStretch['ela'][:,:,1,1] = 1.34
    DepositionStretch['ela'][:,:,0,0] = 0.60
    DepositionStretch['smc'] = 1.1
    DepositionStretch['col'] = 1.062

    # 2D neo-Hookean parameter
    mu2D = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    mu2D[:,:] = 25.0e3

    # Directions for fibres and element orientations
    Directions(mesh,LagrangeGaussCoords,fibre_direction)
    
    # Store some total initial density and elastin
    den0_tot = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    den0_tot[:,:] = total_density
    den0_e = density['ela']

    # Define hyperelastic material for mesh
    material = ArterialWallMixture_1(ndim,
                 mu2D=mu2D,
                 mu3D=72.0e3,
                 c1m=15.2e3,
                 c2m=11.4,
                 c1c=1136.0e3,
                 c2c=11.2,
                 kappa=72.0e4,
                 anisotropic_orientations=fibre_direction,
                 deposition_stretch=DepositionStretch,
                 mix_density=density,
                 Growth=Growth,
                 Remodeling=Remodeling,
                 TotalDeformation=TotalDeformation
                 )

#===============  BOUNDARY CONDITIONS  ====================
    # Dirichlet Boundary Conditions
    def Dirichlet_Function(mesh, DirichletBoundary):
        #num_of_points=mesh.point.shape[0], ndim=3 ==> boundary_data(npoin,ndim)
        boundary_data = np.zeros((mesh.nnode, 3))+np.NAN
        # boundary conditions base on BoundarySurface boolean array
        boundary_data[DirichletBoundary['Bottom'],1] = 0.
        boundary_data[DirichletBoundary['Top'],1] = 0.
        boundary_data[DirichletBoundary['SymmetryZ'],2] = 0.
        boundary_data[DirichletBoundary['SymmetryX'],0] = 0.

        return boundary_data

    # Neumann Boundary Conditions
    def Neumann_Function(mesh, NeumannBoundary):
        boundary_flags = np.zeros(mesh.faces.shape[0],dtype=np.uint8)
        boundary_data = np.zeros((mesh.faces.shape[0], 3))
        # Force magnitud
        mag = 13.332e-3

        face = 0
        for idf in range(mesh.faces.shape[0]):
            if NeumannBoundary['InnerLogic'][idf]:
                coord = mesh.points[NeumannBoundary['InnerFaces'][face,:],:]
                vector = np.sum(coord,axis=0)/mesh.faces.shape[1]
                distance = np.sqrt(np.square(vector[0]) + np.square(vector[2]))
                normal = vector/distance
                mag0 = normal[0]*mag
                mag2 = normal[2]*mag
                boundary_data[idf,0] = mag0
                boundary_data[idf,2] = mag2
                face = face + 1

        boundary_flags[NeumannBoundary['InnerLogic']] = True

        return boundary_flags, boundary_data

#===============  SOLVER DEFINITION  ======================
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
    boundary_condition.SetNeumannCriteria(Neumann_Function, mesh, NeumannBoundary)

    fem_solver = FEMSolver( number_of_load_increments=1,
                            analysis_nature="nonlinear",
                            analysis_type="static",
                            break_at_stagnation=False,
                            maximum_iteration_for_newton_raphson=150,
                            optimise=False,
                            print_incremental_log=True,
                            memory_store_frequency=1)

#=================  HOMEOSTATIC SOLUTION  =======================
    # Homeostatic step solution
    print('=====================================')
    print('==  COMPUTE HOMEOSTATIC STATE  ==')
    print('=====================================')
    # Compute Radial Deposition Stretch at Gauss points
    ElastinDepositionStretch(material=material, R_gauss=R_gauss, nelem=mesh.nelem, NoGauss=NoGauss)
    # Compute neo-Hookean 2D parameter at Gauss points
    CircumferentialStress(material=material, nelem=mesh.nelem, NoGauss=NoGauss)
    # Call the solver for Homeostatic computation
    fem_solver.newton_raphson_tolerance=1.0e-5
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, 
            material=material, boundary_condition=boundary_condition)
    print('=> Homeostatic step finished')
    # Write Homeostatic state to paraview
    solution.WriteVTK('GandR_0',quantity=0)
    # Store Homeostatic Stress
    HomeostaticStress(material=material, Stress_H=Stress_H, nelem=mesh.nelem , NoGauss=NoGauss)

#=================  REMODELING SOLUTION  =======================
    print('=====================================')
    print('==  COMPUTE GROWTH AND REMODELING  ==')
    print('=====================================')
    Y_square = np.square(LagrangeGaussCoords[:,:,1])
    # Growth and Remodeling steps [30 days]
    total_time = 5500
    time = 0.0
    step = 0
    fem_solver.newton_raphson_tolerance=1.0e-5
    while time<total_time:
        # choose a Dt
        Delta_t = 10.0
        # Update Euler coordinates
        TotalDisplacements = solution.sol[:,:,0]
        euler_x = mesh.points + TotalDisplacements
        # initialize time total density
    #**** compute of G&R at t_n **** F_{gr}(t_{n+1})
        GrowthRemodeling(time=time, material=material, Stress_H=Stress_H,
                Y_square=Y_square, den0_tot=den0_tot, den0_e=den0_e, 
                nelem=mesh.nelem, NoGauss=NoGauss)

    #**** compute properties at t_{n+1} **** F_{gr}(t_{n+1})
        step += 1
        time += Delta_t
        print('==== STEP: '+str(step)+' |8===D| TIME: '+str(time)+' days ====')
        print('*** Compute Solution')

    #**** compute mechanical equilibrium at t_{n+1} **** F(t_{n+1})
        solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material, 
                boundary_condition=boundary_condition,Eulerx=euler_x)
        # Check Residual
        if np.isnan(fem_solver.norm_residual) or fem_solver.norm_residual>1e06:
            fem_solver.newton_raphson_failed_to_converge = True
            print('MODEL DID NOT CONVERGE!')
            break

    #**** OUTPUT ****
        solution.WriteVTK('GandR_'+str(step),quantity=0)


if __name__ == "__main__":
    homogenized_CMT()

