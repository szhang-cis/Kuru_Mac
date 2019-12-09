# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Florence
sys.path.append(os.path.join(os.path.expanduser("~/femme"),"florence"))
#import Florence
from Florence import *

def Directions(mesh):
    """
        Routine dedicated to compute the fibre direction of components in integration point for 
        the Material in Florence and for the auxiliar routines in this script.
    """
    ndim = mesh.InferSpatialDimension()
    direction = np.zeros((6,mesh.nelem,ndim),dtype=np.float64)
    # Geometric definitions per element
    center = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    tangential = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    divider = mesh.elements.shape[1]
    directrix = [0.,1.,0.]
    # Loop throught the element in the mesh
    for elem in range(mesh.nelem):
        # Geometric definitions per element
        center[elem,:] = np.sum(mesh.points[mesh.elements[elem,:],:],axis=0)/divider
        tangential[elem,:] = np.cross(directrix,center[elem,:])
        tangential[elem,:] = tangential[elem,:]/np.linalg.norm(tangential[elem,:])
        direction[0][elem,:] = np.cross(tangential[elem,:],directrix)
        direction[0][elem,:] = direction[0][elem,:]/np.linalg.norm(direction[0][elem,:])
        # Define the anisotropic orientations
        direction[1][elem,:]=tangential[elem,:]
        direction[2][elem,:]=tangential[elem,:]
        direction[3][elem,:]=np.multiply(directrix,np.cos(np.pi/4)) + \
            np.multiply(tangential[elem,:],np.sin(np.pi/4))
        direction[4][elem,:]=np.multiply(directrix,np.cos(np.pi/4)) - \
            np.multiply(tangential[elem,:],np.sin(np.pi/4))
        direction[5][elem,:]=directrix

    return direction

def PrestrainGradient(ElastinDepositionStretch, DeformationGradient, mesh):
    """ 
        This routine aim to get the radial deposition stretch from the boundary conditions 
        in radial direction, as stress=-p in inner wall and stress=0 in outer wall. The 
        relation for the stresses along the thickness come from thin-walled tube. So radial
        deposition stretch is depending of the position (radius).
    """

    print('==> Looking for Elastin Deposition Stretch')
    F_average = np.zeros((3,3),dtype=np.float64)
    for node in range(mesh.nnode):
        # Directional vector for element
        directrix = [0.,1.,0.]
        Axial = directrix
        Tangential = np.cross(directrix,mesh.points[node,:])
        Tangential = Tangential/np.linalg.norm(Tangential)
        Normal = np.cross(Tangential,directrix)
        Normal = Normal/np.linalg.norm(Normal)
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = Axial[i]
        F = DeformationGradient[node]
        F = np.dot(Rotation,np.dot(F,Rotation.T))
        F_average += F
    F_average = np.divide(F_average,mesh.nnode)
    for node in range(mesh.nnode):
        ElastinDepositionStretch[node] = np.dot(F_average,ElastinDepositionStretch[node])
    print('==> Elastin Deposition Stretch calculated')
    return ElastinDepositionStretch

def GetGrowthRemodelingRates(time,mesh,Stress_H,FibreStress,Softness,GrowthRemodeling):
    """
        This are the rates of Growth and Remodeling.
        Forward Euler.
    """
    GrowthRemodeling_rate = np.zeros((mesh.nnode,11),dtype=np.float64)
    # Elastin degradation [days]
    L_dam = 0.010
    t_dam = 40.0
    D_max = 0.5
    T_ela = 101.0*365.25
    den0_e = 1050.0*0.23
    # Fribres turnover and remodeling
    turnover = 101.0
    gain = 0.05/turnover
    for node in range(mesh.nnode):
        AxialCoord = mesh.points[node,1]
        # Update elastin density
        GrowthRemodeling_rate[node][0] = -GrowthRemodeling[node][0]/T_ela - \
            (D_max/t_dam)*np.multiply(np.exp(-0.5*(AxialCoord/L_dam)**2 - time/t_dam),den0_e)
        #SMC AND COLLAGEN FIBRES
        for key in [1,2,3,4,5]:
            # fibres density rates
            DeltaStress = FibreStress[node][key-1] - Stress_H[node][key-1]
            GrowthRemodeling_rate[node][key] = gain*np.multiply(GrowthRemodeling[node][key],
                np.divide(DeltaStress,Stress_H[node][key-1]))
            # fibres remodeling rates
            lambda_r_dot = np.divide(GrowthRemodeling_rate[node][key],GrowthRemodeling[node][key]) + 1./turnover
            GrowthRemodeling_rate[node][key+5] = np.multiply(np.multiply(lambda_r_dot,DeltaStress),
                Softness[node][key-1])

    return GrowthRemodeling_rate

def GetGrowthRemodeling(Delta_t,mesh,GrowthRemodeling,GrowthRemodeling_rate):
    """
        This are the rates of Growth and Remodeling.
        Forward Euler.
    """
    den0_tot = 1050.0
    for node in range(mesh.nnode):
        den_tot = 0.0
        for key in range(11):
            GrowthRemodeling[node][key] += Delta_t*GrowthRemodeling_rate[node][key]
            if key<6: den_tot += GrowthRemodeling[node][key]

        GrowthRemodeling[node][11] = den_tot/den0_tot

    return GrowthRemodeling

def UpdateDeposition(DepositionMatrix, DepositionFibre, DeformationGradient, mesh):
    """ 
        This routine aim to get the radial deposition stretch from the boundary conditions 
        in radial direction, as stress=-p in inner wall and stress=0 in outer wall. The 
        relation for the stresses along the thickness come from thin-walled tube. So radial
        deposition stretch is depending of the position (radius).
    """

    print('==> Updating Deposition Stretches')
    directrix = [0.,1.,0.]
    for node in range(mesh.nnode):
        # Directional vector for element
        #Axial = directrix
        Tangential = np.cross(directrix,mesh.points[node,:])
        Tangential = Tangential/np.linalg.norm(Tangential)
        Normal = np.cross(Tangential,directrix)
        Normal = Normal/np.linalg.norm(Normal)
        Rotation = np.eye(3,3,dtype=np.float64)
        for i in range(3):
            Rotation[0,i] = Normal[i]
            Rotation[1,i] = Tangential[i]
            Rotation[2,i] = directrix[i]
        F = DeformationGradient[node]
        # Elastin
        F = np.dot(Rotation,np.dot(F,Rotation.T))
        DepositionMatrix[node] = np.dot(F,DepositionMatrix[node])
        # Smooth Muscle Cells
        FN = np.dot(DeformationGradient[node],Tangential)
        innerFN = einsum('i,i',FN,FN)
        DepositionFibre[node][0] = DepositionFibre[node][0]*np.sqrt(innerFN)
        # Collagen 1
        DepositionFibre[node][1] = DepositionFibre[node][1]*np.sqrt(innerFN)
        # Collagen 2
        N = np.multiply(directrix,np.cos(np.pi/4.)) + np.multiply(Tangential,np.sin(np.pi/4.))
        FN = np.dot(DeformationGradient[node],N)
        innerFN = einsum('i,i',FN,FN)
        DepositionFibre[node][2] = DepositionFibre[node][2]*np.sqrt(innerFN)
        # Collagen 3
        N = np.multiply(directrix,np.cos(np.pi/4.)) - np.multiply(Tangential,np.sin(np.pi/4.))
        FN = np.dot(DeformationGradient[node],N)
        innerFN = einsum('i,i',FN,FN)
        DepositionFibre[node][3] = DepositionFibre[node][3]*np.sqrt(innerFN)
        # Collagen 4
        FN = np.dot(DeformationGradient[node],directrix)
        innerFN = einsum('i,i',FN,FN)
        DepositionFibre[node][4] = DepositionFibre[node][4]*np.sqrt(innerFN)
    
    print('==> Deposition Stretches updated')
    return DepositionMatrix,DepositionFibre

def homogenized_CMT():
    """A hyperelastic implicit static example using ArterialMixture model
        of a cylinder under pression with hexahedral elements
    """
    ProblemPath = PWD(__file__)
    mesh_file = ProblemPath + '/Quarter_Ring.msh'

#===============  MESH PROCESING  ==========================
    # Build mesh with Florence tools from GMSH mesh
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex",read_surface_info=True)
    #mesh.GetHighOrderMesh(p=3)
    ndim = mesh.InferSpatialDimension()
    #Boolean arrays for boundary condition in Dirichlet
    BottomSurface = np.zeros(mesh.nnode,dtype=bool)
    TopSurface = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_Z = np.zeros(mesh.nnode,dtype=bool)
    Symmetry_X = np.zeros(mesh.nnode,dtype=bool)
    #Boolean array for boundary condition in Neumann mesh.faces[id,node]
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

#===============  MATERIAL DEFINITION  ====================
    # Total initial density
    total_density = 1050.0
    # Growth and remodeling variables (densities, remodeling in fibres and growth)
    GrowthRemodeling = np.zeros((mesh.nnode,12),dtype=np.float64)
    # denstities[ela,smc,co1,co2,co3,co4], remodeling[smc,co1,co2,co3,co4], growth
    GrowthRemodeling[:,0] = 0.23*total_density
    GrowthRemodeling[:,1] = 0.15*total_density
    GrowthRemodeling[:,2] = 0.062*total_density
    GrowthRemodeling[:,3] = 0.248*total_density
    GrowthRemodeling[:,4] = 0.248*total_density
    GrowthRemodeling[:,5] = 0.062*total_density
    GrowthRemodeling[:,6:12] = 1.0
    
    # Deposition Stretch
    Deposition = {}
    Deposition['Matrix'] = np.zeros((mesh.nnode,3,3),dtype=np.float64)
    Deposition['Fibre'] = np.ones((mesh.nnode,5),dtype=np.float64)
    Deposition['Matrix'][:,2,2] = 1.25
    Deposition['Matrix'][:,1,1] = 4.17
    Deposition['Matrix'][:,0,0] = 1.0/(1.25*4.17)
    Deposition['Fibre'][:,0] = 1.1
    Deposition['Fibre'][:,1:5] = 1.062

    # fibre directions [thick,sms,co1,co2,co3,co4]
    fibre_direction = Directions(mesh)

    # Define hyperelastic material for mesh
    material = ArterialWallMixture_ex(ndim,
                mu3D=72.0,
                c1m=15.2,
                c2m=11.4,
                c1c=1136.0,
                c2c=11.2,
                kappa=72.0e4,
                anisotropic_orientations=fibre_direction,
                Deposition=Deposition,
                GrowthRemodeling=GrowthRemodeling)

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
        mag = 13.3322e3

        face = 0
        for idf in range(mesh.faces.shape[0]):
            if NeumannBoundary['InnerLogic'][idf]:
                vertex = np.zeros(3)
                distances = np.zeros(mesh.faces.shape[1])
                coord = mesh.points[NeumannBoundary['InnerFaces'][face,:],:]
                # center point in the surface
                midpoint = np.sum(coord,axis=0)/mesh.faces.shape[1]
                midpoint_mag = np.sqrt(np.square(midpoint[0]) + np.square(midpoint[2]))
                # normal to surface
                normal = midpoint/midpoint_mag
                # area of the face
                for i in range(NeumannBoundary['InnerFaces'].shape[1]):
                    vertex = coord[i,:] - midpoint
                    distances[i] = la.norm(vertex)
            
                semi_diagonal = max(distances)
                area = 2.*semi_diagonal**2
                mag0 = normal[0]*mag*area
                mag2 = normal[2]*mag*area
                boundary_data[idf,0] = mag0
                boundary_data[idf,2] = mag2
                face += 1

        boundary_flags[NeumannBoundary['InnerLogic']] = True

        return boundary_flags, boundary_data

#===============  SOLVER DEFINITION  ======================
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
    boundary_condition.SetNeumannCriteria(Neumann_Function, mesh, NeumannBoundary)

    fem_solver0 = FEMSolver( analysis_nature="nonlinear",
                            analysis_type="static",
                            break_at_stagnation=False,
                            maximum_iteration_for_newton_raphson=50,
                            optimise=False,
                            print_incremental_log=True,
                            newton_raphson_tolerance=1.0e-5,
                            number_of_load_increments=1)

    fem_solver1 = FEMSolver( analysis_nature="nonlinear",
                            analysis_type="static",
                            break_at_stagnation=False,
                            maximum_iteration_for_newton_raphson=50,
                            optimise=False,
                            print_incremental_log=True,
                            newton_raphson_tolerance=1.0e-5,
                            number_of_load_increments=2)

#=================  HOMEOSTATIC SOLUTION  =======================
    # Homeostatic step solution
    print('=====================================')
    print('==  COMPUTE HOMEOSTATIC STATE  ==')
    print('=====================================')
    # Call the solver for Homeostatic computation
    for Iter in range(10):
        print('Iterarion ',Iter)
        # Call FEM solver for the current state
        solution = fem_solver0.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)
        # Check the error displacement
        dmesh = Mesh()
        dmesh.points = solution.sol[:,:,-1]
        dmesh_bounds = dmesh.Bounds
        distortion = np.sqrt(dmesh_bounds[0,0]**2+dmesh_bounds[0,1]**2+dmesh_bounds[0,2]**2)/0.010
        print('Distortion: '+str(distortion))
        if distortion<0.02: break
        # GET DEFOMATION GRADIENT AT NODES TO COMPUTE A NEW ELASTIN DEPOSITION STRETCH
        solution.StressRecovery()
        DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
        # Compute Deposition Stretch at Gauss points
        Deposition['Matrix'] = PrestrainGradient(Deposition['Matrix'], DeformationGradient, mesh)
        #print(DepositionStretch['ela'])
        material.Deposition['Matrix']=Deposition['Matrix']
    # Write Homeostatic state to paraview
    solution.WriteVTK('GandR_0',quantity=0)
    print('... Homeostatic step finished')
    # HOMEOSTATIC POSTPROCESS
    solution.StressRecovery()
    DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
    Stress_H = solution.recovered_fields['FibreStress'][-1,:,:]
    FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    Softness = solution.recovered_fields['Softness'][-1,:,:]
    # Update mesh coordinates
    #TotalDisplacements = solution.sol[:,:,-1]
    #mesh.points += TotalDisplacements
    #Deposition['Matrix'],Deposition['Fibre'] = UpdateDeposition(Deposition['Matrix'],Deposition['Fibre'],DeformationGradient,mesh)
    #print(Deposition['Matrix'])
    #print(Deposition['Fibre'])
    # Check the displacement
    #print('Maximum and Minimum coordinates at this point')
    #print(mesh.Bounds)

#=================  REMODELING SOLUTION  =======================
    print('=====================================')
    print('==  COMPUTE GROWTH AND REMODELING  ==')
    print('=====================================')
    #print(mesh.points[mesh.elements,:])
    # Growth and Remodeling steps [10 days]
    total_time = 1000
    time = 0.0
    Delta_t = 10.0
    step = 0
    while time<total_time:
        # prepare for next step
        time += Delta_t
        step += 1
        print('==== STEP: '+str(step)+' |8===D| TIME: '+str(time)+' days ====')
        print('*** Compute Solution')
    #**** compute of G&R at t_n **** F_{gr}(t_{n+1})
        GrowthRemodeling_rate=GetGrowthRemodelingRates(time,mesh,Stress_H,FibreStress,Softness,GrowthRemodeling)
        GrowthRemodeling = GetGrowthRemodeling(Delta_t,mesh,GrowthRemodeling,GrowthRemodeling_rate)
        material.GrowthRemodeling = GrowthRemodeling
        #print(GrowthRemodeling)
    #**** compute mechanical equilibrium at t_{n+1} **** F(t_{n+1})
        solution = fem_solver1.Solve(formulation=formulation, mesh=mesh, material=material, 
                boundary_condition=boundary_condition)
        # Check Residual
        if np.isnan(fem_solver1.norm_residual) or fem_solver1.norm_residual>1e06:
            exit('MODEL DID NOT CONVERGE!')

    #**** STEPS POSTPROCESS ****
        solution.WriteVTK('GandR_'+str(step),quantity=0)
        solution.StressRecovery()
        FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
        Softness = solution.recovered_fields['Softness'][-1,:,:]
        # Update mesh coordinates
        TotalDisplacements = solution.sol[:,:,-1]
        mesh.points += TotalDisplacements
        Deposition['Matrix'],Deposition['Fibre'] = UpdateDeposition(Deposition['Matrix'],Deposition['Fibre'],DeformationGradient,mesh)
        material.Deposition = Deposition
        # Check the displacement
        print('Maximum and Minimum coordinates at this point')
        print(mesh.Bounds)


if __name__ == "__main__":
    homogenized_CMT()

