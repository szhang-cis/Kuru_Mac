# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
# Build a path for python to Florence
sys.path.append(os.path.join(os.path.expanduser("~"),"florence"))
#import Florence
from Florence import *

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

    Y_Gauss = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    for elem in range(mesh.nelem):
        LagrangeCoords = mesh.points[mesh.elements[elem,:],:]
        LagrangeGaussCoords = np.einsum('ij,ik->jk',Bases,LagrangeCoords)
        Y_Gauss[elem,:] = LagrangeGaussCoords[:,1]

#===============  MATERIAL DEFINITION  ====================
    # Geometric definitions per element
    center = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    tangential = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    norma1 = np.zeros((mesh.nelem),dtype=np.float64)
    norma2 = np.zeros((mesh.nelem),dtype=np.float64)
    directrix = [0.,1.,0.]
    divider = mesh.elements.shape[1]
    # Dictionary for fibres orientations
    fibre_direction = {}
    for key in ['co1','co2','co3','co4','smc','thick']:
        fibre_direction[key] = np.zeros((mesh.nelem,ndim),dtype=np.float64)

    # Dictionary for Density mixture and Inelastic gradient
    total_density = 1050.0e-9
    density = {}
    density_rate = {}
    for key in ['ela','co1','co2','co3','co4','smc']:
        density[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
        density_rate[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)

    # Dictionary for Inelastic gradient
    F_gr = {}
    for key in ['ela','co1','co2','co3','co4','smc']:
        F_gr[key] = np.zeros((mesh.nelem,NoGauss,3,3),dtype=np.float64)
        F_gr[key][:,:,0,0] = 1.0
        F_gr[key][:,:,1,1] = 1.0
        F_gr[key][:,:,2,2] = 1.0

    # Dictionary for Directional Stresses and Remodeling 
    Stress_H = {}
    Remodeling = {}
    Remodeling_rate = {}
    for key in ['co1','co2','co3','co4','smc']:
        Stress_H[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
        Remodeling[key] = np.ones((mesh.nelem,NoGauss),dtype=np.float64)
        Remodeling_rate[key] = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)

    # Densities per constituent
    density['ela'][:,:] = 0.23*total_density
    density['smc'][:,:] = 0.15*total_density
    density['co1'][:,:] = 0.062*total_density
    density['co2'][:,:] = 0.248*total_density
    density['co3'][:,:] = 0.248*total_density
    density['co4'][:,:] = 0.062*total_density
    Growth = np.ones((mesh.nelem,NoGauss),dtype=np.float64)
    # Loop throught the element in the mesh
    for i in range(mesh.nelem):
        # Geometric definitions per element
        center[i,:] = np.sum(mesh.points[mesh.elements[i,:],:],axis=0)/divider
        tangential[i,:] = np.cross(directrix,center[i,:])
        norma1 = np.linalg.norm(tangential[i,:])
        tangential[i,:] = tangential[i,:]/norma1
        fibre_direction['thick'][i,:] = np.cross(tangential[i,:],directrix)
        norma2 = np.linalg.norm(fibre_direction['thick'][i,:])
        fibre_direction['thick'][i,:] = fibre_direction['thick'][i,:]/norma2
        # Define the anisotropic orientations
        fibre_direction['smc'][i,:] = np.multiply(tangential[i,:],np.sin(np.pi/2))
        fibre_direction['co1'][i,:] = np.multiply(tangential[i,:],np.sin(np.pi/2))
        fibre_direction['co2'][i,:] = np.multiply(directrix,np.cos(np.pi/4)) + \
                np.multiply(tangential[i,:],np.sin(np.pi/4))
        fibre_direction['co3'][i,:] = np.multiply(directrix,np.cos(np.pi/4)) - \
                np.multiply(tangential[i,:],np.sin(np.pi/4))
        fibre_direction['co4'][i,:] = np.multiply(directrix,np.cos(0.))

    # Deposition Stretch
    DepStr = {}
    DepStr['ela'] = np.eye(ndim,ndim,dtype=np.float64)
    DepStr['ela'][2,2] = 1.25
    DepStr['ela'][1,1] = 1.34
    DepStr['ela'][0,0] = 1./(np.multiply(DepStr['ela'][1,1],DepStr['ela'][0,0]))
    DepStr['smc'] = 1.1
    DepStr['col'] = 1.062

    # Define hyperelastic material for mesh
    material = ArterialWallMixture(ndim,
                 mu2D=0.0e3,
                 mu3D=72.0e3,
                 c1m=15.2e3,
                 c2m=11.4,
                 c1c=1136.0e3,
                 c2c=11.2,
                 kappa=72.0e6,
                 anisotropic_orientations=fibre_direction,
                 mix_density=density,
                 density_rate=density_rate,
                 deposition_stretch=DepStr,
                 Stress_H=Stress_H,
                 Growth=Growth,
                 Remodeling=Remodeling,
                 Remodeling_rate=Remodeling_rate
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
    material.label = 'Homeostatic'
    # Call the solver for Homeostatic computation
    fem_solver.newton_raphson_tolerance=1.0e-5
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, 
            material=material, boundary_condition=boundary_condition)
    # Write Homeostatic state to paraview
    solution.WriteVTK('GandR_0',quantity=0)

#=================  REMODELING SOLUTION  =======================
    print('=====================================')
    print('==  COMPUTE GROWTH AND REMODELING  ==')
    print('=====================================')
    material.label = 'GandR'
    # Store some Homeostatic values for G&R
    den0_tot = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
    den0_tot[:,:] = total_density
    den0_e = density['ela']
    square = np.square(Y_Gauss)
    # Growth and Remodeling steps [30 days]
    total_time = 5500
    time = 0.0
    step = 0
    fem_solver.newton_raphson_tolerance=1.0e-5,
    while time<total_time:
        # choose a Dt
        #if time<300.0:
        #    Delta_t = 20.0
        #elif time>2000.0:
        #    Delta_t = 60.0
        #else:
        Delta_t = 1.0

        # Update coordinates
        TotalDisplacements = solution.sol[:,:,0]
        euler_x = mesh.points + TotalDisplacements
        # initialize time total density
        den_tot = np.zeros((mesh.nelem,NoGauss),dtype=np.float64)
        # Elastin degradation [months]
        L_dam = 10.0
        t_dam = 40.0
        D_max = 0.5
        T_ela = 101.0*365.25
    #**** COMPUTE K1 euler-METHOD **** f(t_n,w_n)
        # Update elastin density
        exponential = np.exp(-0.5*square/L_dam**2 - time/t_dam)
        density_rate['ela']= -density['ela']/T_ela - (D_max/t_dam)*np.multiply(exponential,den0_e)

        step += 1
        time += Delta_t
        print('==== STEP: '+str(step)+' |8===D| TIME: '+str(time)+' days ====')
    #**** COMPUTE THE SOLUTION **** w_n+1=w_n+h*k1
        print('*** Compute Solution')
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

        # Solution algorithm (get f in [n+1])
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

