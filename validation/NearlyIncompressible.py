# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Kuru
#sys.path.append(os.path.join(os.path.expanduser("~/kuru")))
sys.path.append(os.path.expanduser("~/kuru"))
#import Florence
from Kuru import *

#============================================================
#============= ANISOTROPIC FIBRE DIRECTIONS  ================
#============================================================
def Directions(mesh):
    """
        Routine dedicated to compute the fibre direction of components in integration point for 
        the Material in Florence and for the auxiliar routines in this script. First three directions 
        are taken into the code for Rotation matrix, so always it should be present in this order,
        Normal, Tangential, Axial.
    """
    ndim = mesh.InferSpatialDimension()
    nfibre = 6
    direction = np.zeros((mesh.nelem,nfibre,ndim),dtype=np.float64)
    # Geometric definitions per element
    divider = mesh.elements.shape[1]
    directrix = [0.,1.,0.]
    # Loop throught the element in the mesh
    for elem in range(mesh.nelem):
        # Geometric definitions per element
        center = np.sum(mesh.points[mesh.elements[elem,:],:],axis=0)/divider
        tangential = np.cross(directrix,center)
        tangential = tangential/np.linalg.norm(tangential)
        normal = np.cross(tangential,directrix)
        direction[elem][0][:] = normal/np.linalg.norm(normal)
        # Define the anisotropic orientations
        direction[elem][1][:]=tangential
        direction[elem][2][:]=directrix
        direction[elem][3][:]=np.multiply(directrix,np.cos(np.pi/4.)) + np.multiply(tangential,np.sin(np.pi/4.))
        direction[elem][4][:]=np.multiply(directrix,np.cos(np.pi/4.)) - np.multiply(tangential,np.sin(np.pi/4.))
        direction[elem][5][:]=tangential

    return direction

#============================================================
#===============  HOMOGENIZED CMT  ==========================
#============================================================
def NearlyIncompressible(p=1,optimise=True):

    ProblemPath = os.path.dirname(os.getcwd())
    mesh_file = ProblemPath + '/Quarter_Ring.msh'

    #===============  MESH PROCESING  ==========================
    # Build mesh with Florence tools from GMSH mesh
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex",read_surface_info=True)
    ndim = mesh.InferSpatialDimension()
    mesh.GetHighOrderMesh(p=p)
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
    PressureBoundary = {}
    PressureBoundary['InnerLogic'] = InnerSurface
    PressureBoundary['InnerFaces'] = InnerFaces

    #===============  MATERIAL DEFINITION  ====================
    field_variables = np.zeros((mesh.nnode,23),dtype=np.float64)
    # Deposition Stretches
    field_variables[:,0] = 1.0/(1.34*1.25)
    field_variables[:,4] = 1.34
    field_variables[:,8] = 1.25
    field_variables[:,9] = 1.1
    field_variables[:,10] = 1.062
    # Total initial density
    field_variables[:,11] = 241.5
    field_variables[:,12] = 157.5
    field_variables[:,13] = 65.1
    field_variables[:,14] = 260.4
    field_variables[:,15] = 260.4
    field_variables[:,16] = 65.1
    field_variables[:,17:23] = 1.0

    # fibre directions [thick,smc,co1,co2,co3,co4]
    fibre_direction = Directions(mesh)

    # Define hyperelastic material for mesh
    material = ArterialWallMixture(ndim,
            mu=72.0,
            kappa=72.0*33.0,
            k1m=7.6,
            k2m=11.4,
            k1c=568.0,
            k2c=11.2,
            anisotropic_orientations=fibre_direction,
            field_variables=field_variables)

    # kappa/mu=20  => nu=0.475 (Poisson's ratio)
    # kappa/mu=33  => nu=0.485 (Poisson's ratio)
    # kappa/mu=100 => nu=0.495 (Poisson's ratio)

    #==================  FORMULATION  =========================
    formulation = DisplacementFormulation(mesh)

    #===============  BOUNDARY CONDITIONS  ====================
    # Dirichlet Boundary Conditions
    def Dirichlet_Function(mesh, DirichletBoundary):
        boundary_data = np.zeros((mesh.nnode, 3))+np.NAN
        # boundary conditions base on BoundarySurface boolean array
        boundary_data[DirichletBoundary['Bottom'],1] = 0.
        boundary_data[DirichletBoundary['Top'],1] = 0.
        boundary_data[DirichletBoundary['SymmetryZ'],2] = 0.
        boundary_data[DirichletBoundary['SymmetryX'],0] = 0.

        return boundary_data

    # Pressure Boundary Conditions
    def Pressure_Function(mesh, PressureBoundary):
        boundary_flags = np.zeros(mesh.faces.shape[0],dtype=np.uint8)
        boundary_data = np.zeros((mesh.faces.shape[0]))
        # Force magnitud
        mag = 13.3322e3

        for idf in range(mesh.faces.shape[0]):
            if PressureBoundary['InnerLogic'][idf]:
                boundary_data[idf] = mag

        boundary_flags[PressureBoundary['InnerLogic']] = True

        return boundary_flags, boundary_data

    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
    boundary_condition.SetPressureCriteria(Pressure_Function, mesh, PressureBoundary)

    #===============  SOLVER DEFINITION  ======================
    fem_solver = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=optimise,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

    #=================  COMPUTE SOLUTION  =======================
    # Call the solver
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
        material=material, boundary_condition=boundary_condition)
    # Check the error displacement
    print(solution.sol[0,:,-1])
    print(solution.sol[1,:,-1])

if __name__ == "__main__":
    NearlyIncompressible(optimise=True)
