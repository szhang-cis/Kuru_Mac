# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Kuru
sys.path.append(os.path.expanduser("~/kuru"))
from Kuru import *

def Directions(mesh):
    """
        Routine dedicated to compute the fibre direction of components in integration point for 
        the Material in Kuru and for the auxiliar routines in this script.
    """
    ndim = mesh.InferSpatialDimension()
    nfibre = 2
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
        # Define the anisotropic orientations
        direction[elem][0][:]=np.multiply(directrix,np.cos(np.pi/4)) + \
            np.multiply(tangential,np.sin(np.pi/4))
        direction[elem][1][:]=np.multiply(directrix,np.cos(np.pi/4)) - \
            np.multiply(tangential,np.sin(np.pi/4))

    return direction

#============================================================
#===============  FUNG QUADRATIC   ==========================
#============================================================
def FungQuad(optimise=False):
    ProblemPath = os.path.dirname(os.getcwd())
    mesh_file = ProblemPath + '/Quarter_Ring.msh'
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
    # Total initial density
    total_density = 525.0

    # fibre directions [thick,sms,co1,co2,co3,co4]
    fibre_direction = Directions(mesh)

    # Define hyperelastic material for mesh
    material = IncompressibleAnisotropicFungQuadratic(ndim,
            mu=72.0*total_density,
            kappa=72.0*total_density*33.0,
            k1=568.0*total_density,
            k2=11.2,
            anisotropic_orientations=fibre_direction)

    # kappa/mu=33 give nu=0.485 (Poisson's ratio)

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

    #=================  SOLUTION  =======================
    # Call the solver
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
        material=material, boundary_condition=boundary_condition)
    # Write to paraview
    #solution.WriteVTK('FungQ',quantity=0)
    #print(solution.sol[:,:,-1])

if __name__=="__main__":
    FungQuad(optimise=True)
