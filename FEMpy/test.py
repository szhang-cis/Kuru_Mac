# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
from numpy import einsum
# Build a path for python to Florence
sys.path.append(os.path.join(os.path.expanduser("~/femme"),"FEMpy"))
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
#============================================================
#===============  HOMOGENIZED CMT  ==========================
#============================================================
ProblemPath = os.getcwd()
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
NeumannBoundary = {}
NeumannBoundary['InnerLogic'] = InnerSurface
NeumannBoundary['InnerFaces'] = InnerFaces

#===============  MATERIAL DEFINITION  ====================
# fibre directions [thick,sms,co1,co2,co3,co4]
fibre_direction = Directions(mesh)

# Define hyperelastic material for mesh
#material = ArterialWallMixture(ndim,
#            mu3D=72.0,
#            c1m=15.2,
#            c2m=11.4,
#            c1c=1136.0,
#            c2c=11.2,
#            kappa=72.0e3,
#            anisotropic_orientations=fibre_direction)

# Define hyperelastic material for mesh
material = NearlyIncompressibleNeoHookean(ndim,
            mu=72.0*5000.0,
            kappa=72.0e4*5000.0)

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
                distances[i] = np.linalg.norm(vertex)
            
            semi_diagonal = max(distances)
            area = 2.*semi_diagonal**2
            boundary_data[idf,0] = normal[0]*mag*area
            boundary_data[idf,2] = normal[2]*mag*area
            face += 1

    boundary_flags[NeumannBoundary['InnerLogic']] = True

    return boundary_flags, boundary_data

boundary_condition = BoundaryCondition()
boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, DirichletBoundary)
boundary_condition.SetNeumannCriteria(Neumann_Function, mesh, NeumannBoundary)

#===============  SOLVER DEFINITION  ======================
fem_solver = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=False,
                       print_incremental_log=True,
                       number_of_load_increments=2)

#===============  COMPUTE SOLUTION  ======================
# Call FEM solver for the current state
solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
    material=material, boundary_condition=boundary_condition)
# Write Homeostatic state to paraview
solution.WriteVTK('HiperElastic',quantity=0)
