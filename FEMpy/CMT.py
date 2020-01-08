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
PressureBoundary = {}
PressureBoundary['InnerLogic'] = InnerSurface
PressureBoundary['InnerFaces'] = InnerFaces

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
Deposition['Matrix'] = np.zeros((3,3),dtype=np.float64)
Deposition['Fibre'] = np.ones((5),dtype=np.float64)
Deposition['Matrix'][2,2] = 1.25
Deposition['Matrix'][1,1] = 1.34
Deposition['Matrix'][0,0] = 1./(1.25*1.34)
Deposition['Fibre'][0] = 1.1
Deposition['Fibre'][1:5] = 1.062

# fibre directions [thick,sms,co1,co2,co3,co4]
fibre_direction = Directions(mesh)

# Define hyperelastic material for mesh
material = ArterialWallMixture_(ndim,
            mu3D=72.0,
            c1m=15.2,
            c2m=11.4,
            c1c=1136.0,
            c2c=11.2,
            kappa=72.0e6,
            anisotropic_orientations=fibre_direction,
            Deposition=Deposition,
            GrowthRemodeling=GrowthRemodeling)

# Define hyperelastic material for mesh
#material = NearlyIncompressibleNeoHookean(ndim,
#            mu=72.0*1.e3,
#            kappa=72.0e5*1.e3)

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
                       optimise=False,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

#=================  HOMEOSTATIC SOLUTION  =======================
# Homeostatic step solution
print('=====================================')
print('==  COMPUTE HOMEOSTATIC STATE  ==')
print('=====================================')
# Call the solver for Homeostatic computation
for Iter in range(1):
    print('Iterarion ',Iter)
    # Call FEM solver for the current state
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
        material=material, boundary_condition=boundary_condition)
    # Check the error displacement
    dmesh = Mesh()
    dmesh.points = solution.sol[:,:,-1]
    dmesh_bounds = dmesh.Bounds
    distortion = np.sqrt(dmesh_bounds[0,0]**2+dmesh_bounds[0,1]**2+dmesh_bounds[0,2]**2)/0.010
    print('Distortion: '+str(distortion))
    #if distortion<0.05: break
    # GET DEFOMATION GRADIENT AT NODES TO COMPUTE A NEW ELASTIN DEPOSITION STRETCH
    #solution.StressRecovery()
    #DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
    # Compute Deposition Stretch at Gauss points
    #Deposition['Matrix'] = PrestrainGradient(Deposition['Matrix'], DeformationGradient, mesh)
    #print(DepositionStretch['ela'])
    #material.Deposition['Matrix']=Deposition['Matrix']
# Write Homeostatic state to paraview
solution.WriteVTK('GandR_0',quantity=0)
print('... Homeostatic step finished')
# HOMEOSTATIC POSTPROCESS
#solution.StressRecovery()
#DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
#Stress_H = solution.recovered_fields['FibreStress'][-1,:,:]
#FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
#Softness = solution.recovered_fields['Softness'][-1,:,:]
# Update mesh coordinates
#TotalDisplacements = solution.sol[:,:,-1]
#euler_x = mesh.points + TotalDisplacements
"""
#=================  REMODELING SOLUTION  =======================
print('=====================================')
print('==  COMPUTE GROWTH AND REMODELING  ==')
print('=====================================')
#print(mesh.points[mesh.elements,:])
# Growth and Remodeling steps [10 days]
total_time = 5500
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
    GrowthRemodeling = GetGrowthRemodeling(time,Delta_t,mesh,GrowthRemodeling,Stress_H,FibreStress,Softness)
    material.GrowthRemodeling = GrowthRemodeling
    #print(GrowthRemodeling)
    #**** compute mechanical equilibrium at t_{n+1} **** F(t_{n+1})
    solution = fem_solver.Solve(formulation=formulation, mesh=mesh, material=material,
           boundary_condition=boundary_condition, Eulerx=euler_x)
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
    euler_x = mesh.points + TotalDisplacements
"""
