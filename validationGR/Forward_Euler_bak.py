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

def GetRadialStretch(mesh,material):
    print("Computing the Radial deposition stretch")
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    for node in range(mesh.nnode):
        kappa = material.kappa*material.growth_remodeling[node,0]
        mu3D = material.mu*material.growth_remodeling[node,0]
        stretch_r = material.deposition_stretch['Radial'][node]
        stretch_t = material.deposition_stretch['Tangential']
        stretch_z = material.deposition_stretch['Axial']
        radius = np.sqrt(mesh.points[node,0]**2+mesh.points[node,2]**2)
        JMAX=21
        for j in range(JMAX):
            if j==20: exit('Maximum number of iteration reach at Newton method at node: '+str(node))
            func = pressure*(1./2.-(radius-R)/H)+ \
                    kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.)+\
                    mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-\
                    1./3.*(stretch_r**2+stretch_t**2+stretch_z**2))
            df = kappa*(2.*stretch_r*stretch_t**2*stretch_z**2-stretch_t*stretch_z)+\
                    4./3.*mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3)*stretch_r-\
                    2./3.*mu3D*stretch_r**(-5./3.)*(stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-\
                    1./3.*(stretch_r**2+stretch_t**2+stretch_z**2))
            dx = func/df
            stretch_r -= dx
            if np.absolute(dx)<1.e-5:
                material.deposition_stretch['Radial'][node] = stretch_r
                break

    print("Radial deposition stretch is computed")

def GetGrowthRemodeling(time,Delta_t,mesh,GrowthRemodeling,Stress_H,FibreStress,Softness):
    """
        This are the rates of Growth and Remodeling.
        Forward Euler.
    """
    den0_tot = 1050.0
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
        elastin_rate = -GrowthRemodeling[node][0]/T_ela - \
            (D_max/t_dam)*np.multiply(np.exp(-0.5*(AxialCoord/L_dam)**2 - time/t_dam),den0_e)
        GrowthRemodeling[node][0] += Delta_t*elastin_rate
        den_tot = GrowthRemodeling[node][0]
        for key in [1,2,3,4,5]:
            # fibres density rates
            DeltaStress = FibreStress[node][key-1] - Stress_H[node][key-1]
            density_rate = gain*np.multiply(GrowthRemodeling[node][key],
                np.divide(DeltaStress,Stress_H[node][key-1]))
            GrowthRemodeling[node][key] += Delta_t*density_rate
            # fibres remodeling rates
            lambda_r_dot = np.divide(density_rate,GrowthRemodeling[node][key]) + 1./turnover
            remodeling_rate = np.multiply(np.multiply(lambda_r_dot,DeltaStress),Softness[node][key-1])
            GrowthRemodeling[node][key+5] += Delta_t*remodeling_rate
            den_tot += GrowthRemodeling[node][key]
        GrowthRemodeling[node][11] = den_tot/den0_tot

    return GrowthRemodeling

#============================================================
#===============  HOMOGENIZED CMT  ==========================
#============================================================
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
# Deposition Stretches
deposition_stretch = {}
deposition_stretch['Radial'] = np.zeros(mesh.nnode,dtype=np.float64)
deposition_stretch['Axial'] = 1.25
deposition_stretch['Tangential'] = 1.34
deposition_stretch['Radial'][:] = 1.0/(1.34*1.25)
deposition_stretch['Muscle'] = 1.1
deposition_stretch['Collagen'] = 1.062

# Total initial density
growth_remodeling = np.zeros((mesh.nnode,12),dtype=np.float64)
growth_remodeling[:,0] = 241.5
growth_remodeling[:,1] = 157.5
growth_remodeling[:,2] = 65.1
growth_remodeling[:,3] = 260.4
growth_remodeling[:,4] = 260.4
growth_remodeling[:,5] = 62.1
growth_remodeling[:,6:12] = 1.0

# fibre directions [thick,sms,co1,co2,co3,co4]
fibre_direction = Directions(mesh)

# Define hyperelastic material for mesh
material = ArterialWallMixtureGR(ndim,
            mu=72.0,
            kappa=72.0*33.0,
            k1m=7.6,
            k2m=11.4,
            k1c=568.0,
            k2c=11.2,
            anisotropic_orientations=fibre_direction,
            deposition_stretch=deposition_stretch,
            growth_remodeling=growth_remodeling)

# kappa/mu=20  => nu=0.475 (Poisson's ratio)
# kappa/mu=33  => nu=0.485 (Poisson's ratio)
# kappa/mu=100 => nu=0.495 (Poisson's ratio)
GetRadialStretch(mesh,material)
print(material.deposition_stretch['Radial'][:])

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
fem_solver_h = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=False,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

fem_solver_gr = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=False,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

#=================  HOMEOSTATIC SOLUTION  =======================
print('=====================================')
print('==  COMPUTE HOMEOSTATIC STATE  ==')
print('=====================================')
# Call the solver for Homeostatic computation
solution = fem_solver_h.Solve(formulation=formulation, mesh=mesh,
    material=material, boundary_condition=boundary_condition)
# Check the error displacement
dmesh = Mesh()
dmesh.points = solution.sol[:,:,-1]
dmesh_bounds = dmesh.Bounds
distortion = np.sqrt(dmesh_bounds[0,0]**2+dmesh_bounds[0,1]**2+dmesh_bounds[0,2]**2)/0.010
print('Distortion: '+str(distortion))
# Write Homeostatic state to paraview
solution.WriteVTK('ForwardEuler_0',quantity=0)
print('... Homeostatic step finished')
# HOMEOSTATIC POSTPROCESS
solution.StressRecovery()
#DeformationGradient = solution.recovered_fields['F'][-1,:,:,:]
Stress_H = solution.recovered_fields['FibreStress'][-1,:,:]
FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
Softness = solution.recovered_fields['Softness'][-1,:,:]
# Update mesh coordinates
TotalDisplacements = solution.sol[:,:,-1]
euler_x = mesh.points + TotalDisplacements
file_out = open("growth_remodeling.txt","w+")
file_out.write('%3d %f %f %f %f %f %f %f %f %f %f %f %f %f\n'%(0,euler_x[0,0],growth_remodeling[0,0],growth_remodeling[0,1],growth_remodeling[0,2],growth_remodeling[0,3],growth_remodeling[0,4],growth_remodeling[0,5],growth_remodeling[0,6],growth_remodeling[0,7],growth_remodeling[0,8],growth_remodeling[0,9],growth_remodeling[0,10],growth_remodeling[0,11]))
file_out.close()
#=================  REMODELING SOLUTION  =======================
print('=====================================')
print('==  COMPUTE GROWTH AND REMODELING  ==')
print('=====================================')
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
    growth_remodeling = GetGrowthRemodeling(time,Delta_t,mesh,growth_remodeling,Stress_H,FibreStress,Softness)
    material.growth_remodeling = growth_remodeling
    #**** compute mechanical equilibrium at t_{n+1} **** F(t_{n+1})
    solution = fem_solver_gr.Solve(formulation=formulation, mesh=mesh, material=material,
            boundary_condition=boundary_condition, Eulerx=euler_x)
    # Check Residual
    if np.isnan(fem_solver_gr.norm_residual) or fem_solver_gr.norm_residual>1e06:
        file_out.close()
        exit('MODEL DID NOT CONVERGE!')
    #**** STEPS POSTPROCESS ****
    solution.WriteVTK('ForwardEuler_'+str(step),quantity=0)
    solution.StressRecovery()
    FibreStress = solution.recovered_fields['FibreStress'][-1,:,:]
    Softness = solution.recovered_fields['Softness'][-1,:,:]
    # Update mesh coordinates
    TotalDisplacements = solution.sol[:,:,-1]
    euler_x = mesh.points + TotalDisplacements
    file_out = open("growth_remodeling.txt","a")
    file_out.write('%3d %f %f %f %f %f %f %f %f %f %f %f %f %f\n'%(step,euler_x[0,0],growth_remodeling[0,0],growth_remodeling[0,1],growth_remodeling[0,2],growth_remodeling[0,3],growth_remodeling[0,4],growth_remodeling[0,5],growth_remodeling[0,6],growth_remodeling[0,7],growth_remodeling[0,8],growth_remodeling[0,9],growth_remodeling[0,10],growth_remodeling[0,11]))
    file_out.close()