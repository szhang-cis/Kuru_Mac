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
        direction[elem][3][:]=np.multiply(directrix,np.cos(np.pi/4)) + np.multiply(tangential,np.sin(np.pi/4))
        direction[elem][4][:]=np.multiply(directrix,np.cos(np.pi/4)) - np.multiply(tangential,np.sin(np.pi/4))
        direction[elem][5][:]=tangential

    return direction

def GetRadialStretch(mesh,material):
    print("Computing the Radial deposition stretch")
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    for node in range(mesh.nnode):
        kappa = material.kappa*material.growth_remodeling[node,0]
        mu3D = material.mu3D*material.growth_remodeling[node,0]
        stretch_r = material.deposition_stretch['Radial'][node]
        stretch_t = material.deposition_stretch['Tangential']
        stretch_z = material.deposition_stretch['Axial']
        radius = np.sqrt(mesh.points[node,0]**2+mesh.points[node,2]**2)
        JMAX=21
        for j in range(JMAX):
            if j==20: exit('Maximum number of iteration reach at Newton method at node: '+str(node))
            func = pressure*(1.-(radius-R)/H)+ \
                    mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-\
                    1./3.*(stretch_r**2+stretch_t**2+stretch_z**2)) +\
                    kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.)
            df = 4./3.*mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3)*stretch_r-\
                    2./3.*mu3D*stretch_r**(-5./3.)*(stretch_t*stretch_z)**(-2./3.)*(stretch_r**2-\
                    1./3.*(stretch_r**2+stretch_t**2+stretch_z**2)) +\
                    kappa*(2.*stretch_r*stretch_t**2*stretch_z**2-stretch_t*stretch_z)
            dx = func/df
            stretch_r -= dx
            if np.absolute(dx)<1.e-6:
                material.deposition_stretch['Radial'][node] = stretch_r
                break

    print("Radial deposition stretch is ready")

def GetTangentialPenalty(mesh,material):
    print("Computing the Tangential penalty mu2D")
    pressure = 13.3322e3
    R = 0.010
    H = 0.00141
    for node in range(mesh.nnode):
        kappa = material.kappa*material.growth_remodeling[node,0]
        mu3D = material.mu3D*material.growth_remodeling[node,0]
        stretch_r = material.deposition_stretch['Radial'][node]
        stretch_t = material.deposition_stretch['Tangential']
        stretch_z = material.deposition_stretch['Axial']
        k1m = material.k1m*material.growth_remodeling[node,1]
        k2m = material.k2m
        s_act = 54.0e3*growth_remodeling[node,1]
        den0 = 1050.0
        lambda_m = 1.4
        lambda_0 = 0.8
        lambda_a = 1.0
        stretch_m = material.deposition_stretch['Muscle']
        k1c1 = material.k1c*material.growth_remodeling[node,2]
        k1c2 = material.k1c*material.growth_remodeling[node,3]
        k1c3 = material.k1c*material.growth_remodeling[node,4]
        k2c = material.k2c
        stretch_c = material.deposition_stretch['Collagen']
        #radius = np.sqrt(mesh.points[node,0]**2+mesh.points[node,2]**2)
        mu2D = (pressure*R/H - \
                kappa*stretch_r*stretch_t*stretch_z*(stretch_r*stretch_t*stretch_z-1.) - \
                mu3D*(stretch_r*stretch_t*stretch_z)**(-2./3.)*(stretch_t**2 - \
                1./3.*(stretch_r**2+stretch_t**2+stretch_z**2)) - \
                2.*k1m*(stretch_m**2-1.)*np.exp(k2m*(stretch_m**2-1.)**2)*stretch_m**2 - \
                s_act*(1.-((lambda_m-lambda_a)/(lambda_m-lambda_0))**2)/den0 - \
                2.*(k1c1+k1c2/2.+k1c3/2.)*(stretch_c**2-1.)*np.exp(k2c*(stretch_c**2-1.)**2)*stretch_c**2)/\
                (stretch_t**2-1./(stretch_t*stretch_z)**2)
        material.mu2D[node] = mu2D/material.growth_remodeling[node,0]
    print("Tangential penalty mu2D is ready")

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
    den0_e = 241.5
    # Fribres turnover and remodeling
    turnover = 101.0
    gain = 0.09/turnover
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
mesh_file = ProblemPath + '/QuarterRingShell.msh'

#===============  MESH PROCESING  ==========================
# Build mesh with Kuru tools from GMSH mesh
mesh = Mesh()
mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="quad",read_surface_info=True,read_curve_info=True)
ndim = mesh.InferSpatialDimension()
mesh.GetHighOrderMesh(p=1)
#Boolean arrays for boundary condition in Dirichlet
BottomLine = np.zeros(mesh.nnode,dtype=bool)
TopLine = np.zeros(mesh.nnode,dtype=bool)
Symmetry_Z = np.zeros(mesh.nnode,dtype=bool)
Symmetry_X = np.zeros(mesh.nnode,dtype=bool)
#Boolean array for boundary condition in Neumann mesh.faces[id,node]
InnerSurface = np.zeros(mesh.faces.shape[0],dtype=bool)
InnerFaces = []
for idface in range(mesh.faces.shape[0]):
    if mesh.face_to_surface[idface] == 4:
        InnerFaces.append([int(i) for i in mesh.faces[idface,:]])
        InnerSurface[idface] = True

for idedge in range(mesh.edges.shape[0]):
    if mesh.edge_to_curve[idedge] == 2:
        BottomLine[mesh.edges[idedge,:]] = True
    elif mesh.edge_to_curve[idedge] == 3:
        TopLine[mesh.edges[idedge,:]] = True
    elif mesh.edge_to_curve[idedge] == 0:
        Symmetry_Z[mesh.edges[idedge,:]] = True
    elif mesh.edge_to_curve[idedge] == 1:
        Symmetry_X[mesh.edges[idedge,:]] = True

DirichletBoundary = {}
DirichletBoundary['Bottom'] = BottomLine
DirichletBoundary['Top'] = TopLine
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

# Total initial density
#mu2D = np.zeros((mesh.nnode),dtype=np.float64)
#mu2D = 0.0

# Define hyperelastic material for mesh
material = ArterialWallMixture(ndim,
            mu=72.0,
            kappa=72.0*20.0,
            k1m=7.6,
            k2m=11.4,
            k1c=568.0,
            k2c=11.2,
            anisotropic_orientations=fibre_direction,
            field_variables=field_variables)

# kappa/mu=20  => nu=0.475 (Poisson's ratio)
# kappa/mu=33  => nu=0.485 (Poisson's ratio)
# kappa/mu=100 => nu=0.495 (Poisson's ratio)

# Homogenization of stresses to a thin-walled case
#GetRadialStretch(mesh,material)
#print(material.deposition_stretch['Radial'][:])
#GetTangentialPenalty(mesh,material)
#print(material.mu2D[:])

#==================  FORMULATION  =========================
formulation = DisplacementFormulationShell(mesh)

#===============  BOUNDARY CONDITIONS  ====================
# Dirichlet Boundary Conditions
def Dirichlet_Function(mesh, DirichletBoundary):
    boundary_data = np.zeros((mesh.nnode, 5))+np.NAN
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
fem_solver_h = ShellSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=True,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

fem_solver_gr = FEMSolver(analysis_nature="nonlinear",
                       analysis_type="static",
                       break_at_stagnation=False,
                       maximum_iteration_for_newton_raphson=50,
                       optimise=True,
                       print_incremental_log=True,
                       has_moving_boundary=True,
                       number_of_load_increments=1)

#=================  HOMEOSTATIC SOLUTION  =======================
print('=====================================')
print('===    COMPUTE HOMEOSTATIC STEP   ===')
print('=====================================')
# Call the solver for Homeostatic computation
solution = fem_solver_h.Solve(formulation=formulation, mesh=mesh,
    material=material, boundary_condition=boundary_condition)
"""
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
file_out.write('%3d %f %f %f %f %f %f %f %f %f %f %f %f %f \n'%(0,1000.*np.sqrt(euler_x[0,0]**2+euler_x[0,2]**2),field_variables[0,11],field_variables[0,12],field_variables[0,13],field_variables[0,14],field_variables[0,15],field_variables[0,16],field_variables[0,17],field_variables[0,18],field_variables[0,19],field_variables[0,20],field_variables[0,21],field_variables[0,22]))
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
    field_variables[:,11:23] = GetGrowthRemodeling(time,Delta_t,mesh,field_variables[:,11:23],Stress_H,FibreStress,Softness)
    material.field_variables = field_variables
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
    file_out.write('%3d %f %f %f %f %f %f %f %f %f %f %f %f %f \n'%(0,1000.*np.sqrt(euler_x[0,0]**2+euler_x[0,2]**2),field_variables[0,11],field_variables[0,12],field_variables[0,13],field_variables[0,14],field_variables[0,15],field_variables[0,16],field_variables[0,17],field_variables[0,18],field_variables[0,19],field_variables[0,20],field_variables[0,21],field_variables[0,22]))
    file_out.close()
"""
