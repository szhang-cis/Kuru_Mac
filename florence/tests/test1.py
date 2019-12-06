# System libraries and interaction with user
import sys, os
# Mathematics libraries
import numpy as np
# Build a path for python to Florence
sys.path.append(os.path.join(os.path.expanduser("~"),"florence"))
#import Florence
from Florence import *

output= open("results.dat","w")

def homogenized_constrained_mixture():
    """A hyperelastic explicit dynamics example using Mooney Rivlin model
        of a cylinder column under compression with cubic (p=3) hexahedral elements
    """
    ProblemPath = PWD(__file__)
    mesh_file = ProblemPath + '/Half_Cylinder.msh'
    # Build mesh with Florence tools
    mesh = Mesh()
    mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex")
#    mesh.GetHighOrderMesh(p=2)
    ndim = mesh.InferSpatialDimension()

    # Define the anisotropic orientations
    fibre_axial = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    fibre_dia1 = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    fibre_dia2 = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    fibre_circ = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    center = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    tangential = np.zeros((mesh.nelem,ndim),dtype=np.float64)
    norma = np.zeros((mesh.nelem),dtype=np.float64)
    directrix = [0.,1.,0.]
    divider = mesh.points[mesh.elements[0,:],0].shape[0]
    for i in range(mesh.nelem):
        center[i,0] = np.sum(mesh.points[mesh.elements[i,:],0])/divider
        center[i,1] = np.sum(mesh.points[mesh.elements[i,:],1])/divider
        center[i,2] = np.sum(mesh.points[mesh.elements[i,:],2])/divider
        tangential[i,:] = np.cross(directrix,center[i,:])
        norma[i] = np.linalg.norm(tangential[i,:])
        tangential[i,:] = tangential[i,:]/norma[i]
        fibre_axial[i,:] = np.multiply(directrix,np.cos(0.)) + np.multiply(tangential[i,:],np.sin(0.))
        fibre_dia1[i,:] = np.multiply(directrix,np.cos(np.pi/4)) + np.multiply(tangential[i,:],np.sin(np.pi/4))
        fibre_dia2[i,:] = np.multiply(directrix,np.cos(np.pi/4)) - np.multiply(tangential[i,:],np.sin(np.pi/4))
        fibre_circ[i,:] = np.multiply(directrix,np.cos(np.pi/2)) + np.multiply(tangential[i,:],np.sin(np.pi/2))

    # Dictionary with multiple fibre directions
    fibre_directions = {}
    fibre_directions['smc'] = fibre_circ
    fibre_directions['co1'] = fibre_circ
    fibre_directions['co2'] = fibre_dia1
    fibre_directions['co3'] = fibre_dia2
    fibre_directions['co4'] = fibre_axial

    # Dictionary with mass fractions
    phi = {}
    phi['ela'] = 0.23
    phi['smc'] = 0.15
    phi['co1'] = 0.062
    phi['co2'] = 0.248
    phi['co3'] = 0.248
    phi['co4'] = 0.062

    # Define hyperelastic material for mesh
    material = ArterialWallMixture(ndim,
                 mu=72.0e3,
                 c1m=15.2e3,
                 c2m=11.4,
                 c1c=1136.0e3,
                 c2c=11.2,
                 kappa=720.0e3,
                 rho=1050.0e-9,
                 anisotropic_orientations=fibre_directions,
                 mass_fraction=phi
                 )
    # Dirichlet Boundary Conditions
    def Dirichlet_Function(mesh, time_step):
        #num_of_points=mesh.point.shape[0], ndim=3 ==> boundary_data(npoin,ndim,ntime) is third order
        boundary_data = np.zeros((mesh.points.shape[0],3, time_step))+np.NAN
        #boolean array with True where Z=0 (X_0)
        Y_0 = np.isclose(mesh.points[:,1],0)
        Y_1 = np.isclose(mesh.points[:,1],mesh.points[:,1].max())
        Z_0 = np.isclose(mesh.points[:,2],0)
        boundary_data[Y_0,1,:] = 0.
        boundary_data[Y_1,1,:] = 0.
        boundary_data[Z_0,2,:] = 0.

        return boundary_data

    def Neumann_Function(mesh, time_step):
        #external_faces=mesh.faces.shape[0]
        boundary_flags = np.zeros((mesh.faces.shape[0], time_step),dtype=np.uint8)
        boundary_data = np.zeros((mesh.faces.shape[0],3, time_step))
        # Force magnitud
        mag = 13.3e-3

        count = 0
        for i in range(mesh.faces.shape[0]):
            coord = mesh.points[mesh.faces[i,:],:]
            coord2 = np.square(coord)
            x2_plus_z2 = coord2[:,0] + coord2[:,2]
            radius = np.sqrt(x2_plus_z2)
            avg = np.sum(radius)/mesh.faces.shape[1]
            if np.isclose(avg,9.295):
                normal = np.sum(coord,axis=0)/mesh.faces.shape[1]
                largo = np.sqrt(np.square(normal[0]) + np.square(normal[2]))
                normal = normal/largo
                mag0 = normal[0]*mag
                mag2 = normal[2]*mag
                boundary_data[i,0,:] = np.linspace(0,mag0,time_step)
                boundary_data[i,2,:] = np.linspace(0,mag2,time_step)
                boundary_flags[i,:] = True
                count = count + 1

	print count

        return boundary_flags, boundary_data

    time_step = 3
    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(Dirichlet_Function, mesh, time_step)
    boundary_condition.SetNeumannCriteria(Neumann_Function, mesh, time_step)
"""
    formulation = DisplacementFormulation(mesh)
    fem_solver = FEMSolver( total_time=1.,
                            number_of_load_increments=time_step,
                            analysis_nature="linear",
                            analysis_type="dynamic",
                            analysis_subtype="implicit",
#                            mass_type="consistent",
                            optimise=False,
                            print_incremental_log=True,
                            memory_store_frequency=1)

    solution = fem_solver.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)

    # Write to paraview
    solution.WriteVTK("test1",quantity=0)
"""

if __name__ == "__main__":
    homogenized_constrained_mixture()



"""
==================    STATIC EXAMPLES     ============================
    # Quasi-static nonlinear monolithic solver
    static_monolithic_solver = FEMSolver(
        number_of_load_increments=time_step,
        analysis_nature="nonlinear",
        analysis_type="static",
        newton_raphson_tolerance=1e-5,
        optimise=True,
        )
    static_monolithic_solver_results = static_monolithic_solver.Solve(formulation=formulation,
        mesh=mesh, material=material, boundary_condition=boundary_condition)

    nonlinear_static_solver = FEMSolver(total_time=60.,
        number_of_load_increments=25,
        analysis_nature="nonlinear",
        analysis_type="static",
        newton_raphson_tolerance=1e-5,
        newton_raphson_solution_tolerance=1e-11,
        optimise=optimise,
        print_incremental_log=True,
        )
    nonlinear_static_results = nonlinear_static_solver.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)

==================    DYNAMIC EXAMPLES     ============================

    # Dynamic nonlinear monolithic solver with Newmark beta
    dynamic_monolithic_solver = FEMSolver(
        total_time=10.,
        number_of_load_increments=time_step,
        analysis_nature="nonlinear",
        analysis_type="dynamic",
        newton_raphson_tolerance=1e-5,
        optimise=True,
        )
    dynamic_monolithic_solver_results = dynamic_monolithic_solver.Solve(formulation=formulation,
        mesh=mesh, material=material, boundary_condition=boundary_condition)

    nonlinear_dynamic_solver = FEMSolver(total_time=60.,
        number_of_load_increments=250,
        analysis_nature="nonlinear",
        analysis_type="dynamic",
        newton_raphson_tolerance=1e-5,
        newton_raphson_solution_tolerance=1e-11,
        optimise=optimise,
        print_incremental_log=True,
        compute_energy_dissipation=True,
        compute_linear_momentum_dissipation=True,
        )
    nonlinear_dynamic_results = nonlinear_dynamic_solver.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)

    implicit_solver = FEMSolver(total_time=.1,
        number_of_load_increments=time_step,
        analysis_nature="nonlinear",
        analysis_type="dynamic",
        analysis_subtype="implicit",
        newton_raphson_tolerance=1e-10,
        optimise=optimise,
        compute_energy_dissipation=True,
        compute_linear_momentum_dissipation=True
        )
    results_implicit = implicit_solver.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)

    explicit_solver_consistent_mass = FEMSolver(total_time=.1,
        number_of_load_increments=time_step*8,
        analysis_nature="nonlinear",
        analysis_type="dynamic",
        analysis_subtype="explicit",
        mass_type="consistent",
        optimise=optimise,
        compute_energy_dissipation=True,
        compute_linear_momentum_dissipation=True
        )
    results_explicit_consistent_mass = explicit_solver_consistent_mass.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)

    explicit_solver_lumped_mass = FEMSolver(total_time=.1,
        number_of_load_increments=time_step*8,
        analysis_nature="nonlinear",
        analysis_type="dynamic",
        analysis_subtype="explicit",
        mass_type="lumped",
        optimise=optimise
        )
    results_explicit_lumped_mass = explicit_solver_lumped_mass.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)
"""
