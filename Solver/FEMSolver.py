from __future__ import print_function
import gc, os, sys
#import multiprocessing
from copy import deepcopy
from warnings import warn
from time import time
import numpy as np
import numpy.linalg as la
from numpy.linalg import norm
#import scipy as sp
from Kuru.Utils import insensitive

from Kuru.FiniteElements.Assembly import Assemble #, AssembleExplicit
from Kuru.PostProcessing import *
from Kuru.Solver import LinearSolver
#from Florence.TimeIntegrators import LinearImplicitStructuralDynamicIntegrator
#from Florence.TimeIntegrators import NonlinearImplicitStructuralDynamicIntegrator
#from Florence.TimeIntegrators import ExplicitStructuralDynamicIntegrator
#from .LaplacianSolver import LaplacianSolver
from Kuru import Mesh


__all__ = ["FEMSolver"]

class FEMSolver(object):
    """Solver for linear and non-linear finite elements.
        This class is fundamentally different from the LinearSolver, as linear solver
        specifically deals with the solution of linear system of equations, whereas FEM
        solver is essentially responsible for building the discretised linear, linearised
        and nonlinear systems arising from variational formulations
    """

    def __init__(self,
        has_low_level_dispatcher=False,
        optimise=False,
        analysis_type="static",
        analysis_nature="nonlinear",
        analysis_subtype="implicit",
        #linearised_electromechanics_solver="traction_based",
        is_geometrically_linearised=False,
        requires_geometry_update=True,
        requires_line_search=False,
        requires_arc_length=False,
        has_moving_boundary=False,
        has_prestress=True,
        has_growth_remodeling=False,
        number_of_load_increments=1,
        number_of_time_increments=1,
        time_factor=None,
        load_factor=None,
        newton_raphson_tolerance=1.0e-6,
        newton_raphson_solution_tolerance=None,
        maximum_iteration_for_newton_raphson=50,
        nonlinear_iterative_technique="newton_raphson",
        reduce_quadrature_for_quads_hexes=True,
        add_self_weight=False,
        mass_type=None,
        compute_mesh_qualities=False,
        force_not_computing_mesh_qualities=False,
        parallelise=False,
        ncpu=None,
        parallel_model=None,
        memory_model="shared",
        platform="cpu",
        backend="opencl",
        print_incremental_log=False,
        save_incremental_solution=False,
        incremental_solution_filename=None,
        incremental_solution_save_frequency=50,
        break_at_increment=-1,
        break_at_stagnation=True,
        include_physical_damping=False,
        damping_factor=0.1,
        compute_energy=False,
        compute_energy_dissipation=False,
        compute_linear_momentum_dissipation=False,
        total_time=1.,
        user_defined_break_func=None,
        user_defined_stop_func=None,
        save_results=True,
        memory_store_frequency=1,
        has_contact=False,
        activate_explicit_multigrid=False,
        recompute_sparsity_pattern=True,
        squeeze_sparsity_pattern=False,
        ensure_structured_grid=False,
        do_not_reset=True,
        report_log_level=2):

        # ASSUME TRUE IF AT LEAST ONE IS TRUE
        if has_low_level_dispatcher != optimise:
            has_low_level_dispatcher = True
        self.has_low_level_dispatcher = has_low_level_dispatcher

        self.analysis_nature = analysis_nature
        self.analysis_type = analysis_type
        self.analysis_subtype = analysis_subtype

        #self.linearised_electromechanics_solver=linearised_electromechanics_solver # "traction_based" or "potential_based"
        self.is_geometrically_linearised = is_geometrically_linearised
        self.requires_geometry_update = requires_geometry_update
        self.requires_line_search = requires_line_search
        self.requires_arc_length = requires_arc_length
        self.has_moving_boundary = has_moving_boundary
        self.has_prestress = has_prestress
        self.has_growth_remodeling = has_growth_remodeling
        self.is_mass_computed = False
        self.mass_type = mass_type # "consistent" or "lumped"
        self.total_time = float(total_time)
        self.number_of_time_increments = number_of_time_increments
        self.number_of_load_increments = number_of_load_increments
        self.time_factor = time_factor
        self.load_factor = load_factor
        self.save_results = save_results
        # SAVE AT EVERY N TIME STEP WHERE N=save_frequency
        self.save_frequency = int(memory_store_frequency)
        self.incremental_solution_save_frequency = incremental_solution_save_frequency

        self.newton_raphson_tolerance = newton_raphson_tolerance
        self.newton_raphson_solution_tolerance = newton_raphson_solution_tolerance
        self.maximum_iteration_for_newton_raphson = maximum_iteration_for_newton_raphson
        self.newton_raphson_failed_to_converge = False
        self.NRConvergence = None
        self.nonlinear_iterative_technique = nonlinear_iterative_technique
        self.include_physical_damping = include_physical_damping
        self.damping_factor = damping_factor
        self.add_self_weight = add_self_weight
        self.reduce_quadrature_for_quads_hexes = reduce_quadrature_for_quads_hexes

        self.compute_energy = compute_energy
        self.compute_energy_dissipation = compute_energy_dissipation
        self.compute_linear_momentum_dissipation = compute_linear_momentum_dissipation

        self.compute_mesh_qualities = compute_mesh_qualities
        self.force_not_computing_mesh_qualities = force_not_computing_mesh_qualities
        self.is_scaled_jacobian_computed = False

        self.vectorise = True
        self.parallel = parallelise
        self.no_of_cpu_cores = ncpu
        self.parallel_model = parallel_model # "pool", "context_manager", "queue", "joblib"
        self.memory_model = memory_model
        self.platform = platform
        self.is_partitioned = False
        self.backend = backend
        self.debug = False

        self.print_incremental_log = print_incremental_log
        self.save_incremental_solution = save_incremental_solution
        self.incremental_solution_filename = incremental_solution_filename
        self.break_at_increment = break_at_increment
        self.break_at_stagnation = break_at_stagnation
        self.user_defined_break_func = user_defined_break_func
        self.user_defined_stop_func = user_defined_stop_func

        self.has_contact = has_contact
        self.activate_explicit_multigrid = activate_explicit_multigrid

        self.recompute_sparsity_pattern = recompute_sparsity_pattern
        self.squeeze_sparsity_pattern = squeeze_sparsity_pattern
        self.ensure_structured_grid = ensure_structured_grid
        self.is_sparsity_pattern_computed = False

        self.fem_timer = 0.
        self.assembly_time = 0.
        self.is_dask_scheduler_initialised = False

        if self.newton_raphson_solution_tolerance is None:
            self.newton_raphson_solution_tolerance = 10.*self.newton_raphson_tolerance

        self.report_log_level = report_log_level

        # STORE A COPY OF SELF AT THE START TO RESET TO AT THE END
        self.__save_state__()
        # FOR INTERNAL PURPOSES WHEN WE DO NOT WANT TO REST
        self.do_not_reset = do_not_reset

    def __save_state__(self):
        self.__initialdict__ = deepcopy(self.__dict__)

    def __checkdata__(self, materials, boundary_condition, formulation, mesh, function_spaces, 
            solver, contact_formulation=None, growth_remodeling=None):
        """Checks the state of data for FEMSolver"""

        # INITIAL CHECKS
        ###########################################################################
        if mesh is None:
            raise ValueError("No mesh detected for the analysis")
        elif not isinstance(mesh,Mesh):
            raise ValueError("mesh has to be an instance of Florence.Mesh")
        if boundary_condition is None:
            raise ValueError("No boundary conditions detected for the analysis")
        if isinstance(materials,tuple):
            if len(materials)==0:
                raise ValueError("No material model chosen for the analysis")
        else:
            raise ValueError("The material should be insterted in a tuple")
        if formulation is None:
            raise ValueError("No variational formulation specified")
        # MONKEY PATCH CONTACT FORMULATION
        if contact_formulation is not None:
            self.contact_formulation = contact_formulation
            self.has_contact = True
        else:
            self.contact_formulation = None
            self.has_contact = False
        # GROWTH REMODELING FORMULATION
        if growth_remodeling is not None:
            self.has_growth_remodeling = True
        else:
            self.has_growth_remodeling = False

        # GET FUNCTION SPACES FROM THE FORMULATION
        if function_spaces is None:
            if formulation.function_spaces is None:
                raise ValueError("No interpolation functions specified")
            else:
                function_spaces = formulation.function_spaces

        # CHECK IF A SOLVER IS SPECIFIED
        if solver is None:
            solver = LinearSolver(linear_solver="direct", linear_solver_type="umfpack", geometric_discretisation=mesh.element_type)

        # CHECK MESH DIMENSION BEFORE ANY FURTHER ANALYSIS
        if mesh.points.shape[1] == 3 and (mesh.element_type == "tri" or mesh.element_type == "quad"):
            if np.allclose(mesh.points[:,2],0.):
                mesh.points = np.ascontiguousarray(mesh.points[:,:2])
            elif np.allclose(mesh.points[:,1],0.):
                mesh.points = np.ascontiguousarray(mesh.points[:,[0,2]])
            elif np.allclose(mesh.points[:,0],0.):
                mesh.points = np.ascontiguousarray(mesh.points[:,1:])
            else:
                warn("Mesh spatial dimensionality is not correct for 2D. I will ignore Z coordinates")
                mesh.points = np.ascontiguousarray(mesh.points[:,:2])

        # TURN OFF PARALLELISM IF NOT AVAILABLE
        if self.parallel:
            if formulation.fields == "mechanics" and self.has_low_level_dispatcher:
                if self.parallel_model is None:
                    if self.analysis_type == "dynamic" and self.analysis_subtype == "explicit":
                        self.parallel_model = "context_manager"
                    else:
                        self.parallel_model = "pool"
            else:
                warn("Parallelism cannot be activated right now")
                self.parallel = False
                self.no_of_cpu_cores = 1

        if self.parallel:
            if self.no_of_cpu_cores is None:
                self.no_of_cpu_cores = multiprocessing.cpu_count()
            # PARTITION THE MESH AND PREPARE
            self.PartitionMeshForParallelFEM(mesh,self.no_of_cpu_cores,formulation.nvar)

        for imat in range(len(materials)):
            if materials[imat].ndim != mesh.InferSpatialDimension():
                # THIS HAS TO BE AN ERROR BECAUSE OF THE DYNAMIC NATURE OF MATERIAL
                raise ValueError("Material {} model and mesh are incompatible. Change the dimensionality of the material".format(imat))
            if materials[imat].nvar != formulation.nvar:
                # THIS HAS TO BE AN ERROR BECAUSE OF THE DYNAMIC NATURE OF MATERIAL
                raise ValueError("Material {} model and formulation are incompatible.".format(imat))
        ###########################################################################
#        if material.mtype == "LinearElastic" and self.number_of_load_increments > 1 and self.analysis_type=="static":
#            warn("Can not solve a linear elastic model in multiple step. "
#                "The number of load increments is going to be set to one (1). "
#                "Use IncrementalLinearElastic model for incremental soluiton.")
#            self.number_of_load_increments = 1

        self.has_prestress = False
#        if "nonlinear" not in insensitive(self.analysis_nature) and formulation.fields == "mechanics":
#            # RUN THE SIMULATION WITHIN A NONLINEAR ROUTINE
#            if material.mtype != "IncrementalLinearElastic" and \
#                material.mtype != "LinearElastic" and material.mtype != "TranservselyIsotropicLinearElastic":
#                self.has_prestress = True
#            else:
#                self.has_prestress = False


        ##############################################################################
#        if "Explicit" in material.mtype and self.analysis_subtype == "implicit":
#                raise ValueError("Incorrect material model ({}) used for implicit analysis".format(material.mtype))
        if self.mass_type is None:
            if self.analysis_subtype == "explicit":
                self.mass_type = "lumped"
            else:
                self.mass_type = "consistent"
        if self.analysis_type == "dynamic" and self.analysis_subtype == "implicit" and self.mass_type == "lumped":
            warn("Cannot use lumped mass matrix for implicit analysis. Changing to consistent mass matrix")
            self.mass_type = "consistent"
        if self.analysis_type == "static":
            if self.save_frequency != 1:
                warn("memory_store_frequency must be one")
                self.save_frequency = 1
        if self.analysis_type == "dynamic" and self.analysis_subtype=="implicit" and self.analysis_nature=="linear":
            if self.save_frequency != 1:
                warn("memory_store_frequency must be one")
                self.save_frequency = 1
        if self.analysis_type == "dynamic":
            if self.number_of_load_increments < self.save_frequency:
                raise ValueError("Number of load increments cannot be less than memory store frequency")
            if self.number_of_load_increments < 3:
                warn("Number of load steps={} is excessively low for dynamic analysis. "
                    "I will increase it by a bit".format(self.number_of_load_increments))
                self.number_of_load_increments = 3
        ##############################################################################


        # GEOMETRY UPDATE FLAGS
        ###########################################################################
        self.requires_geometry_update = True
        if formulation.fields == "mechanics":
            if self.analysis_nature == "nonlinear":
                self.requires_geometry_update = True
            else:
                self.requires_geometry_update = False

        # CHECK IF MATERIAL MODEL AND ANALYSIS TYPE ARE COMPATIBLE
        #############################################################################
        for imat in range(len(materials)):
            if materials[imat].fields != formulation.fields:
                raise RuntimeError("Incompatible material model and formulation type")
            if materials[imat].is_transversely_isotropic or materials[imat].is_anisotropic:
                if materials[imat].anisotropic_orientations is None:
                    materials[imat].GetFibresOrientation(mesh)
            materials[imat].ConnectivityOfMaterial(mesh)
        ##############################################################################

        ##############################################################################
        if not self.force_not_computing_mesh_qualities:
            if boundary_condition.boundary_type == "nurbs":
                self.compute_mesh_qualities = True
        ##############################################################################

        ##############################################################################
        if self.load_factor is not None:
            self.load_factor = np.array(self.load_factor).ravel()
            if self.load_factor.shape[0] != self.number_of_load_increments:
                raise ValueError("Supplied load factor should have the same length as the number of load increments")
            if not np.isclose(self.load_factor.sum(),1.0):
                raise ValueError("Load factor should sum up to one")
        ##############################################################################

        ##############################################################################
        for imat in range(len(materials)):
            if self.analysis_type == "dynamic" and materials[imat].rho is None:
                raise ValueError("Material does not seem to have density. Mass matrix cannot be computed")
            if self.analysis_type == "static" and materials[imat].rho is None and self.has_low_level_dispatcher:
                # FOR LOW-LEVEL DISPATCHER
                materials[imat].rho = 1.0
            # check material load factor against solver load increments
            if materials[imat].load_factor is not None:
                materials[imat].load_factor = np.array(materials[imat].load_factor).ravel()
                if materials[imat].load_factor.shape[0] != self.number_of_load_increments:
                    raise ValueError("Supplied material load factor should have the same length as the number of load increments")
                if not np.isclose(materials[imat].load_factor.sum(),1.0):
                    raise ValueError("Load factor (material {}) should sum up to one".format(imat))
        ##############################################################################

        ##############################################################################
        if self.include_physical_damping and self.compute_energy_dissipation:
            warn("Energy is not going to be preserved due to physical damping")
        ##############################################################################

        ##############################################################################
        if self.analysis_type == "static" and self.has_contact is True:
            warn("Contact formulation does not get activated under static problems")
        if self.analysis_type == "dynamic" and self.analysis_nature == "nonlinear" \
            and self.analysis_subtype == "implicit" and self.has_contact is True:
            raise ValueError("Implicit contact formulation for nonlinear implicit dynamics is not supported")
        ##############################################################################

        ##############################################################################
        for imat in range(len(materials)):
            if materials[imat].has_state_variables:
                elem = materials[imat].element_set[0]
                materials[imat].MappingStateVariables(mesh,function_spaces[0],elem=elem)
            # AT THE MOMENT ALL HESSIANS SEEMINGLY HAVE THE SAME SIGNATURE SO THIS IS O.K.
            try:
                F = np.eye(materials[imat].ndim,materials[imat].ndim)[None,:,:]
                materials[imat].Hessian(KinematicMeasures(F,self.analysis_nature))
            except TypeError:
                # CATCH ONLY TypeError. OTHER MATERIAL CONSTANT RELATED ERRORS ARE SELF EXPLANATORY
                raise ValueError("Material constants for {} does not seem correct".format(material.mtype))
        ##############################################################################

        # CHANGE MESH DATA TYPE
        mesh.ChangeType()

        # ASSIGN ANALYSIS PARAMTER TO BOUNDARY CONDITION
        boundary_condition.analysis_type = self.analysis_type
        boundary_condition.analysis_nature = self.analysis_nature
        ##############################################################################
        if self.analysis_type == "dynamic" and self.analysis_subtype != "explicit" and\
            self.analysis_nature == "linear":
            boundary_condition.ConvertStaticsToDynamics(mesh, self.number_of_load_increments)
        ##############################################################################

        ##############################################################################
        if boundary_condition.dirichlet_flags is not None:
            # SPECIFIC CHECKS TO AVOID CONFUSING ERRORS OCCURING DOWN STREAM
            ndim = mesh.InferSpatialDimension()
            if boundary_condition.dirichlet_flags.ndim == 2:
                if boundary_condition.dirichlet_flags.shape[1] != formulation.nvar:
                    raise ValueError("Essential boundary conditions are not imposed correctly")
            if formulation.fields == "mechanics":
                if boundary_condition.dirichlet_flags.shape[1] != ndim:
                    raise ValueError("Essential boundary conditions are not imposed correctly")

            if boundary_condition.has_step_wise_dirichlet_loading or \
                boundary_condition.has_step_wise_neumann_loading:
                if self.analysis_type != "dynamic":
                    raise ValueError("User defined step loading only applies for dynamic problems")
                if formulation.fields != "mechanics":
                    raise NotImplementedError("User defined step loading under this setting is not supported")
                if self.analysis_nature == "linear":
                    raise NotImplementedError("User defined step loading under this setting is not supported")
        ##############################################################################

        return function_spaces, solver

    def __makeoutput__(self, mesh, TotalDisp, formulation=None, function_spaces=None, materials=None,
            gr_variables=None):

        post_process = PostProcess(formulation.ndim,formulation.nvar)
        post_process.SetBases(postdomain=function_spaces[1], domain=function_spaces[0], boundary=None)
        post_process.SetAnalysis(analysis_type=self.analysis_type,
            analysis_nature=self.analysis_nature)
        post_process.SetMesh(mesh)
        post_process.SetSolution(TotalDisp)
        post_process.SetFormulation(formulation)
        post_process.SetMaterial(materials)
        post_process.SetFEMSolver(self)
        if gr_variables is not None:
            post_process.SetGrowthRemodeling(gr_variables)
        else:
            post_process.SetGrowthRemodeling(None)

        if self.compute_mesh_qualities and self.is_scaled_jacobian_computed is False:
            # COMPUTE QUALITY MEASURES
            post_process.ScaledFF, post_process.ScaledHH, post_process.ScaledJacobian \
                = post_process.MeshQualityMeasures(mesh,TotalDisp,False,False)[1:4]
        elif self.is_scaled_jacobian_computed:
            post_process.ScaledJacobian=self.ScaledJacobian
            post_process.ScaledHH=self.ScaledHH
            post_process.ScaledFF=self.ScaledFF

        if self.analysis_nature == "nonlinear":
            post_process.newton_raphson_convergence = self.NRConvergence

        if self.analysis_type == "dynamic":
            if self.compute_energy_dissipation:
                post_process.energy_dissipation = formulation.energy_dissipation
                post_process.internal_energy = formulation.internal_energy
                post_process.kinetic_energy = formulation.kinetic_energy
                post_process.external_energy = formulation.external_energy
            if self.compute_linear_momentum_dissipation:
                post_process.power_dissipation = formulation.power_dissipation
                post_process.internal_power = formulation.internal_power
                post_process.kinetic_power = formulation.kinetic_power
                post_process.external_power = formulation.external_power


        post_process.assembly_time = self.assembly_time

        # CLOSE DASK CLIENT
        self.CloseDaskDistributedClient()

        # AT THE END, WE CALL THE __reset_state__ TO RESET TO INITIAL STATE.
        # THIS WAY WE CLEAR MONKEY PATCHED AND OTHER DATA STORED DURING RUN TIME
        if self.do_not_reset is False:
            self.__reset_state__()

        return post_process

    def Solve(self, formulation=None, mesh=None,
        materials=None, boundary_condition=None,
        function_spaces=None, solver=None,
        contact_formulation=None, growth_remodeling=None,
        Eulerx=None):
        """Main solution routine for FEMSolver """

        # LOG LEVEL CONTROLLER
        if self.report_log_level == 0:
            sys.stdout = open(os.devnull, "w")

        # CHECK DATA CONSISTENCY
        #---------------------------------------------------------------------------#
        function_spaces, solver = self.__checkdata__(materials, boundary_condition, formulation, 
                mesh, function_spaces, solver, contact_formulation=contact_formulation,
                growth_remodeling=growth_remodeling)

        # PRINT INFO
        #---------------------------------------------------------------------------#
        self.PrintPreAnalysisInfo(mesh, formulation)
        #---------------------------------------------------------------------------#

        # COMPUTE SPARSITY PATTERN
        #---------------------------------------------------------------------------#
        self.ComputeSparsityFEM(mesh, formulation)
        #---------------------------------------------------------------------------#

        # LAUNCH DISTRIBUTED SCHEDULER
        #---------------------------------------------------------------------------#
        if self.parallel and self.parallel_model=="dask" and self.is_dask_scheduler_initialised is False:
            self.LaunchDaskDistributedClient()
        #---------------------------------------------------------------------------#

        # INITIATE DATA FOR THE ANALYSIS
        NodalForces, Residual = np.zeros((mesh.points.shape[0]*formulation.nvar,1),dtype=np.float64),\
            np.zeros((mesh.points.shape[0]*formulation.nvar,1),dtype=np.float64)
        # SET NON-LINEAR PARAMETERS
        if not self.has_growth_remodeling:
            self.NRConvergence = { 'Increment_'+str(Increment) : [] for Increment in range(\
                    self.number_of_load_increments) }
        else:
            self.NRConvergence = { 'Increment_'+str(Increment) : [] for Increment in range(\
                    self.number_of_time_increments) }

        # ALLOCATE FOR SOLUTION FIELDS
        if self.has_growth_remodeling:
            if self.save_frequency == 1:
                TotalDisp = np.zeros((mesh.points.shape[0],formulation.nvar,
                    self.number_of_time_increments),dtype=np.float64)
            else:
                TotalDisp = np.zeros((mesh.points.shape[0],formulation.nvar,
                    int(self.number_of_time_increments/self.save_frequency)),dtype=np.float64)
        else:
            if self.save_frequency == 1:
                TotalDisp = np.zeros((mesh.points.shape[0],formulation.nvar,
                    self.number_of_load_increments),dtype=np.float64)
            else:
                TotalDisp = np.zeros((mesh.points.shape[0],formulation.nvar,
                    int(self.number_of_load_increments/self.save_frequency)),dtype=np.float64)

        # PRE-ASSEMBLY
        print('Assembling the system and acquiring neccessary information for the analysis...')
        tAssembly=time()

        # APPLY DIRICHELT BOUNDARY CONDITIONS AND GET DIRICHLET RELATED FORCES
        boundary_condition.GetDirichletBoundaryConditions(formulation, mesh, materials, solver, self)

        # ALLOCATE FOR GEOMETRY - GetDirichletBoundaryConditions CHANGES THE MESH
        # SO EULERX SHOULD BE ALLOCATED AFTERWARDS
        Eulerx = np.copy(mesh.points)

        # FIND PURE NEUMANN (EXTERNAL) NODAL FORCE VECTOR
        NeumannForces = boundary_condition.ComputeNeumannForces(mesh, materials, function_spaces,
            compute_traction_forces=True, compute_body_forces=self.add_self_weight)

        # ADOPT A DIFFERENT PATH FOR INCREMENTAL LINEAR ELASTICITY
        if formulation.fields == "mechanics" and self.analysis_nature != "nonlinear":
            if self.analysis_type == "static":
                # DISPATCH INCREMENTAL LINEAR ELASTICITY SOLVER
                TotalDisp = self.IncrementalLinearElasticitySolver(function_spaces, formulation, mesh, materials,
                    boundary_condition, solver, TotalDisp, Eulerx, NeumannForces)
                return self.__makeoutput__(mesh, TotalDisp, formulation, function_spaces, material)

        # ASSEMBLE STIFFNESS MATRIX AND TRACTION FORCES FOR THE FIRST TIME (INTERNAL ENERGY)
        if self.analysis_type == "static":
            K, TractionForces, _ = Assemble(self, function_spaces, formulation, mesh, materials, boundary_condition, Eulerx)
        else:
            if self.reduce_quadrature_for_quads_hexes:
                fspace = function_spaces[0] if (mesh.element_type=="hex" or mesh.element_type=="quad") else function_spaces[1]
            else:
                fspace = function_spaces[1]
            # COMPUTE CONSTANT PART OF MASS MATRIX
            formulation.GetConstantMassIntegrand(fspace, material)

            if self.analysis_subtype != "explicit":
                # COMPUTE BOTH STIFFNESS AND MASS USING HIGHER QUADRATURE RULE
                K, TractionForces, M = Assemble(self, fspace, formulation, mesh, material, boundary_condition, Eulerx)
            else:
                # lmesh = mesh.ConvertToLinearMesh()
                # COMPUTE BOTH STIFFNESS AND MASS USING HIGHER QUADRATURE RULE
                TractionForces, _, M = AssembleExplicit(self, fspace, formulation, mesh, material, Eulerx)

        if self.analysis_nature == 'nonlinear':
            print('Finished all pre-processing stage. Time elapsed was', time()-tAssembly, 'seconds')
        else:
            print('Finished the assembly stage. Time elapsed was', time()-tAssembly, 'seconds')

        if self.analysis_type != 'static':
            if self.analysis_subtype != "explicit":
                if self.analysis_nature == "nonlinear":
                    structural_integrator = NonlinearImplicitStructuralDynamicIntegrator()
                    TotalDisp = structural_integrator.Solver(function_spaces, formulation, solver,
                        K, M, NeumannForces, NodalForces, Residual,
                        mesh, TotalDisp, Eulerx, material, boundary_condition, self)
                elif self.analysis_nature == "linear":
                    boundary_condition.ConvertStaticsToDynamics(mesh, self.number_of_load_increments)
                    structural_integrator = LinearImplicitStructuralDynamicIntegrator()
                    TotalDisp = structural_integrator.Solver(function_spaces, formulation, solver,
                        K, M, NeumannForces, NodalForces, Residual,
                        mesh, TotalDisp, Eulerx, material, boundary_condition, self)
            else:
                structural_integrator = ExplicitStructuralDynamicIntegrator()
                TotalDisp = structural_integrator.Solver(function_spaces, formulation, solver,
                    TractionForces, M, NeumannForces, NodalForces, Residual,
                    mesh, TotalDisp, Eulerx, material, boundary_condition, self)

        else:
            if self.has_growth_remodeling:
                TotalDisp,GRVariables = growth_remodeling.Solver(function_spaces, formulation, 
                        solver, K, NeumannForces, NodalForces, Residual, mesh, TotalDisp, Eulerx, 
                        materials, boundary_condition, self)
                return self.__makeoutput__(mesh, TotalDisp, formulation, function_spaces, materials,
                        gr_variables=GRVariables)
            else:
                if self.nonlinear_iterative_technique == "newton_raphson" or \
                    self.nonlinear_iterative_technique == "modified_newton_raphson" or \
                    self.nonlinear_iterative_technique == "line_search_newton_raphson" or \
                    self.nonlinear_iterative_technique == "snes":
                    TotalDisp = self.StaticSolver(function_spaces, formulation, solver,
                        K, NeumannForces, NodalForces, Residual,
                        mesh, TotalDisp, Eulerx, materials, boundary_condition)
                elif self.nonlinear_iterative_technique == "arc_length":
                    from FEMSolverArcLength import StaticSolverArcLength
                    TotalDisp = StaticSolverArcLength(self,function_spaces, formulation, solver,
                        K, NeumannForces, NodalForces, Residual,
                        mesh, TotalDisp, Eulerx, materials, boundary_condition)
                else:
                    raise RuntimeError("Iterative technique for nonlinear solver not understood")

        # LOG LEVEL CONTROLLER
        if self.report_log_level == 0:
            sys.stdout = sys.__stdout__

        return self.__makeoutput__(mesh, TotalDisp, formulation, function_spaces, materials)

    def StaticSolver(self, function_spaces, formulation, solver, K,
            NeumannForces, NodalForces, Residual,
            mesh, TotalDisp, Eulerx, materials, boundary_condition):

        LoadIncrement = self.number_of_load_increments
        LoadFactor = 1./LoadIncrement
        AppliedDirichletInc = np.zeros(boundary_condition.applied_dirichlet.shape[0],dtype=np.float64)

        LoadFactorInc = 0.0
        for Increment in range(LoadIncrement):

            # CHECK ADAPTIVE LOAD FACTOR
            if self.load_factor is not None:
                LoadFactor = self.load_factor[Increment]

            LoadFactorInc += LoadFactor
            # APPLY ROBIN BOUNDARY CONDITIONS - STIFFNESS(_) AND FORCES
            boundary_condition.pressure_increment = LoadFactorInc
            K, RobinForces = boundary_condition.ComputeRobinForces(mesh, materials, function_spaces,
                self, Eulerx, K, np.zeros_like(Residual))
            # APPLY NEUMANN BOUNDARY CONDITIONS
            DeltaF = LoadFactor*NeumannForces
            NodalForces += DeltaF
            # OBRTAIN INCREMENTAL RESIDUAL - CONTRIBUTION FROM BOTH NEUMANN AND DIRICHLET
            Residual = -boundary_condition.ApplyDirichletGetReducedMatrices(K,np.zeros_like(Residual),
                boundary_condition.applied_dirichlet,LoadFactor=LoadFactor,only_residual=True)
            Residual += RobinForces - DeltaF
            # GET THE INCREMENTAL DISPLACEMENT
            AppliedDirichletInc = LoadFactor*boundary_condition.applied_dirichlet

            t_increment = time()

            # LET NORM OF THE FIRST RESIDUAL BE THE NORM WITH RESPECT TO WHICH WE
            # HAVE TO CHECK THE CONVERGENCE OF NEWTON RAPHSON. TYPICALLY THIS IS
            # NORM OF NODAL FORCES
            if Increment==0:
                self.NormForces = np.linalg.norm(Residual)
                # AVOID DIVISION BY ZERO
                if np.isclose(self.NormForces,0.0):
                    self.NormForces = 1e-14

            self.norm_residual = np.linalg.norm(Residual)/self.NormForces

            if self.nonlinear_iterative_technique == "newton_raphson":
                Eulerx, K, Residual = self.NewtonRaphson(function_spaces, formulation, solver,
                    Increment, K, NodalForces, Residual, mesh, Eulerx,
                    materials, boundary_condition, AppliedDirichletInc)
            elif self.nonlinear_iterative_technique == "modified_newton_raphson":
                Eulerx, K, Residual = self.ModifiedNewtonRaphson(function_spaces, formulation, solver,
                    Increment, K, NodalForces, Residual, mesh, Eulerx,
                    materials, boundary_condition, AppliedDirichletInc)
            elif self.nonlinear_iterative_technique == "line_search_newton_raphson":
                Eulerx, K, Residual = self.NewtonRaphsonLineSearch(function_spaces, formulation, solver,
                    Increment, K, NodalForces, Residual, mesh, Eulerx,
                    materials, boundary_condition, AppliedDirichletInc)
            else:
                raise RuntimeError("Iterative technique for nonlinear solver not understood")

            # UPDATE DISPLACEMENTS FOR THE CURRENT LOAD INCREMENT
            TotalDisp[:,:formulation.ndim,Increment] = Eulerx - mesh.points

            # PRINT LOG IF ASKED FOR
            self.LogSave(formulation, TotalDisp, Increment)

            print('\nFinished Load increment', Increment, 'in', time()-t_increment, 'seconds')
            try:
                print('Norm of Residual is',
                    np.abs(la.norm(Residual[boundary_condition.columns_in])/self.NormForces), '\n')
            except RuntimeWarning:
                print("Invalid value encountered in norm of Newton-Raphson residual")

            # STORE THE INFORMATION IF NEWTON-RAPHSON FAILS
            if self.newton_raphson_failed_to_converge:
                solver.condA = np.NAN
                Increment = Increment if Increment!=0 else 1
                TotalDisp = TotalDisp[:,:,:Increment]
                self.number_of_load_increments = Increment
                break

            # BREAK AT A SPECIFICED LOAD INCREMENT IF ASKED FOR
            if self.break_at_increment != -1 and self.break_at_increment is not None:
                if self.break_at_increment == Increment:
                    if self.break_at_increment < LoadIncrement - 1:
                        print("\nStopping at increment {} as specified\n\n".format(Increment))
                        TotalDisp = TotalDisp[:,:,:Increment]
                        self.number_of_load_increments = Increment
                    break

            # UPDATE. MATERIAL ADAPTATIVE LOAD FACTOR. FOR DEPOSITION STRETCH
            for imat in range(len(materials)):
                if materials[imat].load_factor is not None and Increment is not (LoadIncrement-1):
                    materials[imat].factor_increment += materials[imat].load_factor[Increment+1]


        return TotalDisp

    def NewtonRaphson(self, function_spaces, formulation, solver,
        Increment, K, NodalForces, Residual, mesh, Eulerx, materials,
        boundary_condition, AppliedDirichletInc):

        Tolerance = self.newton_raphson_tolerance
        LoadIncrement = self.number_of_load_increments
        Iter = 0
        self.iterative_norm_history = []

        # APPLY INCREMENTAL DIRICHLET PER LOAD STEP (THIS IS INCREMENTAL NOT ACCUMULATIVE)
        IncDirichlet = boundary_condition.UpdateFixDoFs(AppliedDirichletInc,
            K.shape[0],formulation.nvar)
        # UPDATE EULERIAN COORDINATE
        Eulerx += IncDirichlet[:,:formulation.ndim]

        while self.norm_residual > Tolerance or Iter==0:
            # GET THE REDUCED SYSTEM OF EQUATIONS
            K_b, F_b = boundary_condition.GetReducedMatrices(K,Residual)[:2]

            # SOLVE THE SYSTEM
            sol = solver.Solve(K_b,-F_b)

            # GET ITERATIVE SOLUTION
            dU = boundary_condition.UpdateFreeDoFs(sol,K.shape[0],formulation.nvar)

            # UPDATE THE EULERIAN COMPONENTS
            # UPDATE THE GEOMETRY
            Eulerx += dU[:,:formulation.ndim]

            # RE-ASSEMBLE - COMPUTE STIFFNESS AND INTERNAL TRACTION FORCES
            K, TractionForces = Assemble(self, function_spaces, formulation, mesh, materials,
                boundary_condition, Eulerx)[:2]
            # COMPUTE ROBIN STIFFNESS AND FORCES (EXTERNAL)
            K, TractionForces = boundary_condition.ComputeRobinForces(mesh, materials, function_spaces,
                self, Eulerx, K, TractionForces)

            # FIND THE RESIDUAL
            Residual[boundary_condition.columns_in] = TractionForces[boundary_condition.columns_in] - \
                NodalForces[boundary_condition.columns_in]

            # SAVE THE NORM
            self.abs_norm_residual = la.norm(Residual[boundary_condition.columns_in])
            if Iter==0:
                self.NormForces = la.norm(Residual[boundary_condition.columns_in])
            self.norm_residual = np.abs(la.norm(Residual[boundary_condition.columns_in])/self.NormForces)

            # SAVE THE NORM
            self.NRConvergence['Increment_'+str(Increment)] = np.append(self.NRConvergence['Increment_'+str(Increment)],\
                self.norm_residual)

            print("Iteration {} for increment {}.".format(Iter, Increment) +\
                " Residual (abs) {0:>16.7g}".format(self.abs_norm_residual),
                "\t Residual (rel) {0:>16.7g}".format(self.norm_residual))

            # BREAK BASED ON RELATIVE NORM
            if np.abs(self.abs_norm_residual) < Tolerance:
                break

            # BREAK BASED ON INCREMENTAL SOLUTION - KEEP IT AFTER UPDATE
            if norm(dU) <=  self.newton_raphson_solution_tolerance: # and Iter!=0:
                print("Incremental solution within tolerance i.e. norm(dU): {}".format(norm(dU)))
                break

            # UPDATE ITERATION NUMBER
            Iter +=1

            if Iter==self.maximum_iteration_for_newton_raphson:
                self.newton_raphson_failed_to_converge = True
                break
            if np.isnan(self.norm_residual) or self.norm_residual>1e06:
                self.newton_raphson_failed_to_converge = True
                break

            # IF BREAK WHEN NEWTON RAPHSON STAGNATES IS ACTIVATED
            if self.break_at_stagnation:
                self.iterative_norm_history.append(self.norm_residual)
                if Iter >= 5:
                    if np.mean(self.iterative_norm_history) < 1.:
                        break

            # USER DEFINED CRITERIA TO BREAK OUT OF NEWTON-RAPHSON
            if self.user_defined_break_func != None:
                if self.user_defined_break_func(Increment,Iter,self.norm_residual,self.abs_norm_residual, Tolerance):
                    break

            # USER DEFINED CRITERIA TO STOP NEWTON-RAPHSON AND THE WHOLE ANALYSIS
            if self.user_defined_stop_func != None:
                if self.user_defined_stop_func(Increment,Iter,self.norm_residual,self.abs_norm_residual, Tolerance):
                    self.newton_raphson_failed_to_converge = True
                    break

        return Eulerx, K, Residual

    def LogSave(self, formulation, TotalDisp, Increment):

            # PRINT LOG IF ASKED FOR
            if self.print_incremental_log:
                dmesh = Mesh()
                dmesh.points = TotalDisp[:,:formulation.ndim,Increment]
                dmesh_bounds = dmesh.Bounds
                print("\nMinimum and maximum incremental solution values at increment {} are \n".format(Increment),dmesh_bounds)

            # SAVE INCREMENTAL SOLUTION IF ASKED FOR
            if self.save_incremental_solution:
                # FOR BIG MESHES
                if Increment % self.incremental_solution_save_frequency !=0:
                    return
                from scipy.io import savemat
                filename = self.incremental_solution_filename
                if filename is not None:
                    if ".mat" in filename:
                        filename = filename.split(".")[0]
                    savemat(filename+"_"+str(Increment),
                        {'solution':TotalDisp[:,:,Increment]},do_compression=True)
                else:
                    raise ValueError("No file name provided to save incremental solution")

    def PrintPreAnalysisInfo(self, mesh, formulation):

        print('Pre-processing the information. Getting paths, solution parameters, mesh info, interpolation info etc...')
        print('Number of nodes is',mesh.points.shape[0], 'and number of DoFs is', mesh.points.shape[0]*formulation.nvar)
        if formulation.ndim==2:
            print('Number of elements is', mesh.elements.shape[0], \
                 'and number of boundary nodes is', np.unique(mesh.edges).shape[0], \
                 ' and number of element sets is', mesh.nset_elem)
        elif formulation.ndim==3:
            print('Number of elements is', mesh.elements.shape[0], \
                 'and number of boundary nodes is', np.unique(mesh.faces).shape[0], \
                 ' and number of element sets is', np.unique(mesh.element_to_set).shape[0])

    def ComputeSparsityFEM(self, mesh, formulation):

        if self.is_sparsity_pattern_computed is False:
            if self.recompute_sparsity_pattern is False and self.mass_type != "lumped":
                if self.parallel:
                    raise ValueError("Parallel model cannot use precomputed sparsity pattern due to partitioning. Turn this off")
                t_sp = time()
                from Kuru.FiniteElements.Assembly.ComputeSparsityPattern import ComputeSparsityPattern
                if self.squeeze_sparsity_pattern:
                    self.indices, self.indptr = ComputeSparsityPattern(mesh, formulation.nvar, self.squeeze_sparsity_pattern)
                    mesh.element_sorter = np.argsort(mesh.elements,axis=1)
                    mesh.sorted_elements = mesh.elements[np.arange(mesh.nelem)[:,None], mesh.element_sorter]
                    mesh.sorted_elements = mesh.sorted_elements.astype(np.uint64)

                    self.data_local_indices = self.data_global_indices = np.zeros(1,np.int32)
                else:
                    self.indices, self.indptr, self.data_local_indices, \
                        self.data_global_indices = ComputeSparsityPattern(mesh, formulation.nvar)

                self.is_sparsity_pattern_computed = True
                print("Computed sparsity pattern for the mesh. Time elapsed is {} seconds".format(time()-t_sp))

            else:
                self.indices, self.indptr, self.data_local_indices,\
                self.data_global_indices = np.array([0],dtype=np.int32), np.array([0],dtype=np.int32),\
                np.array([0],dtype=np.int32), np.array([0],dtype=np.int32)

    def CloseDaskDistributedClient(self):
        if self.parallel and self.parallel_model == "dask" and self.is_dask_scheduler_initialised:
            self.dask_client.close()
