from __future__ import print_function
import platform #,os
from time import time
from warnings import warn
import numpy as np
#import scipy as sp
from scipy.sparse import issparse, isspmatrix_csr #, isspmatrix_coo, isspmatrix_csc
from scipy.sparse.linalg import spsolve, cg, gmres, lgmres, bicgstab#, cgs, splu, spilu, LinearOperator, onenormest
#from subprocess import call

__all__ = ["LinearSolver"]

class LinearSolver(object):
    """Base class for all linear sparse direct and iterative solvers"""

    def __init__(self,linear_solver="direct", linear_solver_type="umfpack",
        apply_preconditioner=False, preconditioner="smoothed_aggregation",
        iterative_solver_tolerance=1.0e-12, reduce_matrix_bandwidth=False,
        out_of_core=False, geometric_discretisation=None, dont_switch_solver=False):
        """
            input:
                linear_solver:          [str] type of solver either "direct",
                                        "iterative", "petsc" or "amg"
                linear_solver_type      [str] type of direct or linear solver to
                                        use, for instance "umfpack", "superlu" or
                                        "mumps" for direct solvers, or "cg", "gmres"
                                        etc for iterative solvers or "amg" for algebraic
                                        multigrid solver. See WhichSolvers method for
                                        the complete set of available linear solvers
                preconditioner:         [str] either "smoothed_aggregation",
                                        or "ruge_stuben" or "rootnode" for
                                        a preconditioner based on algebraic multigrid
                                        or "ilu" for scipy's spilu linear
                                        operator
                geometric_discretisation:
                                        [str] type of geometric discretisation used, for
                                        instance for FEM discretisations this would correspond
                                        to "tri", "quad", "tet", "hex" etc
                dont_switch_solver:     Do not switch between solvers automatically
        """

        self.is_sparse = True
        self.solver_type = linear_solver
        self.solver_subtype = linear_solver_type
        self.requires_cuthill_mckee = reduce_matrix_bandwidth
        self.iterative_solver_tolerance = iterative_solver_tolerance
        self.apply_preconditioner = apply_preconditioner
        self.preconditioner_type = preconditioner
        self.out_of_core = False
        self.geometric_discretisation = geometric_discretisation
        self.dont_switch_solver = dont_switch_solver
        self.reuse_factorisation = False
        self.solver_context_manager = None

        self.has_amg_solver = True
        if platform.python_implementation() == "PyPy":
            self.has_amg_solver = False
        else:
            try:
                import pyamg
            except ImportError:
                self.has_amg_solver = False

        self.has_umfpack = True
        try:
            from scikits.umfpack import spsolve
        except (ImportError, AttributeError) as umfpack_error:
            self.has_umfpack = False

        self.has_mumps = False
        try:
            from _dMUMPS_ import dMUMPS_solve
            self.has_mumps = True
        except ImportError:
            self.has_mumps = False

        self.has_pardiso = False
        try:
            import pypardiso
            self.has_pardiso = True
        except ImportError:
            self.has_pardiso = False

        self.has_petsc = False
        try:
            import petsc4py
            self.has_petsc = True
        except ImportError:
            self.has_petsc = False

        self.switcher_message = False

    def Solve(self, A, b, reuse_factorisation=False):
        """Solves the linear system of equations"""

        if not issparse(A):
            raise ValueError("Linear system is not of sparse type")

        if A.shape == (0,0) and b.shape[0] == 0:
            warn("Empty linear system!!! Nothing to solve!!!")
            return np.copy(b)


        self.reuse_factorisation = reuse_factorisation
        if self.solver_type != "direct" and self.reuse_factorisation is True:
            warn("Re-using factorisation for non-direct solvers is not possible. The pre-conditioner is going to be reused instead")

        # DECIDE IF THE SOLVER TYPE IS APPROPRIATE FOR THE PROBLEM
        if self.switcher_message is False and self.dont_switch_solver is False:
            # PREFER PARDISO OR MUMPS OVER AMG IF AVAILABLE
            if self.has_pardiso:
                self.solver_type = "direct"
                self.solver_subtype = "pardiso"
            elif self.has_mumps:
                self.solver_type = "direct"
                self.solver_subtype = "mumps"
            elif b.shape[0] > 100000 and self.has_amg_solver:
                self.solver_type = "amg"
                self.solver_subtype = "gmres"
                print('Large system of equations. Switching to algebraic multigrid solver')
                self.switcher_message = True
            # elif mesh.points.shape[0]*MainData.nvar > 50000 and MainData.C < 4:
                # self.solver_type = "direct"
                # self.solver_subtype = "MUMPS"
                # print 'Large system of equations. Switching to MUMPS solver'
            elif b.shape[0] > 70000 and self.geometric_discretisation=="hex" and self.has_amg_solver:
                self.solver_type = "amg"
                self.solver_subtype = "gmres"
                print('Large system of equations. Switching to algebraic multigrid solver')
                self.switcher_message = True
            else:
                self.solver_type = "direct"
                self.solver_subtype = "umfpack"


        if self.solver_type == 'direct':
            # CALL DIRECT SOLVER
            if self.solver_subtype=='umfpack' and self.has_umfpack:
                if A.dtype != np.float64:
                    A = A.astype(np.float64)

                t_solve = time()
                if self.solver_context_manager is None:
                    if self.reuse_factorisation is False:
                        sol = spsolve(A,b,permc_spec='MMD_AT_PLUS_A',use_umfpack=True)
                        # from scikits import umfpack
                        # sol = umfpack.spsolve(A, b)

                        # SuperLU
                        # from scipy.sparse.linalg import splu
                        # lu = splu(A.tocsc())
                        # sol = lu.solve(b)
                    else:
                        from scikits import umfpack
                        lu = umfpack.splu(A)
                        sol = lu.solve(b)
                        self.solver_context_manager = lu
                else:
                    sol = self.solver_context_manager.solve(b)

                # print("UmfPack solver time is {}".format(time() - t_solve))


            elif self.solver_subtype=='mumps' and self.has_mumps:

                from _dMUMPS_ import dMUMPS_solve
                t_solve = time()
                A = A.tocoo()
                sol=dMUMPS_solve(A.data,A.row,A.col,b)
                # False means non-symmetric - Do not change it to True. True means symmetric pos def
                # which is not the case for electromechanics
                #if self.solver_context_manager is None:
                    #context = MUMPSContext((A.shape[0], A.row, A.col, A.data, False), verbose=False)
                    #context.analyze()
                    #context.factorize()
                    #sol = context.solve(rhs=b)

                    #if self.reuse_factorisation:
                    #    self.solver_context_manager = context
                #else:
                    #sol = self.solver_context_manager.solve(rhs=b)

                print("MUMPS solver time is {}".format(time() - t_solve))

                return sol


            elif self.solver_subtype == "pardiso" and self.has_pardiso:
                # NOTE THAT THIS PARDISO SOLVER AUTOMATICALLY SAVES THE RIGHT FACTORISATION
                import pypardiso
                from pypardiso.scipy_aliases import pypardiso_solver as ps
                A = A.tocsr()
                t_solve = time()
                sol = pypardiso.spsolve(A,b)
                if self.reuse_factorisation is False:
                    ps.remove_stored_factorization()
                    ps.free_memory()
                print("Pardiso solver time is {}".format(time() - t_solve))


            else:
                # FOR 'super_lu'
                if A.dtype != np.float64:
                    A = A.astype(np.float64)
                A = A.tocsc()

                if self.solver_context_manager is None:
                    if self.reuse_factorisation is False:
                        sol = spsolve(A,b,permc_spec='MMD_AT_PLUS_A',use_umfpack=True)
                    else:
                        lu = splu(A)
                        sol = lu.solve(b)
                        self.solver_context_manager = lu
                else:
                    sol = self.solver_context_manager.solve(b)



        elif self.solver_type == "iterative":
            t_solve = time()
            # CALL ITERATIVE SOLVER
            if self.solver_subtype == "gmres":
                sol = gmres(A,b,tol=self.iterative_solver_tolerance)[0]
            if self.solver_subtype == "lgmres":
                sol = lgmres(A,b,tol=self.iterative_solver_tolerance)[0]
            elif self.solver_subtype == "bicgstab":
                sol = bicgstab(A,b,tol=self.iterative_solver_tolerance)[0]
            else:
                sol = cg(A,b,tol=self.iterative_solver_tolerance)[0]

            # PRECONDITIONED ITERATIVE SOLVER - CHECK
            # P = spilu(A.tocsc(), drop_tol=1e-5)
            # M_x = lambda x: P.solve(x)
            # m = A.shape[1]
            # n = A.shape[0]
            # M = LinearOperator((n * m, n * m), M_x)
            # sol = cg(A, b, tol=self.iterative_solver_tolerance, M=M)[0]
            print("Iterative solver time is {}".format(time() - t_solve))

        elif self.solver_type == "amg":
            if self.has_amg_solver is False:
                raise ImportError('Algebraic multigrid solver was not found. Please install it using "pip install pyamg"')
            from pyamg import ruge_stuben_solver, rootnode_solver, smoothed_aggregation_solver

            if A.dtype != b.dtype:
                # DOWN-CAST
                b = b.astype(A.dtype)

            if not isspmatrix_csr(A):
                A = A.tocsr()

            t_solve = time()

            if self.iterative_solver_tolerance > 1e-9:
                self.iterative_solver_tolerance = 1e-10

            # AMG METHOD
            amg_func = None
            if self.preconditioner_type=="smoothed_aggregation":
                # THIS IS TYPICALLY FASTER BUT THE TOLERANCE NEED TO BE SMALLER, TYPICALLY 1e-10
                amg_func = smoothed_aggregation_solver
            elif self.preconditioner_type == "ruge_stuben":
                amg_func = ruge_stuben_solver
            elif self.preconditioner_type == "rootnode":
                amg_func = rootnode_solver
            else:
                amg_func = rootnode_solver

            ml = amg_func(A)
            # ml = amg_func(A, smooth=('energy', {'degree':2}), strength='evolution' )
            # ml = amg_func(A, max_levels=3, diagonal_dominance=True)
            # ml = amg_func(A, coarse_solver=spsolve)
            # ml = amg_func(A, coarse_solver='cholesky')

            if self.solver_context_manager is None:
                # M = ml.aspreconditioner(cycle='V')
                M = ml.aspreconditioner()
                if self.reuse_factorisation:
                    self.solver_context_manager = M
            else:
                M = self.solver_context_manager

            # EXPLICIT CALL TO KYROLOV SOLVERS WITH AMG PRECONDITIONER
            # sol, info = bicgstab(A, b, M=M, tol=self.iterative_solver_tolerance)
            # sol, info = cg(A, b, M=M, tol=self.iterative_solver_tolerance)
            # sol, info = gmres(A, b, M=M, tol=self.iterative_solver_tolerance)

            # IMPLICIT CALL TO KYROLOV SOLVERS WITH AMG PRECONDITIONER
            residuals = []
            sol = ml.solve(b, tol=self.iterative_solver_tolerance, accel=self.solver_subtype, residuals=residuals)

            print("AMG solver time is {}".format(time() - t_solve))



        elif self.solver_type == "petsc" and self.has_petsc:
            if self.solver_subtype != "gmres" and self.solver_subtype != "minres" and self.solver_subtype != "cg":
                self.solver_subtype == "cg"
            if self.iterative_solver_tolerance < 1e-9:
                self.iterative_solver_tolerance = 1e-7

            from petsc4py import PETSc
            t_solve = time()
            pA = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
            pb = PETSc.Vec().createWithArray(b)

            ksp = PETSc.KSP()
            ksp.create(PETSc.COMM_WORLD)
            # ksp.create()
            ksp.setType(self.solver_subtype)
            ksp.setTolerances(atol=self.iterative_solver_tolerance,
                    rtol=self.iterative_solver_tolerance)
            # ILU
            ksp.getPC().setType('icc')

            # CREATE INITIAL GUESS
            psol = PETSc.Vec().createWithArray(np.ones(b.shape[0]))
            # SOLVE
            ksp.setOperators(pA)
            ksp.setFromOptions()
            ksp.solve(pb, psol)
            sol = psol.getArray()

            # print('Converged in', ksp.getIterationNumber(), 'iterations.')
            print("Petsc linear iterative solver time is {}".format(time() - t_solve))

        else:
            warn("{} solver is not available. Default solver is going to be used".format(self.solver_type))
            # FOR 'super_lu'
            if A.dtype != np.float64:
                A = A.astype(np.float64)
            A = A.tocsc()

            if self.solver_context_manager is None:
                if self.reuse_factorisation is False:
                    sol = spsolve(A,b,permc_spec='MMD_AT_PLUS_A',use_umfpack=True)
                else:
                    lu = splu(A)
                    sol = lu.solve(b)
                    self.solver_context_manager = lu
            else:
                sol = self.solver_context_manager.solve(b)


        return sol
