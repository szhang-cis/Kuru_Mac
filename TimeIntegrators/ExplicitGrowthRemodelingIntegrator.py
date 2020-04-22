from __future__ import print_function
import gc, os, sys
import numpy as np
import scipy as sp
import numpy.linalg as la
import scipy.linalg as sla
from numpy.linalg import norm
from time import time
from copy import deepcopy
from warnings import warn
from time import time

#from Florence.FiniteElements.Assembly import Assemble, AssembleExplicit
#from Florence import Mesh
#from Florence.PostProcessing import PostProcess
from .GrowthRemodelingIntegrator import GrowthRemodelingIntegrator

__all__ = ["ExplicitGrowthRemodelingIntegrator"]


class ExplicitGrowthRemodelingIntegrator(GrowthRemodelingIntegrator):
    """Generic explicit structural time integerator based on forward Euler"""

    def __init__(self, gain, turnover, **kwargs):
        super(ExplicitGrowthRemodelingIntegrator, self).__init__(gain, turnover, **kwargs)

    def Solver(self, function_spaces, formulation, solver,
        K, NeumannForces, NodalForces, Residual,
        mesh, TotalDisp, Eulerx, material, boundary_condition, fem_solver):

        GRVariables = np.zeros((mesh.points.shape[0],12,fem_solver.number_of_time_increments),
                dtype=np.float64)

        TimeIncrements = fem_solver.number_of_time_increments
        TimeFactor = fem_solver.total_time/TimeIncrements

        LoadIncrements = fem_solver.number_of_load_increments
        LoadFactor = 1./LoadIncrements

        AppliedDirichletInc = np.zeros(boundary_condition.applied_dirichlet.shape[0],dtype=np.float64)

        # HOMEOSTATIC STATE
        #TotalDisp = self.GetHomeostaticState(function_spaces, formulation, solver, 
        #        K, NeumannForces, NodalForces, Residual, mesh, TotalDisp, Eulerx, 
        #        material, boundary_condition, fem_solver, AppliedDirichletInc)
        
        IncrementalTime = 0.0
        # TIME LOOP
        for TIncrement in range(TimeIncrements):

            # CHECK ADAPTIVE LOAD FACTOR
            if fem_solver.time_factor is not None:
                TimeFactor = fem_solver.time_factor[TIncrement]
            Delta_t = TimeFactor

            # COMPUTE THE GROWTH AND REMODELING
            if TIncrement != 0:
                Rates = self.RatesGrowthRemodeling(mesh, material, FibreStress, Softness)
                GRVariables[:,:,TIncrement] = self.ExplicitGrowthRemodeling(mesh, material, 
                    IncrementalTime, Delta_t, Rates)
            else:
                GRVariables[:,:,0] = material.state_variables[:,9:21]
            IncrementalLoad = 1.0

            print("=============================")
            print("== Time elapsed, {:4.0f} days ==".format(IncrementalTime))
            print("=============================")

            boundary_condition.pressure_increment = IncrementalLoad
            # APPLY NEUMANN BOUNDARY CONDITIONS
            DeltaF = LoadFactor*NeumannForces
            NodalForces += DeltaF
            # OBRTAIN INCREMENTAL RESIDUAL - CONTRIBUTION FROM BOTH NEUMANN AND DIRICHLET
            Residual = -boundary_condition.ApplyDirichletGetReducedMatrices(K,np.zeros_like(Residual),
                boundary_condition.applied_dirichlet,LoadFactor=LoadFactor,only_residual=True)
            Residual -= DeltaF
            # GET THE INCREMENTAL DISPLACEMENT
            AppliedDirichletInc = LoadFactor*boundary_condition.applied_dirichlet

            t_increment = time()

            # LET NORM OF THE FIRST RESIDUAL BE THE NORM WITH RESPECT TO WHICH WE
            # HAVE TO CHECK THE CONVERGENCE OF NEWTON RAPHSON. TYPICALLY THIS IS
            # NORM OF NODAL FORCES
            if TIncrement==0:
                fem_solver.NormForces = np.linalg.norm(Residual)
                # AVOID DIVISION BY ZERO
                if np.isclose(fem_solver.NormForces,0.0):
                    fem_solver.NormForces = 1e-14
            fem_solver.norm_residual = np.linalg.norm(Residual)/fem_solver.NormForces

            if fem_solver.nonlinear_iterative_technique == "newton_raphson":
                Eulerx, K, Residual = self.NewtonRaphson(function_spaces, formulation, solver,
                        TIncrement, K, NodalForces, Residual, mesh, Eulerx, material,
                        boundary_condition, AppliedDirichletInc, fem_solver)
            else:
                raise RuntimeError("Iterative technique for nonlinear solver not understood")

            # UPDATE DISPLACEMENTS FOR THE CURRENT LOAD INCREMENT
            TotalDisp[:,:formulation.ndim,TIncrement] = Eulerx - mesh.points

            #CHECK HOMEOSTATIC DISTORTION
            if TIncrement==0:
                self.HomeostaticDistortion(fem_solver, formulation, TotalDisp, TIncrement)

            # COMPUTE THE FIBRE-STRESS AND SOFTNESS
            FibreStress, Softness = self.GetFibreStressAndSoftness(mesh, formulation, material,
                    fem_solver, Eulerx)
            # SET HOMEOSTATIC FIBRE-STRESS
            if TIncrement==0:
                self.HomeostaticStress = FibreStress

            # PRINT LOG IF ASKED FOR
            self.LogSave(fem_solver, formulation, TotalDisp, TIncrement, material, FibreStress)

            # UPDATE THE TIME
            IncrementalTime += TimeFactor

            print('\nFinished Time increment', TIncrement, 'in', time()-t_increment, 'seconds')
            try:
                print('Norm of Residual is',
                    np.abs(la.norm(Residual[boundary_condition.columns_in])/fem_solver.NormForces), '\n')
            except RuntimeWarning:
                print("Invalid value encountered in norm of Newton-Raphson residual")

            # STORE THE INFORMATION IF NEWTON-RAPHSON FAILS
            if fem_solver.newton_raphson_failed_to_converge:
                solver.condA = np.NAN
                TIncrement = TIncrement if TIncrement!=0 else 1
                TotalDisp = TotalDisp[:,:,:TIncrement]
                fem_solver.number_of_time_increments = TIncrement
                break

            # BREAK AT A SPECIFICED TIME INCREMENT IF ASKED FOR
            if fem_solver.break_at_increment != -1 and fem_solver.break_at_increment is not None:
                if fem_solver.break_at_increment == TIncrement:
                    if fem_solver.break_at_increment < TimeIncrements - 1:
                        print("\nStopping at time increment {} as specified\n\n".format(TIncrement))
                        TotalDisp = TotalDisp[:,:,:TIncrement]
                        fem_solver.number_of_time_increments = TIncrement
                    break


        return TotalDisp,GRVariables

    def ExplicitGrowthRemodeling(self, mesh, material, IncrementalTime, Delta_t, Rates):
        """ Routine to get the evolution of Growth and Remodeling parameters
        """

        # Elastin degradation
        den0_tot = 1050.0
        den0_e = 241.5
        D_max = 0.5
        L_dam = 0.010
        t_dam = 40.0
        T_ela = 101.0*365.25

        # Loop on nodes
        for node in range(mesh.nnode):
            # Elastin function density in time f(t), analytic solution
            AxialCoord = mesh.points[node,1]
            material.state_variables[node,14] = den0_e*np.exp(-IncrementalTime/T_ela) + \
                den0_e*(D_max/t_dam)*(T_ela*t_dam/(t_dam-T_ela))*np.exp(-0.5*(AxialCoord/L_dam)**2)*\
                (np.exp(-IncrementalTime/T_ela)-np.exp(-IncrementalTime/t_dam))
            # Time Integration of fibre densities
            for fibre in range(5):
                material.state_variables[node,15+fibre] = material.state_variables[node,15+fibre] + \
                        Delta_t*Rates[node,fibre+5]
                material.state_variables[node,9+fibre] = material.state_variables[node,9+fibre] + \
                        Delta_t*Rates[node,fibre]
            material.state_variables[node,20] = 0.0
            for n in range(6):
                material.state_variables[node,20] += material.state_variables[node,14+n]/den0_tot

        return material.state_variables[:,9:21]

