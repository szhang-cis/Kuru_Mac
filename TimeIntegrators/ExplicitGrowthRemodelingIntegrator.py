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

from Kuru.FiniteElements.Assembly import Assemble #, AssembleExplicit
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
        mesh, TotalDisp, Eulerx, Eulerx0, materials, boundary_condition, fem_solver):

        if len(materials)>1 and self.density_turnover=="muscle":
            warn("More than one material. I will assume the first material is the Media")

        # check materials with growth and remodeling
        gr_materials = []
        for imat in range(len(materials)):
            if materials[imat].has_growth_remodeling:
                gr_materials.append(imat)
        gr_materials = np.array(gr_materials, dtype=np.int64, copy=True).flatten()

        GRVariables = [[] for i in range(gr_materials.shape[0])]
        for imat in range(gr_materials.shape[0]):
            GRVariables[imat] = np.zeros((materials[gr_materials[imat]].node_set.shape[0],12,fem_solver.number_of_time_increments),
                    dtype=np.float64)

        TimeIncrements = fem_solver.number_of_time_increments
        TimeFactor = fem_solver.total_time/TimeIncrements

        LoadIncrements = fem_solver.number_of_load_increments
        LoadFactor = 1./LoadIncrements
        AppliedDirichletInc = np.zeros(boundary_condition.applied_dirichlet.shape[0],dtype=np.float64)
        #print(boundary_condition.applied_dirichlet.shape[0])
        #exit()
        # HOMEOSTATIC STATE
        #TotalDisp = self.GetHomeostaticState(function_spaces, formulation, solver, 
        #        K, NeumannForces, NodalForces, Residual, mesh, TotalDisp, Eulerx, 
        #        material, boundary_condition, fem_solver, AppliedDirichletInc)

        IncrementalTime = 0.0
        # TIME LOOP
        for TIncrement in range(TimeIncrements):
            #temprary test of gain parameter after stenting
            #if (TIncrement >= 100):
            #    self.gain=0.15

            #pressure flags for this increment (always true for all inner surfaces)
            pressure_flags_inc = np.atleast_2d(boundary_condition.pressure_flags[:, TIncrement]).T
            # Block for stent control
            if (TIncrement >= 100000000): #number of increment to apply spring-like stent forces
                for face in np.where(pressure_flags_inc == True)[0]:
                    coord = Eulerx[mesh.faces[face, :], :]
                    avg = np.mean(coord, axis=0)
                    if (avg[2] <= 75):
                        boundary_condition.spring_flags[face,:] = True
                        boundary_condition.applied_spring[face,:] = 0.02
                    else:
                        boundary_condition.spring_flags[face,:] = False
                        boundary_condition.applied_spring[face,:] = 0.0
            # Block for pressure control
            #??? relocated to assembly
            #
            #update the dirichlet boundary conditions in stent contact zones
            #t1 = boundary_condition.columns_out
            #t2 = boundary_condition.applied_dirichlet
            #t3 = boundary_condition.columns_in
            #
            #for face in range(mesh.faces.shape[0]):
                #if boundary_condition.neumann_flags[face][0] == True:
                    #for node in mesh.faces[face,:]:
                        #coord = Eulerx[node, :]
                        #r = np.linalg.norm(coord[0:2])
                        #r0 = 7.88  # simulation results defined
                        #os = 0.2  # user defined
                        #r_free = r0 * (1 + os)  # 11 = 10 * (1+ 10%)  10% oversizing
                        #if (r < r_free) and (TIncrement>=110):
                        #if (coord[2]<75) and (coord[2]>55):
                            #if (TIncrement>=120000): #actuellement pas active
                                #ind = 3*node-1
                                #if(~(np.any(t1[:] == ind))):
                                    #t1=np.insert(t1,t1.shape[0],ind)
                                    #t2=np.insert(t2,t2.shape[0],0)
            #t3 = np.delete(np.arange(0, 3 * mesh.points.shape[0]), t1)
            #
            #boundary_condition.columns_out = t1
            #boundary_condition.applied_dirichlet=t2
            #boundary_condition.columns_in=t3

            #
            # CHECK ADAPTIVE STEP TIME FACTOR
            if fem_solver.time_factor is not None:
                TimeFactor = fem_solver.time_factor[TIncrement]
            Delta_t = TimeFactor

            # COMPUTE THE GROWTH AND REMODELING
            if TIncrement != 0:
                if self.has_dissection:
                    self.ComputeDissectionSpread(boundary_condition, mesh, IncrementalTime, Eulerx, TIncrement)
                for imat in range(gr_materials.shape[0]):
                    Rates = self.RatesGrowthRemodeling(mesh, materials[gr_materials[imat]], FibreStress, Softness, imat)
                    GRVariables[imat][:,:,TIncrement] = self.ExplicitGrowthRemodeling(mesh, materials[gr_materials[imat]],
                        IncrementalTime, Delta_t, Rates)
            else:
                for imat in range(gr_materials.shape[0]):
                    materials[gr_materials[imat]].den0_e = np.copy(materials[gr_materials[imat]].state_variables[:,14])
                    GRVariables[imat][:,:,0] = np.copy(materials[gr_materials[imat]].state_variables[:,9:21])
            IncrementalLoad = 1.0

            print("=======================================================")
            print("== Time elapsed, {:4.3f} days".format(IncrementalTime))
            print("== Current increment, {:4.3f}".format(TIncrement))
            print("== Current time increment, {:4.3f} days".format(Delta_t))
            print("=======================================================")

            # MATERIALS CHANGE AT EACH STEP
            if TIncrement != 0:
                tAssembly = time()
                K, TractionForces = Assemble(fem_solver, function_spaces, formulation, mesh, materials,
                    boundary_condition, Eulerx)[:2]

            # APPLY ROBIN BOUNDARY CONDITIONS - STIFFNESS(_) AND FORCES
            boundary_condition.pressure_increment = IncrementalLoad
            K, RobinForces = boundary_condition.ComputeRobinForces(mesh, materials, function_spaces,
                fem_solver, Eulerx, Eulerx0, K, np.zeros_like(Residual),TIncrement)

            #if (TIncrement >=20):
            #    boundary_condition.neumann_flags[:] = False
            # APPLY NEUMANN BOUNDARY CONDITIONS
            #if (TIncrement == 100):
                #NeumannForces = boundary_condition.ComputeNeumannForces(mesh, materials, function_spaces, Eulerx, TIncrement, compute_traction_forces=True, compute_body_forces=False)
            #
            #if (TIncrement >= 110):
                #NeumannForces = boundary_condition.ComputeNeumannForces(mesh, materials, function_spaces, Eulerx, TIncrement, compute_traction_forces=True, compute_body_forces=False)
            #
            #DeltaF = LoadFactor*np.array([NeumannForces]).T
            DeltaF = LoadFactor * NeumannForces
            #
            NodalForces += DeltaF
            #print(np.linalg.norm(NodalForces))
            #
            # OBRTAIN INCREMENTAL RESIDUAL - CONTRIBUTION FROM BOTH NEUMANN AND DIRICHLET
            Residual = -boundary_condition.ApplyDirichletGetReducedMatrices(K,np.zeros_like(Residual),
                boundary_condition.applied_dirichlet,LoadFactor=LoadFactor,only_residual=True)

            Residual += RobinForces - DeltaF
            # GET THE INCREMENTAL DISPLACEMENT
            AppliedDirichletInc = LoadFactor*boundary_condition.applied_dirichlet

            if TIncrement != 0:
                print('Finished the assembly stage. Time elapsed was', time()-tAssembly, 'seconds')
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
                        TIncrement, K, NodalForces, Residual, mesh, Eulerx, Eulerx0, materials,
                        boundary_condition, AppliedDirichletInc, fem_solver)
            else:
                raise RuntimeError("Iterative technique for nonlinear solver not understood")

            # UPDATE DISPLACEMENTS FOR THE CURRENT LOAD INCREMENT
            TotalDisp[:,:formulation.ndim,TIncrement] = Eulerx - mesh.points
            #
            #Eulerx0 = Eulerx
            #CHECK HOMEOSTATIC DISTORTION
            #if TIncrement==0:
            #    self.HomeostaticDistortion(fem_solver, formulation, TotalDisp, TIncrement)

            # COMPUTE THE FIBRE-STRESS AND SOFTNESS
            FibreStress = [[] for i in range(gr_materials.shape[0])]
            Softness = [[] for i in range(gr_materials.shape[0])]
            for imat in range(gr_materials.shape[0]):
                FibreStress[imat], Softness[imat] = self.GetFibreStressAndSoftness(mesh, formulation, 
                        materials[gr_materials[imat]], fem_solver, Eulerx)
            # SET HOMEOSTATIC FIBRE-STRESS
            if TIncrement==0:
                self.HomeostaticStress = [[] for i in range(gr_materials.shape[0])]
                for imat in range(gr_materials.shape[0]):
                    self.HomeostaticStress[imat] = FibreStress[imat]

            # PRINT LOG IF ASKED FOR
            self.LogSave(fem_solver, formulation, TotalDisp, TIncrement, materials, FibreStress, gr_materials)

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
                for imat in range(gr_materials.shape[0]):
                    GRVariables[imat] = GRVariables[imat][:,:,:TIncrement]
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

#            # UPDATE. MATERIAL ADAPTATIVE LOAD FACTOR. FOR DEPOSITION STRETCH
#            if Increment is not (LoadIncrement-1):
#                for imat in range(len(materials)):
#                    if materials[imat].load_factor is not None:
#                        materials[imat].factor_increment += materials[imat].load_factor[Increment+1]

        return TotalDisp, GRVariables

    def ExplicitGrowthRemodeling(self, mesh, material, IncrementalTime, Delta_t, Rates):
        """ Routine to get the evolution of Growth and Remodeling parameters
        """

        # Elastin degradation
        den0_tot = material.rho
        D_max = 0.5 #0.5
        L_dam = self.damage_spread_space
        t_dam = self.damage_spread_time
        T_ela = 101.0*365.25 #101.0*365.25

        # Loop on nodes
        for node in range(material.node_set.shape[0]):
            # Elastin function density in time f(t), analytic solution
            # degradation at line
            if self.degradation_at_line:
                AxialCoord = mesh.points[material.node_set[node],self.damage_axis]
                material.state_variables[node,14] = material.den0_e[node]*np.exp(-IncrementalTime/T_ela) + \
                    material.den0_e[node]*(D_max/t_dam)*(T_ela*t_dam/(t_dam-T_ela))*np.exp(-0.5*(AxialCoord/L_dam)**2)*\
                    (np.exp(-IncrementalTime/T_ela)-np.exp(-IncrementalTime/t_dam))
            # degradation at point
            elif self.degradation_at_point:
                vec_dam = mesh.points[material.node_set[node],:] - self.degradation_point
                R_dam = np.linalg.norm(vec_dam)
                material.state_variables[node,14] = material.den0_e[node]*np.exp(-IncrementalTime/T_ela) + \
                    material.den0_e[node]*(D_max/t_dam)*(T_ela*t_dam/(t_dam-T_ela))*np.exp(-0.5*(R_dam/L_dam)**2)*\
                    (np.exp(-IncrementalTime/T_ela)-np.exp(-IncrementalTime/t_dam))
            # just natural rate of elastin degradation
            elif self.aging_only:
                material.state_variables[node,14] = material.den0_e[node]*np.exp(-IncrementalTime/T_ela)
            # there is not degradation of elastin
            elif self.not_degradation:
                material.state_variables[node,14] = material.den0_e[node]
            else:
                raise ValueError("Degradation type of elastin not undesrtood. Insert eather at point or at line.")
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

    def ComputeDissectionSpread(self, boundary_condition, mesh, IncrementalTime, Eulerx, TIncrement):
        """ Routine to compute the spread of the dissection"""

        # a ramp function for the dissection spread wich stop after a time threshold
        #dissection_speed = self.maximum_dissection_spread/self.dissection_time_threshold
        #if IncrementalTime < self.dissection_time_threshold:
        #    dissection_spread = dissection_speed*IncrementalTime
        #else:
        #    dissection_spread = self.maximum_dissection_spread
        dissection_spread = self.maximum_dissection_spread
        #text_file = open(str(TIncrement)+"Output.txt","w")

        divider = boundary_condition.connector_elements.shape[1]
        # Loop into the connector elements
        for elem in range(boundary_condition.connector_elements.shape[0]):
            center = np.sum(mesh.points[boundary_condition.connector_elements[elem,:],:],axis=0)/divider
            #print(mesh.points[boundary_condition.connector_elements[elem,:],:])
            masterface = boundary_condition.connector_faces[elem][0]
            slaveface = boundary_condition.connector_faces[elem][1]
            masterAvgCoord = np.mean(Eulerx[mesh.faces[masterface, :], :],axis=0)
            slaveAvgCoord =  np.mean(Eulerx[mesh.faces[slaveface, :], :],axis=0)
            lspring = np.linalg.norm(masterAvgCoord-slaveAvgCoord)
            #text_file.write("%s,%s\n" %(elem,lspring))
            #exit() lbreak depending on z position
            lbreak = 0.0125 #should be defined according to
            #if center[2] < dissection_spread:
            if center[2] < dissection_spread or lspring > lbreak:
                boundary_condition.connector_flags[elem] = False
                free_faces = boundary_condition.connector_faces[elem]
                boundary_condition.pressure_flags[free_faces] = True
                boundary_condition.applied_pressure[free_faces] = self.dissection_pressure
        #text_file.close
        return