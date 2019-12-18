from __future__ import print_function
import gc #, os, sys
#from copy import deepcopy
from warnings import warn
from time import time
import numpy as np
from scipy.sparse import coo_matrix #, csc_matrix, csr_matrix

#from ._LowLevelAssembly_ import _LowLevelAssembly_, _LowLevelAssemblyExplicit_, _LowLevelAssemblyLaplacian_
#from ._LowLevelAssembly_ import _LowLevelAssembly_Par_, _LowLevelAssemblyExplicit_Par_

from .SparseAssemblyNative import SparseAssemblyNative #, SparseAssemblyNativeCSR, SparseAssemblyNativeCSR_RecomputeDataIndex
from .RHSAssemblyNative import RHSAssemblyNative
#from .ComputeSparsityPattern import ComputeSparsityPattern

# PARALLEL PROCESSING ROUTINES
#import multiprocessing
#import Florence.ParallelProcessing.parmap as parmap

__all__ = ['Assemble', 'AssembleForces', 'AssembleExplicit', 'AssembleMass', 'AssembleForm']

def Assemble(fem_solver, function_space, formulation, mesh, material, Eulerx, Eulerp):

    if fem_solver.has_low_level_dispatcher:
        return LowLevelAssembly(fem_solver, function_space, formulation, mesh, material, Eulerx, Eulerp)
    else:
        if mesh.nelem <= 600000:
            return AssemblySmall(fem_solver, function_space, formulation, mesh, material, Eulerx, Eulerp)
        elif mesh.nelem > 600000:
            print("Larger than memory system. Dask on disk parallel assembly is turned on")
            return OutofCoreAssembly(fem_solver, function_space, formulation, mesh, material, Eulerx, Eulerp)

def AssemblySmall(fem_solver, function_space, formulation, mesh, material, Eulerx, Eulerp):

    t_assembly = time()

    # GET MESH DETAILS
    C = mesh.InferPolynomialDegree() - 1
    nvar = formulation.nvar
    ndim = formulation.ndim
    nelem = mesh.nelem
    nodeperelem = mesh.elements.shape[1]
    ndof = nodeperelem*nvar
    local_capacity = ndof*ndof

    if fem_solver.recompute_sparsity_pattern is False:
        indices, indptr = fem_solver.indices, fem_solver.indptr
        if fem_solver.squeeze_sparsity_pattern is False:
            data_global_indices = fem_solver.data_global_indices
            data_local_indices = fem_solver.data_local_indices

    if fem_solver.recompute_sparsity_pattern:
        # ALLOCATE VECTORS FOR SPARSE ASSEMBLY OF STIFFNESS MATRIX - CHANGE TYPES TO INT64 FOR DoF > 1e09
        I_stiffness=np.zeros(int((nvar*nodeperelem)**2*nelem),dtype=np.int32)
        J_stiffness=np.zeros(int((nvar*nodeperelem)**2*nelem),dtype=np.int32)
        V_stiffness=np.zeros(int((nvar*nodeperelem)**2*nelem),dtype=np.float64)

        I_mass=[]; J_mass=[]; V_mass=[]
        if fem_solver.analysis_type !='static' and fem_solver.is_mass_computed is False:
            # ALLOCATE VECTORS FOR SPARSE ASSEMBLY OF MASS MATRIX - CHANGE TYPES TO INT64 FOR DoF > 1e09
            I_mass=np.zeros(int((nvar*nodeperelem)**2*nelem),dtype=np.int32)
            J_mass=np.zeros(int((nvar*nodeperelem)**2*nelem),dtype=np.int32)
            V_mass=np.zeros(int((nvar*nodeperelem)**2*nelem),dtype=np.float64)
    else:
        V_stiffness=np.zeros(indices.shape[0],dtype=np.float64)
        if fem_solver.analysis_type !='static' and fem_solver.is_mass_computed is False:
            V_mass=np.zeros(indices.shape[0],dtype=np.float64)


    T = np.zeros((mesh.points.shape[0]*nvar,1),np.float64)

    mass, F = [], []
    if fem_solver.has_moving_boundary:
        F = np.zeros((mesh.points.shape[0]*nvar,1),np.float64)


    if fem_solver.parallel:
        # COMPUATE ALL LOCAL ELEMENTAL MATRICES (STIFFNESS, MASS, INTERNAL & EXTERNAL TRACTION FORCES)
        ParallelTuple = parmap.map(formulation,np.arange(0,nelem,dtype=np.int32),
            function_space, mesh, material, fem_solver, Eulerx, Eulerp, processes= int(multiprocessing.cpu_count()/2))

    for elem in range(nelem):

        if fem_solver.parallel:
            # UNPACK PARALLEL TUPLE VALUES
            I_stiff_elem = ParallelTuple[elem][0]; J_stiff_elem = ParallelTuple[elem][1]; V_stiff_elem = ParallelTuple[elem][2]
            t = ParallelTuple[elem][3]; f = ParallelTuple[elem][4]
            I_mass_elem = ParallelTuple[elem][5]; J_mass_elem = ParallelTuple[elem][6]; V_mass_elem = ParallelTuple[elem][6]

        else:
            # COMPUATE ALL LOCAL ELEMENTAL MATRICES (STIFFNESS, MASS, INTERNAL & EXTERNAL TRACTION FORCES )
            I_stiff_elem, J_stiff_elem, V_stiff_elem, t, f, \
            I_mass_elem, J_mass_elem, V_mass_elem = formulation.GetElementalMatrices(elem,
                function_space, mesh, material, fem_solver, Eulerx, Eulerp)

        if fem_solver.recompute_sparsity_pattern:
            # SPARSE ASSEMBLY - STIFFNESS MATRIX
            SparseAssemblyNative(I_stiff_elem,J_stiff_elem,V_stiff_elem,I_stiffness,J_stiffness,V_stiffness,
                elem,nvar,nodeperelem,mesh.elements)

            if fem_solver.analysis_type != 'static' and fem_solver.is_mass_computed==False:
                # SPARSE ASSEMBLY - MASS MATRIX
                SparseAssemblyNative(I_mass_elem,J_mass_elem,V_mass_elem,I_mass,J_mass,V_mass,
                    elem,nvar,nodeperelem,mesh.elements)

        else:
            if fem_solver.squeeze_sparsity_pattern:
                # SPARSE ASSEMBLY - STIFFNESS MATRIX
                SparseAssemblyNativeCSR_RecomputeDataIndex(mesh,V_stiff_elem,indices,indptr,V_stiffness,elem,nvar)

                if fem_solver.analysis_type != 'static' and fem_solver.is_mass_computed==False:
                    # SPARSE ASSEMBLY - MASS MATRIX
                    SparseAssemblyNativeCSR_RecomputeDataIndex(mesh,V_mass_elem,indices,indptr,V_mass,elem,nvar)
            else:
                # SPARSE ASSEMBLY - STIFFNESS MATRIX
                V_stiffness[data_global_indices[elem*local_capacity:(elem+1)*local_capacity]] \
                    += V_stiff_elem[data_local_indices[elem*local_capacity:(elem+1)*local_capacity]]

                if fem_solver.analysis_type != 'static' and fem_solver.is_mass_computed==False:
                    # SPARSE ASSEMBLY - MASS MATRIX
                    V_mass[data_global_indices[elem*local_capacity:(elem+1)*local_capacity]] \
                    += V_mass_elem[data_local_indices[elem*local_capacity:(elem+1)*local_capacity]]

        if fem_solver.has_moving_boundary:
            # RHS ASSEMBLY
            # for iterator in range(0,nvar):
            #     F[mesh.elements[elem,:]*nvar+iterator,0]+=f[iterator::nvar,0]
            RHSAssemblyNative(F,f,elem,nvar,nodeperelem,mesh.elements)

        # INTERNAL TRACTION FORCE ASSEMBLY
        # for iterator in range(0,nvar):
            # T[mesh.elements[elem,:]*nvar+iterator,0]+=t[iterator::nvar,0]
        RHSAssemblyNative(T,t,elem,nvar,nodeperelem,mesh.elements)


    if fem_solver.parallel:
        del ParallelTuple
        gc.collect()

    if fem_solver.recompute_sparsity_pattern:
        stiffness = coo_matrix((V_stiffness,(I_stiffness,J_stiffness)),
            shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])),dtype=np.float64).tocsr()

        # GET STORAGE/MEMORY DETAILS
        fem_solver.spmat = stiffness.data.nbytes/1024./1024.
        fem_solver.ijv = (I_stiffness.nbytes + J_stiffness.nbytes + V_stiffness.nbytes)/1024./1024.

        del I_stiffness, J_stiffness, V_stiffness
        gc.collect()

        if fem_solver.analysis_type != 'static' and fem_solver.is_mass_computed==False:
            mass = csr_matrix((V_mass,(I_mass,J_mass)),shape=((nvar*mesh.points.shape[0],
                nvar*mesh.points.shape[0])),dtype=np.float64)

            fem_solver.is_mass_computed = True

    else:
        stiffness = csr_matrix((V_stiffness,indices,indptr),
            shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])))

        # GET STORAGE/MEMORY DETAILS
        fem_solver.spmat = stiffness.data.nbytes/1024./1024.
        fem_solver.ijv = (indptr.nbytes + indices.nbytes + V_stiffness.nbytes)/1024./1024.

        if fem_solver.analysis_type != 'static' and fem_solver.is_mass_computed==False:
            mass = csr_matrix((V_mass,indices,indptr),
                shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])))

            fem_solver.is_mass_computed = True

    fem_solver.assembly_time = time() - t_assembly

    return stiffness, T, F, mass

#------------------------------- ASSEMBLY ROUTINE FOR EXTERNAL TRACTION FORCES ----------------------------------#
#----------------------------------------------------------------------------------------------------------------#

def AssembleForces(boundary_condition, mesh, material, function_spaces, compute_traction_forces=True, compute_body_forces=False):

    Ft = np.zeros((mesh.points.shape[0]*material.nvar,1))
    Fb = np.zeros((mesh.points.shape[0]*material.nvar,1))

    if compute_traction_forces:
        Ft = AssembleExternalTractionForces(boundary_condition, mesh, material, function_spaces[-1])
    if compute_body_forces:
        Fb = AssembleBodyForces(boundary_condition, mesh, material, function_spaces[0])

    return Ft + Fb

def AssembleExternalTractionForces(boundary_condition, mesh, material, function_space):


    nvar = material.nvar
    ndim = material.ndim
    ngauss = function_space.AllGauss.shape[0]

    if ndim == 2:
        faces = mesh.edges
        nodeperelem = mesh.edges.shape[1]
    else:
        faces = mesh.faces
        nodeperelem = mesh.faces.shape[1]

    if boundary_condition.is_applied_neumann_shape_functions_computed is False:
        N = np.zeros((nodeperelem*nvar,nvar,ngauss))
        for i in range(nvar):
            N[i::nvar,i,:] = function_space.Bases
        boundary_condition.__Nt__ = N
        boundary_condition.is_applied_neumann_shape_functions_computed = True
    else:
        N = boundary_condition.__Nt__


    F = np.zeros((mesh.points.shape[0]*nvar,1))
    for face in range(faces.shape[0]):
        if boundary_condition.neumann_flags[face] == True:
            ElemTraction = boundary_condition.applied_neumann[face,:]
            external_traction = np.einsum("ijk,j,k->ik",N,ElemTraction,function_space.AllGauss[:,0]).sum(axis=1)
            RHSAssemblyNative(F,np.ascontiguousarray(external_traction[:,None]),face,nvar,nodeperelem,faces)


    return F
