from __future__ import print_function
import gc #, os, sys
#from copy import deepcopy
from warnings import warn
from time import time
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix

#from ._LowLevelAssembly_ import _LowLevelAssembly_, _LowLevelAssemblyExplicit_, _LowLevelAssemblyLaplacian_
#from ._LowLevelAssembly_ import _LowLevelAssembly_Par_, _LowLevelAssemblyExplicit_Par_

from .SparseAssemblyNative import SparseAssemblyNative #, SparseAssemblyNativeCSR, SparseAssemblyNativeCSR_RecomputeDataIndex
from .RHSAssemblyNative import RHSAssemblyNative
#from .ComputeSparsityPattern import ComputeSparsityPattern

# PARALLEL PROCESSING ROUTINES
#import multiprocessing
#import Florence.ParallelProcessing.parmap as parmap

__all__ = ['Assemble', 'AssembleForces']

#----------------------------------------------------------------------------------------------------------------#
#------------------------------- ASSEMBLY ROUTINE FOR INTERNAL TRACTION FORCES ----------------------------------#
#----------------------------------------------------------------------------------------------------------------#
def Assemble(fem_solver, function_spaces, formulation, mesh, material, boundary_condition, Eulerx):

    if fem_solver.has_low_level_dispatcher:
        return LowLevelAssembly(fem_solver, function_spaces[0], formulation, mesh, material, Eulerx)
    else:
        if mesh.nelem <= 600000:
            return AssemblySmall(fem_solver, function_spaces, formulation, mesh, material, boundary_condition, Eulerx)
        elif mesh.nelem > 600000:
            print("Larger than memory system. Dask on disk parallel assembly is turned on")
            return OutofCoreAssembly(fem_solver, function_spaces[0], formulation, mesh, material, Eulerx)

def AssemblySmall(fem_solver, function_spaces, formulation, mesh, material, boundary_condition, Eulerx):

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

    mass = []

    if fem_solver.parallel:
        # COMPUATE ALL LOCAL ELEMENTAL MATRICES (STIFFNESS, MASS, INTERNAL & EXTERNAL TRACTION FORCES)
        ParallelTuple = parmap.map(formulation,np.arange(0,nelem,dtype=np.int32),
            function_spaces[0], mesh, material, fem_solver, Eulerx, processes= int(multiprocessing.cpu_count()/2))

    for elem in range(nelem):

        if fem_solver.parallel:
            # UNPACK PARALLEL TUPLE VALUES
            I_stiff_elem = ParallelTuple[elem][0]; J_stiff_elem = ParallelTuple[elem][1]; V_stiff_elem = ParallelTuple[elem][2]
            t = ParallelTuple[elem][3]; f = ParallelTuple[elem][4]
            I_mass_elem = ParallelTuple[elem][5]; J_mass_elem = ParallelTuple[elem][6]; V_mass_elem = ParallelTuple[elem][6]

        else:
            # COMPUATE ALL LOCAL ELEMENTAL MATRICES (STIFFNESS, MASS, INTERNAL TRACTION FORCES )
            I_stiff_elem, J_stiff_elem, V_stiff_elem, t, \
            I_mass_elem, J_mass_elem, V_mass_elem = formulation.GetElementalMatrices(elem,
                function_spaces[0], mesh, material, fem_solver, Eulerx)

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

    if fem_solver.has_moving_boundary:
        K_pressure, f_pressure = AssemblyFollowerForces(boundary_condition, mesh, material, function_spaces, Eulerx, fem_solver)
        stiffness -= K_pressure
        T -= f_pressure

    fem_solver.assembly_time = time() - t_assembly
    return stiffness, T, mass

#----------------------------------------------------------------------------------------------------------------#
#------------------------------- ASSEMBLY ROUTINE FOR EXTERNAL TRACTION FORCES ----------------------------------#
#----------------------------------------------------------------------------------------------------------------#

def AssemblyFollowerForces(boundary_condition, mesh, material, function_spaces, Eulerx, fem_solver):
    """Compute/assemble traction (follower)"""

    if boundary_condition.pressure_flags is None:
        raise RuntimeError("Pressure boundary conditions are not set for the analysis")

    nvar = material.nvar
    ndim = mesh.InferSpatialDimension()

    if boundary_condition.pressure_flags.shape[0] == mesh.points.shape[0]:
        boundary_condition.pressure_data_applied_at = "node"
        raise ValueError("Follower forces applied at nodes")

    if not isinstance(function_spaces,tuple):
        raise ValueError("Boundary functional spaces not available for computing pressure stiffness")
    else:
        # CHECK IF A FUNCTION SPACE FOR BOUNDARY EXISTS - SAFEGAURDS AGAINST FORMULATIONS THAT DO NO PROVIDE ONE
        has_boundary_spaces = False
        for fs in function_spaces:
            if ndim == 3 and fs.ndim == 2:
                has_boundary_spaces = True
                break
            elif ndim == 2 and fs.ndim == 1:
                has_boundary_spaces = True
                break
        if not has_boundary_spaces:
            from Kuru import QuadratureRule, FunctionSpace
            # COMPUTE BOUNDARY FUNCTIONAL SPACES
            p = mesh.InferPolynomialDegree()
            bquadrature = QuadratureRule(optimal=3, norder=2*p+1,
                mesh_type=mesh.boundary_element_type, is_flattened=False)
            bfunction_space = FunctionSpace(mesh.CreateDummyLowerDimensionalMesh(),
                bquadrature, p=p, equally_spaced=mesh.IsEquallySpaced, use_optimal_quadrature=False)
            function_spaces = (function_spaces[0],bfunction_space)

    local_size = function_spaces[-1].Bases.shape[0]*nvar
    boundary_condition.local_rows = np.repeat(np.arange(0,local_size),local_size,axis=0)
    boundary_condition.local_columns = np.tile(np.arange(0,local_size),local_size)
    boundary_condition.local_size = local_size

    if boundary_condition.analysis_type == "static":

        if ndim == 2:
            faces = mesh.edges
            nodeperface = mesh.edges.shape[1]
        else:
            faces = mesh.faces
            nodeperface = mesh.faces.shape[1]

        ngauss = function_spaces[-1].AllGauss.shape[0]
        ndof = nodeperface*nvar
        local_capacity = ndof*ndof

        if fem_solver.recompute_sparsity_pattern is False:
            indices, indptr = fem_solver.indices, fem_solver.indptr
            if fem_solver.squeeze_sparsity_pattern is False:
                data_global_indices = fem_solver.data_global_indices
                data_local_indices = fem_solver.data_local_indices

        if fem_solver.recompute_sparsity_pattern:
            # ALLOCATE VECTORS FOR SPARSE ASSEMBLY OF STIFFNESS MATRIX - CHANGE TYPES TO INT64 FOR DoF > 1e09
            I_stiffness=np.zeros(int((nvar*nodeperface)**2*faces.shape[0]),dtype=np.int32)
            J_stiffness=np.zeros(int((nvar*nodeperface)**2*faces.shape[0]),dtype=np.int32)
            V_stiffness=np.zeros(int((nvar*nodeperface)**2*faces.shape[0]),dtype=np.float64)
        else:
            V_stiffness=np.zeros(indices.shape[0],dtype=np.float64)

        F = np.zeros((mesh.points.shape[0]*nvar,1))

        for face in range(faces.shape[0]):
            if boundary_condition.pressure_flags[face] == True:
                # Compute stiffness_matrix and force_vector
                I_stiff_face, J_stiff_face, V_stiff_face, f = GetFaceMatrices(face,function_spaces[-1],faces,boundary_condition,material,Eulerx)

                if fem_solver.recompute_sparsity_pattern:
                    # SPARSE ASSEMBLY - STIFFNESS MATRIX
                    SparseAssemblyNative(I_stiff_face,J_stiff_face,V_stiff_face,I_stiffness,J_stiffness,V_stiffness,
                        face,nvar,nodeperface,faces)
                else:
                    if fem_solver.squeeze_sparsity_pattern:
                        # SPARSE ASSEMBLY - STIFFNESS MATRIX
                        raise ValueError("This way of follower load is not well done yet")
                        #SparseAssemblyNativeCSR_RecomputeDataIndex(mesh,V_stiff_face,indices,indptr,V_tangent,face,ndim)
                    else:
                        # SPARSE ASSEMBLY - STIFFNESS MATRIX
                        V_stiffness[data_global_indices[face*local_capacity:(face+1)*local_capacity]] \
                            += V_stiff_face[data_local_indices[face*local_capacity:(face+1)*local_capacity]]

                RHSAssemblyNative(F,np.ascontiguousarray(f[:,None]),face,nvar,nodeperface,faces)

        if fem_solver.recompute_sparsity_pattern:
            stiffness = coo_matrix((V_stiffness,(I_stiffness,J_stiffness)),
                shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])),dtype=np.float64).tocsr()
        else:
            stiffness = csr_matrix((V_stiffness,indices,indptr),
                shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])))

    elif boundary_condition.analysis_type == "dynamic":
        raise ValueError("Not implemented yet")

    return stiffness, F

def GetFaceMatrices(face, function_space, faces, boundary_condition, material, Eulerx):

    nvar = material.nvar
    ndim = material.ndim
    ngauss = function_space.AllGauss.shape[0]

    nodeperface = faces.shape[1]
    EulerElemCoords = Eulerx[faces[face,:],:]

    alternating = np.zeros((3,3,3),dtype=np.float64)
    alternating[0,1,2]=1.
    alternating[2,0,1]=1.
    alternating[1,2,0]=1.
    alternating[2,1,0]=-1.
    alternating[0,2,1]=-1.
    alternating[1,0,2]=-1.

    N = np.zeros((nodeperface*ndim,ndim,ngauss))
    gNx = np.zeros((nodeperface*ndim,ndim,ngauss))
    gNy = np.zeros((nodeperface*ndim,ndim,ngauss))
    for i in range(ndim):
        N[i::ndim,i,:] = function_space.Bases
        gNx[i::ndim,i,:] = function_space.gBasesx
        gNy[i::ndim,i,:] = function_space.gBasesy

    # mapping tangential vector [\partial\vec{x}/ \partial\zeta (ngauss x ndim)]
    tangential1 = np.einsum("ij,ik->jk",function_space.gBasesx,EulerElemCoords)
    # mapping tangential vector [\partial\vec{x}/ \partial\eta (ngauss x ndim)]
    tangential2 = np.einsum("ij,ik->jk",function_space.gBasesy,EulerElemCoords)
    # mapping normal (ngauss x ndim)
    normal = np.einsum("ijk,lj,lk->li",alternating,tangential1,tangential2)

    # Gauss quadrature of follower load (traction)
    force = np.einsum("ijk,kj,k->ik",N,normal,function_space.AllGauss[:,0]).sum(axis=1)
    force = boundary_condition.pressure_increment*boundary_condition.applied_pressure[face]*force

    # Gauss quadrature of follower load (stiffness)
    cross1 = np.einsum("ijk,ljm,nkm->lnim",alternating,gNy,N)-np.einsum("ijk,ljm,nkm->lnim",alternating,N,gNy)
    cross2 = np.einsum("ijk,ljm,nkm->lnim",alternating,gNx,N)-np.einsum("ijk,ljm,nkm->lnim",alternating,N,gNx)
    quadrature1 = np.einsum("ij,klji,i->kli",tangential1,cross1,function_space.AllGauss[:,0]).sum(axis=2)
    quadrature2 = np.einsum("ij,klji,i->kli",tangential2,cross2,function_space.AllGauss[:,0]).sum(axis=2)
    stiffnessfa = boundary_condition.pressure_increment*0.5*boundary_condition.applied_pressure[face]*(quadrature1-quadrature2)

    I_stiff_face = boundary_condition.local_rows
    J_stiff_face = boundary_condition.local_columns
    V_stiff_face = stiffnessfa.ravel()

    return I_stiff_face, J_stiff_face, V_stiff_face, force

def AssembleFollowerForces(boundary_condition, mesh, material, function_space, Eulerx):

    nvar = material.nvar
    ndim = material.ndim
    ngauss = function_space.AllGauss.shape[0]

    if ndim == 2:
        faces = mesh.edges
        nodeperface = mesh.edges.shape[1]
    else:
        faces = mesh.faces
        nodeperface = mesh.faces.shape[1]

    if boundary_condition.is_applied_pressure_shape_functions_computed is False:
        N = np.zeros((nodeperface*ndim,ndim,ngauss))
        for i in range(ndim):
            N[i::ndim,i,:] = function_space.Bases
        boundary_condition.__Nt__ = N
        boundary_condition.is_applied_pressure_shape_functions_computed = True
    else:
        N = boundary_condition.__Nt__

    F = np.zeros((mesh.points.shape[0]*ndim,1))
    for face in range(faces.shape[0]):
        if boundary_condition.neumann_flags[face] == True:
            # mapping tangential vector [\partial\vec{x}/ \partial\zeta (ngauss x ndim)]
            tangential1 = np.einsum("ij,ik->jk",function_space.gBasesx,Eulerx[faces[face,:],:])
            # mapping tangential vector [\partial\vec{x}/ \partial\eta (ngauss x ndim)]
            tangential2 = np.einsum("ij,ik->jk",function_space.gBasesy,Eulerx[faces[face,:],:])
            # mapping normal (ngauss x ndim)
            normal = np.cross(tangential1,tangential2)
            # Gauss quadrature of follower load (traction)
            f = np.einsum("ijk,kj,k->ik",N,normal,function_space.AllGauss[:,0]).sum(axis=1)
            f = boundary_condition.applied_neumann[face]*f
            RHSAssemblyNative(F,np.ascontiguousarray(f[:,None]),face,nvar,nodeperface,faces)

    return F

def AssembleForces(boundary_condition, mesh, material, function_spaces,
        compute_traction_forces=True, compute_body_forces=False):

    Ft = np.zeros((mesh.points.shape[0]*material.nvar,1))
    Fb = np.zeros((mesh.points.shape[0]*material.nvar,1))

    if compute_body_forces:
        Fb = AssembleBodyForces(boundary_condition, mesh, material, function_spaces[0])
    if compute_traction_forces:
        Ft = AssembleExternalTractionForces(boundary_condition, mesh, material, function_spaces[-1])

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

