import numpy as np
cimport numpy as np
from libc.stdint cimport int64_t, uint64_t
from libc.stdlib cimport malloc, free

ctypedef int64_t Integer
ctypedef uint64_t UInteger
ctypedef double Real

from .SparseAssemblyNative import SparseAssemblyNative, SparseAssemblyNativeCSR_RecomputeDataIndex #, SparseAssemblyNativeCSR
from .RHSAssemblyNative import RHSAssemblyNative


def StaticForces(boundary_condition, mesh, material, function_space, fem_solver, Real[:,::1] Eulerx):

    # GET VARIABLES FOR DISPATCHING TO C
    cdef np.ndarray[int,ndim=1, mode='c'] pressure_flags    = boundary_condition.pressure_flags
    cdef np.ndarray[Real,ndim=1, mode='c'] applied_pressure = boundary_condition.applied_pressure
    cdef Real pressure_increment                            = boundary_condition.pressure_increment

    cdef Integer ndim    = material.ndim
    cdef Integer nvar    = material.nvar

    cdef np.ndarray[Real,ndim=2, mode='c'] Bases         = function_space.Bases
    cdef np.ndarray[Real,ndim=2, mode='c'] gBasesx       = function_space.gBasesx
    cdef np.ndarray[Real,ndim=2, mode='c'] gBasesy       = function_space.gBasesy
    cdef np.ndarray[Real,ndim=1, mode='c'] AllGauss      = function_space.AllGauss.flatten()
    cdef Integer ngauss                                  = function_space.AllGauss.shape[0]
    cdef Integer local_size                              = function_space.Bases.shape[0]*nvar

    cdef int squeeze_sparsity_pattern       = fem_solver.squeeze_sparsity_pattern
    cdef int recompute_sparsity_pattern     = fem_solver.recompute_sparsity_pattern

    cdef np.ndarray[UInteger,ndim=2, mode='c'] faces = np.zeros((1,1),np.uint64)
    cdef Integer nface                               = 0
    cdef Integer nodeperface                         = 0

    cdef np.ndarray[int,ndim=1,mode='c'] I_stiff    = np.zeros(1,np.int32)
    cdef np.ndarray[int,ndim=1,mode='c'] J_stiff    = np.zeros(1,np.int32)
    cdef np.ndarray[Real,ndim=1,mode='c'] V_stiff   = np.zeros(1,np.float64)

    cdef np.ndarray[int,ndim=1,mode='c'] data_global_indices   = np.zeros(1,np.int32)
    cdef np.ndarray[int,ndim=1,mode='c'] data_local_indices    = np.zeros(1,np.int32)
    cdef np.ndarray[Integer,ndim=2, mode='c'] sorter           = np.zeros((1,1),np.int64)
    cdef np.ndarray[UInteger,ndim=2, mode='c'] sorted_elements = np.zeros((1,1),np.uint64)

    cdef np.ndarray[Real,ndim=1, mode='c'] F = np.zeros(mesh.points.shape[0]*nvar,np.float64)

    if ndim == 2:
        faces = mesh.edges
        nface = mesh.edges.shape[0]
        nodeperface = mesh.edges.shape[1]
    else:
        faces = mesh.faces
        nface = mesh.faces.shape[0]
        nodeperface = mesh.faces.shape[1]

    if fem_solver.recompute_sparsity_pattern:
        # ALLOCATE VECTORS FOR SPARSE ASSEMBLY OF STIFFNESS MATRIX - CHANGE TYPES TO INT64 FOR DoF > 1e09
        I_stiff = np.zeros(int((nvar*nodeperface)**2*nface),np.int32)
        J_stiff = np.zeros(int((nvar*nodeperface)**2*nface),np.int32)
        V_stiff = np.zeros(int((nvar*nodeperface)**2*nface),np.float64)
    else:
        I_stiff = fem_solver.indices
        J_stiff = fem_solver.indptr
        V_stiff = np.zeros(indices.shape[0],dtype=np.float64)
        data_global_indices = fem_solver.data_global_indices
        data_local_indices = fem_solver.data_local_indices
        if fem_solver.squeeze_sparsity_pattern:
            sorter = mesh.element_sorter
            sorted_elements = mesh.sorted_elements

    AssemblerStaticForces(  &faces[0,0],
                            &Eulerx[0,0],
                            &Bases[0,0],
                            &gBasesx[0,0],
                            &gBasesy[0,0],
                            &AllGauss[0],
                            &pressure_flags[0],
                            &applied_pressure[0],
                            pressure_increment,
                            recompute_sparsity_pattern,
                            squeeze_sparsity_pattern,
                            &data_global_indices[0],
                            &data_local_indices[0],
                            &sorted_elements[0,0],
                            &sorter[0,0],
                            &I_stiff[0],
                            &J_stiff[0],
                            &V_stiff[0],
                            &F[0],
                            nface,
                            nodeperface,
                            ngauss,
                            local_size,
                            nvar,
                            ndim)

    if fem_solver.recompute_sparsity_pattern:
        return I_stiff, J_stiff, V_stiff, F
    else:
        return V_stiff, F

cdef StaticForcesAssembler(const UInteger *faces,
                           const Real *Eulerx,
                           const Real *Bases,
                           const Real *gBasesx,
                           const Real *gBasesy,
                           const Real *AllGauss,
                           const int *pressure_flags,
                           const Real *applied_pressure,
                           const Real pressure_increment,
                           int recompute_sparsity_pattern,
                           int squeeze_sparsity_pattern,
                           const int *data_global_indices,
                           const int *data_local_indices,
                           cosnt UInteger *sorted_elements,
                           const Integer *sorter,
                           int *I_stiff,
                           int *J_stiff,
                           Real *V_stiff,
                           Real *F,
                           Integer nface,
                           Integer nodeperface,
                           Integer ngauss,
                           Integer local_size,
                           Integer nvar,
                           Integer ndim):

    cdef Integer ndof = nodeperface*nvar
    cdef Integer local_capacity = ndof*ndof

    # usar calloc para iniciar en zeros
    cdef Real *N = <Real*>malloc(sizeof(Real)*ndof*nvar*ngauss)
    cdef Real *gNx = <Real*>malloc(sizeof(Real)*ndof*nvar*ngauss)
    cdef Real *gNy = <Real*>malloc(sizeof(Real)*ndof*nvar*ngauss)
    cdef Real *alternating = <Real*>malloc(sizeof(Real)*3*3*3) #alternating[i*3*3 + j*3 + k]
    alternating[0*3*3 + 1*3 + 2]=1.                            #              j k     k
    alternating[0*3*3 + 2*3 + 1]=-1.
    alternating[1*3*3 + 0*3 + 2]=-1.
    alternating[1*3*3 + 2*3 + 0]=1.
    alternating[2*3*3 + 0*3 + 1]=1.
    alternating[2*3*3 + 1*3 + 0]=-1.

    #0::3=>0,3,6,9 #1::3=>1,4,7,10 #2::3=>2,5,8,11
    cdef Integer ivar
    cdef Integer inode
    cdef Integer idof
    cdef Integer igauss
    for ivar in range(nvar):
        for inode in range(nodeperface):
            idof = inode*nvar + ivar
            for igauss in range(ngauss):
                N[idof*nvar*ngauss + ivar*ngauss + igauss] = Bases[inode*ngauss + igauss]
                gNx[idof*nvar*ngauss + ivar*ngauss + igauss] = gBasesx[inode*ngauss + igauss]
                gNy[idof*nvar*ngauss + ivar*ngauss + igauss] = gBasesy[inode*ngauss + igauss]

    # Compute stiffness_matrix and force_vector
    cdef Integer i
    cdef Integer j
    cdef Integer k
    cdef Real *EulerElemCoords = <Real*>malloc(sizeof(Real)*nodeperface*ndim)
    for face in range(nface):
        if pressure_flags[face]:
            pressure = pressure_increment*applied_pressure[face]

            for i in range(nodeperface):
                inode = faces[face*nodeperface + i]
                for j in range(ndim):
                    EulerElemCoords[i*ndim + j] = Eulerx[inode*ndim + j]

            # mapping tangential vectors [\partial\vec{x}/ \partial\zeta (ngauss x ndim)]
            tangential1 = np.einsum("ij,ik->jk",gBasesx,EulerElemCoords)
            tangential2 = np.einsum("ij,ik->jk",gBasesy,EulerElemCoords)
            for j in range(ngauss):
                for k in range(ndim):
                    for i in range(nodeperface):
                        tangential1[j*ndim + k] += gBasesx[i*ngauss + j]*EulerElemCoords[i*ndim + k]
                        tangential2[j*ndim + k] += gBasesy[i*ngauss + j]*EulerElemCoords[i*ndim + k]
            # mapping normal (ngauss x ndim)
            normal = np.einsum("ijk,lj,lk->li",alternating,tangential1,tangential2)
            for l in range(ngauss):
                for i in range(ndim):
                    for j in range(ndim):
                        for k in range(ndim):
                            normal[l*ndim+i] += alternating[i*ndim*ndim+j*ndim+k]*tangential1[l*ndim+j]*tangential2[l*ndim+k]
            # Gauss quadrature of follower load (traction)
            force = np.einsum("ijk,kj,k->ik",N,normal,AllGauss[:,0]).sum(axis=1)
            for i in range(ndof):
                for k in range(ngauss):
                    for j in range(nvar):
                        force[i*ngauss+k] += N[i*nvar*ngauss+j*ngauss+k]*normal[k*nvar+j]*AllGauss[k]
            force=pressure*force
            # Gauss quadrature of follower load (stiffness)
            cross1 = np.einsum("ijk,ljm,nkm->lnim",alternating,gNy,N)-np.einsum("ijk,ljm,nkm->lnim",alternating,N,gNy)
            cross2 = np.einsum("ijk,ljm,nkm->lnim",alternating,gNx,N)-np.einsum("ijk,ljm,nkm->lnim",alternating,N,gNx)
            for l in range(ndof):
                for n in range(ndof):
                    for i in range(ndim):
                        for m in range(ngauss):
                            for j in range(ndim):
                                for k in range(ndim):
                                    cross1[l*ndof*ndim*ngauss+n*ndim*ngauss+i*ngauss+m] += 
                                        alternating[i*ndim*ndim+j*ndim+k]*gNy[l*ndim*ngauss+j*ngauss+m]*N[n*ndim*ngauss+k*ngauss+m]
                                        -alternating[i*ndim*ndim+j*ndim+k]*N[l*ndim*ngauss+j*ngauss+m]*gNy[n*ndim*ngauss+k*ngauss+m]
                                    cross2[l*ndof*ndim*ngauss+n*ndim*ngauss+i*ngauss+m] += 
                                        alternating[i*ndim*ndim+j*ndim+k]*gNx[l*ndim*ngauss+j*ngauss+m]*N[n*ndim*ngauss+k*ngauss+m]
                                        -alternating[i*ndim*ndim+j*ndim+k]*N[l*ndim*ngauss+j*ngauss+m]*gNx[n*ndim*ngauss+k*ngauss+m]
            quadrature1 = np.einsum("ij,klji,i->kli",tangential1,cross1,AllGauss[:,0]).sum(axis=2)
            quadrature2 = np.einsum("ij,klji,i->kli",tangential2,cross2,AllGauss[:,0]).sum(axis=2)
            for k in range(ndof):
                for l in range(ndof):
                    for i in range(ngauss):
                        for j in range(ndim):
                            quadrature1[k*ndof*ngauss+l*ngauss+i] += tangential1[i*ndim+j]*cross1[k*ndof*ndim*ngauss+l*ndim*ngauss+j*ngauss+i]*AllGauss[i]
                            quadrature2[k*ndof*ngauss+l*ngauss+i] += tangential2[i*ndim+j]*cross2[k*ndof*ndim*ngauss+l*ndim*ngauss+j*ngauss+i]*AllGauss[i]
            stiffnessfa = pressure*0.5*(quadrature1-quadrature2)

            I_stiff_face = np.repeat(np.arange(0,local_size),local_size,axis=0)
            J_stiff_face = np.tile(np.arange(0,local_size),local_size)
            V_stiff_face = stiffnessfa.ravel()

            if recompute_sparsity_pattern:
                # SPARSE ASSEMBLY - STIFFNESS MATRIX
                SparseAssemblyNative(I_stiff_face,J_stiff_face,V_stiff_face,I_stiff,J_stiff,V_stiff,
                    face,nvar,nodeperface,faces)
            else:
                if squeeze_sparsity_pattern:
                    # SPARSE ASSEMBLY - STIFFNESS MATRIX
                    SparseAssemblyNativeCSR_RecomputeDataIndex(mesh,V_stiff_face,I_stiff,J_stiff,V_stiff,face,ndim)
                else:
                    # SPARSE ASSEMBLY - STIFFNESS MATRIX
                    V_stiff[data_global_indices[face*local_capacity:(face+1)*local_capacity]] \
                        += V_stiff_face[data_local_indices[face*local_capacity:(face+1)*local_capacity]]

            RHSAssemblyNative(F,np.ascontiguousarray(f[:,None]),face,nvar,nodeperface,faces)



