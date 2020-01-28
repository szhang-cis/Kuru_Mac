#cython: profile=False
#cython: infer_types=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False

import numpy as np
cimport numpy as np

ctypedef double Real


cdef extern from "_IncompressibleAnisotropicFungQuadratic_.h" nogil:
    cdef cppclass _IncompressibleAnisotropicFungQuadratic_[Real]:
        _IncompressibleAnisotropicFungQuadratic_() except +
        _IncompressibleAnisotropicFungQuadratic_(Real mu, Real k1, Real k2) except +
        void SetParameters(Real mu, Real k1, Real k2) except +
        void KineticMeasures(Real *Snp, Real* Hnp, int ndim, int ngauss, const Real *Fnp, int nfibre, const Real *Nnp, const Real pressure) except +



def KineticMeasures(material, np.ndarray[Real,ndim=3,mode='c'] F, np.ndarray[Real,ndim=2,mode='c'] N):
    
    cdef int ndim = F.shape[2]
    cdef int ngauss = F.shape[0]
    cdef int nfibre = N.shape[0]
    cdef np.ndarray[Real, ndim=3, mode='c'] stress, hessian

    stress = np.zeros((ngauss,ndim,ndim),dtype=np.float64)
    if ndim==3:
        hessian = np.zeros((ngauss,6,6),dtype=np.float64)
    elif ndim==2:
        hessian = np.zeros((ngauss,3,3),dtype=np.float64)

    cdef Real pressure = material.pressure

    cdef _IncompressibleAnisotropicFungQuadratic_[Real] mat_obj = _IncompressibleAnisotropicFungQuadratic_()
    mat_obj.SetParameters(material.mu,material.k1,material.k2)
    mat_obj.KineticMeasures(&stress[0,0,0], &hessian[0,0,0], ndim, ngauss, &F[0,0,0], nfibre, &N[0,0], pressure)

    return stress, hessian
