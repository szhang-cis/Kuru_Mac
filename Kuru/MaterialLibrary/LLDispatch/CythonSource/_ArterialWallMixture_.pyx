#cython: profile=False
#cython: infer_types=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False

import numpy as np
cimport numpy as np

ctypedef double Real


cdef extern from "_ArterialWallMixture_.h" nogil:
    cdef cppclass _ArterialWallMixture_[Real]:
        _ArterialWallMixture_() except +
        _ArterialWallMixture_(Real mu, Real kappa, Real k1m, Real k2m, Real k1c, Real k2c) except +
        void SetParameters(Real mu, Real kappa, Real k1m, Real k2m, Real k1c, Real k2c) except +
        void KineticMeasures(Real *Snp, Real* Hnp, int ndim, int ngauss, const Real *Fnp, int nfibre,
                const Real *Nnp, int nfield, const Real *FVnp) except +



def KineticMeasures(material, np.ndarray[Real,ndim=3,mode='c'] F,
        np.ndarray[Real,ndim=2,mode='c'] anisotropic_orientations,
        np.ndarray[Real,ndim=2,mode='c'] field_variables):
    
    cdef int ndim = F.shape[2]
    cdef int ngauss = F.shape[0]
    cdef int nfibre = anisotropic_orientations.shape[0]
    cdef int nfield = field_variables.shape[1]
    cdef np.ndarray[Real, ndim=3, mode='c'] stress, hessian

    stress = np.zeros((ngauss,ndim,ndim),dtype=np.float64)
    if ndim==3:
        hessian = np.zeros((ngauss,6,6),dtype=np.float64)
    elif ndim==2:
        hessian = np.zeros((ngauss,3,3),dtype=np.float64)

    cdef _ArterialWallMixture_[Real] mat_obj = _ArterialWallMixture_()
    mat_obj.SetParameters(material.mu,material.kappa,material.k1m,material.k2m,material.k1c,material.k2c)
    mat_obj.KineticMeasures(&stress[0,0,0], &hessian[0,0,0], ndim, ngauss, &F[0,0,0], nfibre, 
            &anisotropic_orientations[0,0], nfield, &field_variables[0,0])

    return stress, hessian
