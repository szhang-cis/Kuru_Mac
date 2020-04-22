#cython: profile=False
#cython: infer_types=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: nonecheck=False

import numpy as np
cimport numpy as np

ctypedef double Real


cdef extern from "_VolumetricStiffness_.h":
    void _VolumetricStiffnessFiller_(Real *volumetric_stiffness, Real *mean_volume,
        const Real *SpatialGradient, const Real *detJ, const Real *dV,
        const Real *density, const Real *Growth, int has_growth_remodeling,
        const int ndim, const int nvar, const int nodeperelem, const int nguass)

def _VolumetricStiffnessIntegrand_(material, np.ndarray[Real, ndim=3, mode='c'] SpatialGradient,
    np.ndarray[Real, ndim=1] detJ, np.ndarray[Real, ndim=1] dV, int nvar):

    cdef int ngauss = SpatialGradient.shape[0]
    cdef int nodeperelem = SpatialGradient.shape[1]
    cdef int ndim = SpatialGradient.shape[2]
    cdef int local_size = nvar*nodeperelem

    cdef int has_growth_remodeling                    = material.has_growth_remodeling
    cdef np.ndarray[Real, ndim=1, mode='c'] density   = np.zeros(1,dtype=np.float64)
    cdef np.ndarray[Real, ndim=1, mode='c'] Growth    = np.zeros(1,dtype=np.float64)
    if material.has_growth_remodeling:
        density = np.ascontiguousarray(material.StateVariables[:,14])
        Growth  = np.ascontiguousarray(material.StateVariables[:,20])

    cdef np.ndarray[Real, ndim=2, mode='c'] volumetric_stiffness = np.zeros((local_size,
        local_size),dtype=np.float64)
    cdef np.ndarray[Real, ndim=1, mode='c'] mean_volume = np.zeros(1,dtype=np.float64)

    _VolumetricStiffnessFiller_(&volumetric_stiffness[0,0], &mean_volume[0], &SpatialGradient[0,0,0], 
        &detJ[0], &dV[0], &density[0], &Growth[0], has_growth_remodeling, ndim, nvar, nodeperelem, ngauss)

    return volumetric_stiffness, mean_volume[0]
