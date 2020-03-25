from warnings import warn
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix

try:
    #from ._LLADF__NearlyIncompressibleNeoHookean_ import _LLADF__NearlyIncompressibleNeoHookean_
    #from ._LLADF__AnisotropicFungQuadratic_ import _LLADF__AnisotropicFungQuadratic_
    from ._LLADF__ArterialWallMixture_ import _LLADF__ArterialWallMixture_
    #from ._LLADF__IncompressibleNeoHookean_ import _LLADF__IncompressibleNeoHookean_
    #from ._LLADF__IncompressibleAnisotropicFungQuadratic_ import _LLADF__IncompressibleAnisotropicFungQuadratic_
    #from ._LLADF__IncompressibleArterialWallMixture_ import _LLADF__IncompressibleArterialWallMixture_
    has_low_level_dispatcher = True
except ImportError:
    has_low_level_dispatcher = False
    warn("Cannot use low level dispatchers for Assembly")

__all__ = ['_LowLevelAssembly_']

def _LowLevelAssembly_(fem_solver, function_space, formulation, mesh, material, Eulerx):

    prefix = "_LLADF__"

    assembly_func = prefix + type(material).__name__ + "_"
    # CHECK IF LOW LEVEL ASSEMBLY EXISTS FOR MATERIAL
    ll_exists = False
    for key in globals().keys():
        if assembly_func == key:
            ll_exists = True
            break
    if ll_exists is False:
        raise NotImplementedError("Turning optimise option on for {} material is not supported yet. "
            "Consider 'optimise=False' for now".format(type(material).__name__))

    nvar = formulation.nvar

    # MAKE MESH DATA CONTIGUOUS
    mesh.ChangeType()
    # CALL LOW LEVEL ASSEMBLER
    if fem_solver.recompute_sparsity_pattern:
        I_stiffness, J_stiffness, V_stiffness, T = eval(assembly_func)(fem_solver,
            function_space, formulation, mesh, material, Eulerx)

        stiffness = csr_matrix((V_stiffness,(I_stiffness,J_stiffness)),
            shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])),dtype=np.float64)
    else:
        V_stiffness, T = eval(assembly_func)(fem_solver,
            function_space, formulation, mesh, material, Eulerx)

        stiffness = csr_matrix((V_stiffness,fem_solver.indices,fem_solver.indptr),
            shape=((nvar*mesh.points.shape[0],nvar*mesh.points.shape[0])),dtype=np.float64)


    mass = []

    return stiffness, T, mass

