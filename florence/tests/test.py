from __future__ import print_function
import sys, os, imp, time, gc, importlib
from sys import exit
from datetime import datetime
from warnings import warn

import numpy as np
import scipy as sp
# from scipy.io import loadmat, savemat

# AVOID WRITING .pyc OR .pyo FILES
sys.dont_write_bytecode
# SET NUMPY'S LINEWIDTH PRINT OPTION
np.set_printoptions(linewidth=300)


# APPEND ALL REQUIRED PATHS
sys.path.append(os.path.join(os.path.expanduser("~"),"florence"))

sys.path.append('../examples/simple_laplace/')
sys.path.append('../examples/linear_elastic_dynamics/')
sys.path.append('../examples/curved_mesh_generation/')
sys.path.append('../examples/car_crash_analysis/')
sys.path.append('../examples/hyperelastic_explicit_dynamics/')
sys.path.append('../examples/electro_hyperelastic_explicit_dynamics')
sys.path.append('../examples/wrinkling_of_soft_dielectric_film/')
sys.path.append('../examples/staggered_multiphysics_solver')
sys.path.append('../examples/mixed_fem_multiphysics_strain_gradient_solvers')
sys.path.append('../examples/nonlinear_elastodynamics')
sys.path.append('../examples/nonlinear_electroelastodynamics')


sys.path.append('./test_basics')
sys.path.append('./test_BEM')


# IMPORT FLORENCE
from Florence import *

# IMPORT EXAMPLES
from simple_laplace import simple_laplace
from crash_analysis_with_explicit_contact import crash_analysis
from high_order_curved_mesh_generation import high_order_curved_mesh_generation
from hyperelastic_explicit_dynamics import explicit_dynamics_mechanics
from wrinkling_of_soft_dielectric_film import dielectric_wrinkling
from linear_elastic_dynamics import linear_elastic_dynamics
from electro_hyperelastic_explicit_dynamics import electro_hyperelastic_explicit_dynamics
from staggered_multiphysics_solver import staggered_multiphysics_solver
from mixed_fem_multiphysics_strain_gradient_solvers import strain_gradient_elastodynamics
from mixed_fem_multiphysics_strain_gradient_solvers import strain_gradient_electroelastodynamics
from nonlinear_elastodynamics import nonlinear_elastodynamics
from nonlinear_electroelastodynamics import nonlinear_electroelastodynamics


from test_basics import test_quadrature_functionspace, test_mesh_postprocess_material, test_material
from test_BEM import test_BEM

tick  = u'\u2713'.encode('utf8')  + b' : '
cross = u'\u2717'.encode('utf8')  + b' : '


class run_tests(object):
    def __init__(self):
        pass


def entity_checker(x,y,tol=1e-10):
    """x,y both being ndarrays, x being new solution and y being pre-computed solution.
        tol parameter is a very important aspect of the comparison"""

    if x.shape != y.shape:
        raise TypeError("shape of entity does not match with pre-computed solution")

    # safe-guard against unsigned overflow
    # if x.dtype != y.dtype:
    #     if x.dtype != float:
    #         if x.dtype == np.uint64 or x.dtype == np.uint32:
    #             x = x.astype(np.int64)
    #     if y.dtype != float:
    #         if y.dtype == np.uint64 or y.dtype == np.uint32:
    #             y = y.astype(np.int64)

    if x.dtype == np.uint64 or x.dtype == np.uint32:
        x = x.astype(np.int64)
    if y.dtype == np.uint64 or y.dtype == np.uint32:
        y = y.astype(np.int64)

    # print(np.sum(x[:,:4]-y[:,:4]))

    if np.isclose(x-y,0.,atol=tol).all():
        return True
    elif np.isclose(np.sum(x-y),0.,atol=tol):
        return True
    else:
        return False


def mesh_checker(mesh,Dict):
    """Give a mesh and a Dict loaded from HDF5 to compare"""

    # print((mesh.elements - Dict['elements']).max())
    # print(mesh.elements.dtype, Dict['elements'].dtype)
    print("Checking higher order mesh generators results")
    if entity_checker(mesh.elements,Dict['elements']):
        print(tick, "mesh elements match")
    else:
        print(cross, "mesh elements do not match")
        exit()

    if entity_checker(mesh.points,Dict['points']):
        print(tick, "mesh points match")
    else:
        print(cross, "mesh points do not match")
        exit()

    if entity_checker(mesh.edges,Dict['edges']):
        print(tick, "mesh edges match")
    else:
        print(cross, "mesh edges do not match")
        exit()

    if mesh.element_type == "tet" or mesh.element_type == "hex":
        if entity_checker(mesh.faces,Dict['faces']):
            print(tick, "mesh faces match")
        else:
            print(cross, "mesh faces do not match")
            exit()


def dirichlet_checker(ColumnsOut,AppliedDirichlet,Dict):

    if ColumnsOut.ndim > 1:
        ColumnsOut = ColumnsOut.flatten()
    if Dict['ColumnsOut'].ndim > 1:
        Dict['ColumnsOut'] = Dict['ColumnsOut'].flatten()

    if AppliedDirichlet.ndim > 1:
        AppliedDirichlet = AppliedDirichlet.flatten()
    if Dict['AppliedDirichlet'].ndim > 1:
        Dict['AppliedDirichlet'] = Dict['AppliedDirichlet'].flatten()

    print("Checking for projection data from OpenCascade wrapper")
    if entity_checker(ColumnsOut,Dict['ColumnsOut']):
        print(tick, "Dirichlet degrees of freedom match")
    else:
        print(cross, "Dirichlet degrees of freedom do not match")
        exit()

    if entity_checker(AppliedDirichlet,Dict['AppliedDirichlet']):
        print(tick, "Dirichlet data for degrees of freedom match")
    else:
        print(cross, "Dirichlet data for degrees of freedom do not match")
        exit()

def final_solution_checker(material,solver,fem_solver,TotalDisp,Dict):

    print("Checking for final solution")
    if not np.isclose(material.nu, float(Dict['PoissonRatio'])):
        raise ValueError("Analysis with different material parameters are being compared")
    if not np.isclose(material.E,float(Dict['YoungsModulus'])):
        raise ValueError("Analysis with different material parameters are being compared")
    if material.is_transversely_isotropic:
        if not np.isclose(material.E_A,float(Dict['E_A'])):
            raise ValueError("Analysis with different material parameters are being compared")
        if not np.isclose(material.G_A,float(Dict['G_A'])):
            raise ValueError("Analysis with different material parameters are being compared")

    if solver.solver_type != Dict['SolverType']:
        raise ValueError("Results from different solvers are being compared")
    elif solver.solver_type == "multigrid" or solver.solver_type == "amg":
        if solver.solver_subtype == "multigrid" or solver.solver_subtype == "amg":
            if solver.iterative_solver_tolerance != Dict['SolverTol']:
                raise ValueError("Solver results with different tolerances are being compared")


    tol = 1e-05
    if entity_checker(TotalDisp,Dict['TotalDisp'],tol):
        print(tick, "Final solution is correct")
    else:
        print(np.linalg.norm(TotalDisp - Dict['TotalDisp']))
        print(cross, "Final solution does not match with pre-computed solution")
        exit()

    Dict['ScaledJacobian'] = Dict['ScaledJacobian'].flatten()
    fem_solver.ScaledJacobian = fem_solver.ScaledJacobian.flatten()
    if np.abs((fem_solver.ScaledJacobian.min() - Dict['ScaledJacobian'].min())<tol):
        print(tick,"Final mesh quality is correct")
    else:
        print(cross,"Final mesh quality does not match")
        exit()



def test_examples():

    print("Statrting to run the test suite")

    # RUN EXAMPLES AT TEST CASES
    simple_laplace(optimise=False, recompute_sparsity_pattern=True, squeeze_sparsity_pattern=False)
    simple_laplace(optimise=False, recompute_sparsity_pattern=False, squeeze_sparsity_pattern=False)
    simple_laplace(optimise=False, recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)
    simple_laplace(optimise=True, recompute_sparsity_pattern=True, squeeze_sparsity_pattern=False)
    simple_laplace(optimise=True, recompute_sparsity_pattern=False, squeeze_sparsity_pattern=False)
    simple_laplace(optimise=True, recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)

    high_order_curved_mesh_generation(p=2, analysis_nature="linear", optimise=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="linear", optimise=True,
        recompute_sparsity_pattern=False)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=True,
        recompute_sparsity_pattern=False)
    high_order_curved_mesh_generation(p=2, analysis_nature="linear", optimise=True,
        recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=True,
        recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="linear", optimise=False)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=False)
    high_order_curved_mesh_generation(p=2, analysis_nature="linear", optimise=False,
        recompute_sparsity_pattern=False)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=False,
        recompute_sparsity_pattern=False)
    high_order_curved_mesh_generation(p=2, analysis_nature="linear", optimise=False,
        recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=False,
        recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=False, parallelise=True)
    high_order_curved_mesh_generation(p=2, analysis_nature="nonlinear", optimise=True, parallelise=True)

    linear_elastic_dynamics()
    crash_analysis()
    explicit_dynamics_mechanics()

    electro_hyperelastic_explicit_dynamics()
    electro_hyperelastic_explicit_dynamics(recompute_sparsity_pattern=False, squeeze_sparsity_pattern=False)
    electro_hyperelastic_explicit_dynamics(recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)

    dielectric_wrinkling()
    dielectric_wrinkling(recompute_sparsity_pattern=False, squeeze_sparsity_pattern=False)
    dielectric_wrinkling(recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)

    staggered_multiphysics_solver()
    strain_gradient_elastodynamics()
    strain_gradient_electroelastodynamics()
    nonlinear_elastodynamics(optimise=False)
    nonlinear_elastodynamics(optimise=True)
    nonlinear_electroelastodynamics(optimise=False)
    nonlinear_electroelastodynamics(optimise=True)

    # RUN BASICS TESTSUITE
    test_quadrature_functionspace()
    test_mesh_postprocess_material()
    test_material()
    test_BEM()

    print("Successfully finished running all tests")


# RUN EXAPLES AS TEST CASES
if __name__ == "__main__":
    test_examples()
