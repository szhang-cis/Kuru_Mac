import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"),"florence"))
import numpy as np
from Florence import *


def dielectric_wrinkling(recompute_sparsity_pattern=True, squeeze_sparsity_pattern=False):
    """ Implicit quasi-static analysis of large deformation in a soft dielectric elastomer
        undergoing potential wrinkling using the couple electromechanics formulation
    """

    # Create a cylindrical disc
    radius = 20
    mesh = Mesh()
    # mesh.Cylinder(radius=radius,length=0.1,nlong=1, nrad=15, ncirc=30)
    mesh.Cylinder(radius=radius,length=0.1,nlong=1, nrad=20, ncirc=45)

    # Material constants
    e0 = 8.8541e-12
    mu = 1.0e5
    mu1 = mu
    mu2 = 0.
    eps_2 = 4.0*e0
    v = 0.4
    lamb = 2.*mu*v/(1-2.*v)

    # Use one of the ideal dielectric models
    material = IsotropicElectroMechanics_108(3, mu1=mu1, mu2=mu2, lamb=lamb, eps_2=eps_2, rho=1000.)


    def DirichletFunc(mesh):

        boundary_data = np.zeros((mesh.points.shape[0],4))+np.NAN

        # Constrain (mechanically) the perimeter of disc at the base
        r = np.sqrt(mesh.points[:,0]**2 + mesh.points[:,1]**2)
        Z_0 = np.logical_and(np.isclose(r,radius),np.isclose(mesh.points[:,2],0.))
        boundary_data[Z_0,:3] = 0.

        # Closed circuit condition [electric potential dofs]
        Z_0 = np.isclose(mesh.points[:,2],0.)
        boundary_data[Z_0,3] = 0.
        Z_0 = np.isclose(mesh.points[:,2],mesh.points[:,2].max())
        boundary_data[Z_0,3] = 5e6

        return boundary_data

    boundary_condition = BoundaryCondition()
    boundary_condition.SetDirichletCriteria(DirichletFunc, mesh)

    formulation = DisplacementPotentialFormulation(mesh)

    fem_solver = FEMSolver(number_of_load_increments=50,
        analysis_nature="nonlinear",
        analysis_type="static",
        newton_raphson_tolerance=1e-5,
        maximum_iteration_for_newton_raphson=200,
        optimise=True,
        recompute_sparsity_pattern=recompute_sparsity_pattern,
        squeeze_sparsity_pattern=squeeze_sparsity_pattern,
        print_incremental_log=True
        )

    results = fem_solver.Solve(formulation=formulation, mesh=mesh,
            material=material, boundary_condition=boundary_condition)

    # Check results
    norm = lambda s: np.linalg.norm(s[:,:3,:])
    assert norm(results.GetSolutionVectors()) > 900.
    assert norm(results.GetSolutionVectors()) < 910.

    # Plot the deformation process - requires mayavi
    results.Plot(quantity=0, configuration='deformed')



if __name__ == "__main__":
    dielectric_wrinkling()
    dielectric_wrinkling(recompute_sparsity_pattern=False, squeeze_sparsity_pattern=False)
    dielectric_wrinkling(recompute_sparsity_pattern=False, squeeze_sparsity_pattern=True)
