#ifndef _LLADF__ARTERIALWALLMIXTURE_H
#define _LLADF__ARTERIALWALLMIXTURE_H

#include "assembly_helper.h"
#include "_ConstitutiveStiffnessDF_.h"
#include "_ArterialWallMixture_.h"

void _GlobalAssemblyDF__ArterialWallMixture_(const Real *points,
                        const UInteger* elements,
                        const Real* Eulerx,
                        const Real* Bases,
                        const Real* Jm,
                        const Real* AllGauss,
                        Integer ndim,
                        Integer nvar,
                        Integer ngauss,
                        Integer nelem,
                        Integer nodeperelem,
                        Integer nnode,
                        Integer H_VoigtSize,
                        Integer requires_geometry_update,
                        Integer* local_rows_stiffness,
                        Integer* local_cols_stiffness,
                        int *I_stiff,
                        int *J_stiff,
                        Real *V_stiff,
                        Real *T,
                        int recompute_sparsity_pattern,
                        int squeeze_sparsity_pattern,
                        const int *data_local_indices,
                        const int *data_global_indices,
                        const UInteger *sorted_elements,
                        const Integer *sorter,
                        Real rho,
                        Real mu,
                        Real kappa,
            		Real k1m,
			Real k2m,
            		Real k1c,
			Real k2c,
			const Real *anisotropic_orientations,
			Integer nfibre,
			const Real *field_variables,
			Integer nfield
                        ) {

    Integer ndof = nvar*nodeperelem;
    Integer local_capacity = ndof*ndof;

    Real *LagrangeElemCoords        = allocate<Real>(nodeperelem*ndim);
    Real *EulerElemCoords           = allocate<Real>(nodeperelem*ndim);
    Real *ElemFieldVariables        = allocate<Real>(nodeperelem*nfield);

    Real *F                         = allocate<Real>(ngauss*ndim*ndim);
    Real *SpatialGradient           = allocate<Real>(ngauss*nodeperelem*ndim);
    Real *detJ                      = allocate<Real>(ngauss);
    Real *dV                        = allocate<Real>(ngauss);
    Real *GaussFieldVariables       = allocate<Real>(ngauss*nfield);

    Real *stress                    = allocate<Real>(ngauss*ndim*ndim);
    Real *hessian                   = allocate<Real>(ngauss*H_VoigtSize*H_VoigtSize);

    Real *traction                  = allocate<Real>(ndof);
    Real *stiffness                 = allocate<Real>(local_capacity);
    Real *geometric_stiffness       = allocate<Real>(local_capacity);

    auto mat_obj = _ArterialWallMixture_<Real>(mu,kappa,k1m,k2m,k1c,k2c);

    // LOOP OVER ELEMETNS
    for (Integer elem=0; elem < nelem; ++elem) {

        // GET THE FIELDS AT THE ELEMENT LEVEL
        for (Integer i=0; i<nodeperelem; ++i) {
            const Integer inode = elements[elem*nodeperelem+i];
            for (Integer j=0; j<ndim; ++j) {
                LagrangeElemCoords[i*ndim+j] = points[inode*ndim+j];
                EulerElemCoords[i*ndim+j] = Eulerx[inode*ndim+j];
            }
            for (Integer k=0; k<nfield; ++k) {
                ElemFieldVariables[i*nfield+k] = field_variables[inode*nfield+k];
            }
        }

        // COMPUTE KINEMATIC MEASURES
        std::fill(F,F+ngauss*ndim*ndim,0.);
        std::fill(SpatialGradient,SpatialGradient+ngauss*nodeperelem*ndim,0.);
        std::fill(detJ,detJ+ngauss,0.);
        KinematicMeasures(  SpatialGradient,
                            F,
                            detJ,
                            dV,
                            Jm,
                            AllGauss,
                            LagrangeElemCoords,
                            EulerElemCoords,
                            ngauss,
                            ndim,
                            nodeperelem,
                            requires_geometry_update
                            );

        // MAP FIELD VARIABLES FROM NODES TO GAUSS POINTS
        //GaussFieldVariables = np.einsum('ki,kj->ij',Bases,ElemFieldVariables)
        std::fill(GaussFieldVariables,GaussFieldVariables+ngauss*nfield,0.);
        for (Integer i=0; i<ngauss; ++i) {
            for (Integer j=0; j<nfield; ++j) {
                for (Integer k=0; k<nodeperelem; ++k) {
                GaussFieldVariables[i*nfield+j] += Bases[k*ngauss+i]*ElemFieldVariables[k*nfield+j];
                }
            }
        }

        // COMPUTE KINETIC MEASURES
        mat_obj.KineticMeasures(stress, hessian, ndim, ngauss, F, nfibre, &anisotropic_orientations[elem*nfibre*ndim], nfield, GaussFieldVariables);

        // COMPUTE CONSTITUTIVE STIFFNESS AND TRACTION
        std::fill(stiffness,stiffness+local_capacity,0.);
        std::fill(traction,traction+ndof,0.);
        _ConstitutiveStiffnessIntegrandDF_Filler_(
            stiffness,
            traction,
            SpatialGradient,
            stress,
            hessian,
            detJ,
            ngauss,
            nodeperelem,
            ndim,
            nvar,
            H_VoigtSize,
            requires_geometry_update);

        // COMPUTE GEOMETRIC STIFFNESS
        std::fill(geometric_stiffness,geometric_stiffness+local_capacity,0.);
        _GeometricStiffnessFiller_( geometric_stiffness,
                                    SpatialGradient,
                                    stress,
                                    detJ,
                                    ndim,
                                    nvar,
                                    nodeperelem,
                                    ngauss);

        for (Integer i=0; i<local_capacity; ++i) {
            stiffness[i] += geometric_stiffness[i];
        }

        // ASSEMBLE CONSTITUTIVE STIFFNESS
        fill_global_data(
                nullptr,
                nullptr,
                stiffness,
                I_stiff,
                J_stiff,
                V_stiff,
                elem,
                nvar,
                nodeperelem,
                elements,
                local_capacity,
                local_capacity,
                recompute_sparsity_pattern,
                squeeze_sparsity_pattern,
                data_local_indices,
                data_global_indices,
                sorted_elements,
                sorter);

        // ASSEMBLE TRACTIONS
        {
            for (Integer i = 0; i<nodeperelem; ++i) {
                UInteger T_idx = elements[elem*nodeperelem+i]*nvar;
                for (Integer iterator = 0; iterator < nvar; ++iterator) {
                    T[T_idx+iterator] += traction[i*nvar+iterator];
                }
            }
        }

    }

    deallocate(LagrangeElemCoords);
    deallocate(EulerElemCoords);

    deallocate(F);
    deallocate(SpatialGradient);
    deallocate(detJ);
    deallocate(dV);
    deallocate(stress);
    deallocate(hessian);
    deallocate(traction);
    deallocate(stiffness);
    deallocate(geometric_stiffness);

    deallocate(GaussFieldVariables);
}


#endif // _LLADF__H
