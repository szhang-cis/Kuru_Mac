#include "_MaterialBase_.h"

template<typename U>
class _ArterialWallMixture_ : public _MaterialBase_<U> {
public:
    U mu;
    U kappa;
    U k1tm;
    U k2tm;
    U k1tc;
    U k2tc;
    U k1cm;
    U k2cm;
    U k1cc;
    U k2cc;
    U rho;
    U s_act;
    U stretch_m;
    U stretch_0;
    U stretch_a;
    U pressure;

    FASTOR_INLINE _ArterialWallMixture_() = default;

    FASTOR_INLINE
    _ArterialWallMixture_(U mu,U kappa,U k1tm,U k2tm,U k1tc,U k2tc,U k1cm,U k2cm,U k1cc,U k2cc,U rho,U s_act,U stretch_m,U stretch_0,U stretch_a,U pressure) {
    this->mu = mu;
    this->kappa = kappa;
    this->k1tm = k1tm;
	this->k2tm = k2tm;
	this->k1tc = k1tc;
    this->k2tc = k2tc;
	this->k1cm = k1cm;
    this->k2cm = k2cm;
	this->k1cc = k1cc;
    this->k2cc = k2cc;
    this->rho = rho;
    this->s_act = s_act;
    this->stretch_m = stretch_m;
    this->stretch_0 = stretch_0;
    this->stretch_a = stretch_a;
    this->pressure = pressure;
    }

    FASTOR_INLINE
    void SetParameters(U mu,U kappa,U k1tm,U k2tm,U k1tc,U k2tc,U k1cm,U k2cm,U k1cc,U k2cc,U rho,U s_act,U stretch_m,U stretch_0,U stretch_a,U pressure){
    this->mu = mu;
    this->kappa = kappa;
	this->k1tm = k1tm;
    this->k2tm = k2tm;
    this->k1tc = k1tc;
    this->k2tc = k2tc;
    this->k1cm = k1cm;
    this->k2cm = k2cm;
    this->k1cc = k1cc;
    this->k2cc = k2cc;
    this->rho = rho;
    this->s_act = s_act;
    this->stretch_m = stretch_m;
    this->stretch_0 = stretch_0;
    this->stretch_a = stretch_a;
    this->pressure = pressure;
    }


    template<typename T=U, size_t ndim>
    FASTOR_INLINE
    std::tuple<Tensor<T,ndim,ndim>, typename MechanicsHessianType<T,ndim>::return_type>
    _KineticMeasures_(const T *Fnp, const T *Nnp, int nfibre, const T *SVnp, int near_incomp)
    {

    // CREATE FASTOR TENSORS
	// deformation gradient
    Tensor<T,ndim,ndim> F;
    // local rotation tensor
    Tensor<T,ndim,ndim> R;
    Tensor<T,ndim> Normal;
    // elastin deposition stretch in cylindrcal coordinates
    Tensor<T,ndim,ndim> F_r;
    // COPY NUMPY ARRAY TO FASTOR TENSORS
    copy_numpy(F,Fnp);
    copy_numpy(R,Nnp);
	copy_numpy(Normal,Nnp);
    copy_numpy(F_r,SVnp);

    // FIND THE KINEMATIC MEASURES
    Tensor<Real,ndim,ndim> I; I.eye2();
    auto J = determinant(F);

	// Growth tensor definition
	auto outerNormal = outer(Normal,Normal);
	auto outerTangential = I - outerNormal;
	Tensor<T,ndim,ndim> F_g = SVnp[20]*outerNormal + outerTangential;
	// Remodeling tensor for elastin
	F_r = matmul(transpose(R),F_r);
	F_r = matmul(F_r,R);
	// Elastin Growth and Remodeling
	Tensor<T,ndim,ndim> F_gr = matmul(F_r,F_g);
	Tensor<T,ndim,ndim> F_gr_inv = inverse(F_gr);

	// COMPUTE CAUCHY STRESS TENSOR FOR ELASTIN
	auto F_e = matmul(F,F_gr_inv);
	auto J_e = determinant(F_e);
    auto b_e = matmul(F_e,transpose(F_e));
	T trb = trace(b_e);
	if (ndim==2) { trb += 1.; }
	T coeff0 = std::pow(J_e,-2./3.);

    // ISOCHORIC PART OF CAUCHY STRESS
    Tensor<T,ndim,ndim> stress = mu*coeff0*(b_e-1./3.*trb*I)*SVnp[14]/J;

    // FIND ELASTICITY TENSOR FOR ELASTIN
    auto II_ijkl = einsum<Index<i,j>,Index<k,l>>(I,I);
    auto II_ikjl = permutation<Index<i,k,j,l>>(II_ijkl);
    auto II_iljk = permutation<Index<i,l,j,k>>(II_ijkl);

	auto Ib_ijkl = einsum<Index<i,j>,Index<k,l>>(I,b_e);
	auto bI_ijkl = einsum<Index<i,j>,Index<k,l>>(b_e,I);

    // ISOCHORIC PART OF ELASTICITY TENSOR
    Tensor<T,ndim,ndim,ndim,ndim> elasticity = 2.*mu*coeff0*(1./9.*trb*II_ijkl-
		1./3.*(Ib_ijkl+bI_ijkl)+1./6.*trb*(II_ikjl+II_iljk))*SVnp[14]/J;

    // VOLUMETRIC PARTO OF ELASTICITY TENSOR
    if (near_incomp) {
       stress += J_e*pressure*I*SVnp[14]/J;
       elasticity += pressure*J_e*(II_ijkl-(II_ikjl+II_iljk))*SVnp[14]/J;
    }
    else {
       stress += kappa*J_e*(J_e-1.)*I*SVnp[14]/J;
       elasticity += kappa*J_e*((2.*J_e-1.)*II_ijkl-(J_e-1.)*(II_ikjl+II_iljk))*SVnp[14]/J;
    }

	// LOOP OVER FIBRES ORIENTATIONS
	for (int n=1; n<nfibre; ++n) {
       Tensor<T,ndim> N;
       copy_numpy(N,Nnp+3*n);
	   auto FN = matmul(F,N);
	   // Remodeling stretch along the fibre
	   T lambda_r = SVnp[8+n];
	   // TOTAL deformation
	   auto innerFN = inner(FN,FN);
	   auto outerFN = outer(FN,FN);
	   // ELASTIC deformation
	   T coeff1 = std::pow(lambda_r,2);
	   auto innerFN_e = innerFN/coeff1;
	   auto outerFN_e = outerFN/coeff1;
	   T coeff2 = std::pow((innerFN_e-1.),2);
	   auto FNFN_ijkl_e = einsum<Index<i,j>,Index<k,l>>(outerFN_e,outerFN_e);
       T k1=0.; T k2=0.;
	   if (n==1) {
          if (innerFN_e>=1.) {k1=k1tm; k2=k2tm;} else {k1=k1cm; k2=k2cm;}
	      // Passive component of muscles
	      T coeff3 = std::exp(k2*coeff2);
	      stress += 2.*k1*SVnp[15]/J*(innerFN_e-1.)*coeff3*outerFN_e;
	      elasticity += 4.*k1*SVnp[15]/J*(1.+2.*k2*coeff2)*coeff3*FNFN_ijkl_e;
	      // Active component of muscles
	      auto FNFN_ijkl = einsum<Index<i,j>,Index<k,l>>(outerFN,outerFN);
	      T coeff4 = std::pow(((stretch_m-stretch_a)/(stretch_m-stretch_0)),2);
	      T coeff5 = std::pow(innerFN,2);
	      stress += (s_act*SVnp[15]/(rho*innerFN*J))*(1.-coeff4)*outerFN;
	      elasticity += -2.*(s_act*SVnp[15]/(rho*coeff5*J))*(1.-coeff4)*FNFN_ijkl;
	   }
	   else if (n>1) {
          if (innerFN_e>=1.) {k1=k1tc; k2=k2tc;} else {k1=k1cc; k2=k2cc;}
	      // Component from the collagen
	      T coeff3 = std::exp(k2*coeff2);
	      stress += 2.*k1*SVnp[14+n]/J*(innerFN_e-1.)*coeff3*outerFN_e;
	      elasticity += 4.*k1*SVnp[14+n]/J*(1.+2.*k2*coeff2)*coeff3*FNFN_ijkl_e;
	    }
	}

    auto hessian = voigt(elasticity);

    auto kinetics = std::make_tuple(stress,hessian);
    return kinetics;
    }


    template<typename T>
    void KineticMeasures(T *Snp, T* Hnp, int ndim, int ngauss, const T *Fnp, int nfibre,
		    const T *Nnp, int nstatv, const T *SVnp, int near_incomp);

};

template<> template<>
void _ArterialWallMixture_<Real>::KineticMeasures<Real>(Real *Snp, Real* Hnp, int ndim, int ngauss, 
		const Real *Fnp, int nfibre, const Real *Nnp, int nstatv, const Real *SVnp, int near_incomp)
{

    if (ndim==3) {
        Tensor<Real,3,3> stress;
        Tensor<Real,6,6> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,3>(Fnp+9*g, Nnp, nfibre, SVnp+nstatv*g, near_incomp);
            copy_fastor(Snp,stress,g*9);
            copy_fastor(Hnp,hessian,g*36);
        }
    }
    else if (ndim==2) {
        Tensor<Real,2,2> stress;
        Tensor<Real,3,3> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,2>(Fnp+4*g, Nnp, nfibre, SVnp+nstatv*g, near_incomp);
            copy_fastor(Snp,stress,g*4);
            copy_fastor(Hnp,hessian,g*9);
        }
    }
}

