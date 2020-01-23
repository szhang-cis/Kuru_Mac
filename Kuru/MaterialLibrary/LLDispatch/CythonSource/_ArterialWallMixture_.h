#include "_MaterialBase_.h"

template<typename U>
class _ArterialWallMixture_ : public _MaterialBase_<U> {
public:
    U mu;
    U kappa;
    U k1m;
    U k2m;
    U k1c;
    U k2c;

    FASTOR_INLINE _ArterialWallMixture_() = default;

    FASTOR_INLINE
    _ArterialWallMixture_(U mu, U kappa, U k1m, U k2m, U k1c, U k2c) {
        this->mu = mu;
        this->kappa = kappa;
	this->k1m = k1m;
	this->k2m = k2m;
	this->k1c = k1c;
	this->k2c = k2c;
    }

    FASTOR_INLINE
    void SetParameters(U mu, U kappa, U k1m, U k2m, U k1c, U k2c){
        this->mu = mu;
        this->kappa = kappa;
	this->k1m = k1m;
	this->k2m = k2m;
	this->k1c = k1c;
	this->k2c = k2c;
    }


    template<typename T=U, size_t ndim>
    FASTOR_INLINE
    std::tuple<Tensor<T,ndim,ndim>, typename MechanicsHessianType<T,ndim>::return_type>
    _KineticMeasures_(const T *Fnp, const T *Nnp, int nfibre, const T *GRnp, const T *Gnp)
    {

        // CREATE FASTOR TENSORS
	// deformation gradient
        Tensor<T,ndim,ndim> F;
	// local rotation tensor
        Tensor<T,ndim,ndim> R;
	Tensor<T,ndim> Normal;
	// elastin deposition stretch in cylindrcal coordinates
        Tensor<T,ndim,ndim> Ge;
        // COPY NUMPY ARRAY TO FASTOR TENSORS
        copy_numpy(F,Fnp);
        copy_numpy(R,Nnp);
	copy_numpy(Normal,Nnp);
        copy_numpy(Ge,Gnp);

	// Fibres deposition stretches
	T Gm = Gnp[9], Gc = Gnp[10];

        // FIND THE KINEMATIC MEASURES
        Tensor<Real,ndim,ndim> I; I.eye2();
        auto J = determinant(F);
	auto outerNormal = einsum<Index<i>,Index<j>>();

	// COMPUTE CAUCHY STRESS TENSOR FOR ELASTIN
	Ge = matmul(transpose(R),Ge);
	Ge = matmul(Ge,R);
	auto F_ela = matmul(F,Ge);
	auto J_ela = determinant(F_ela);
        auto b_ela = matmul(F_ela,transpose(F_ela));
	T trb = trace(b_ela);
	if (ndim==2) { trb += 1.; }
	T coeff0 = std::pow(J_ela,-2./3.);
        Tensor<T,ndim,ndim> stress = (mu*coeff0/J*(b_ela-1./3.*trb*I) + kappa*J_ela/J*(J-1.)*I)*GRnp[0];

        // FIND ELASTICITY TENSOR
        auto II_ijkl = einsum<Index<i,j>,Index<k,l>>(I,I);
        auto II_ikjl = permutation<Index<i,k,j,l>>(II_ijkl);
        auto II_iljk = permutation<Index<i,l,j,k>>(II_ijkl);

	auto Ib_ijkl = einsum<Index<i,j>,Index<k,l>>(I,b_ela);
	auto bI_ijkl = einsum<Index<i,j>,Index<k,l>>(b_ela,I);

        Tensor<T,ndim,ndim,ndim,ndim> elasticity = (2.*mu*coeff0/J*(1./9.*trb*II_ijkl-
			1./3.*(Ib_ijkl+bI_ijkl)+1./6.*trb*(II_ikjl+II_iljk))+
			kappa*J_ela/J*((2.*J-1.)*II_ijkl-(J-1.)*(II_ikjl+II_iljk)))*GRnp[0];

	// LOOP OVER FIBRES ORIENTATIONS
	for (int n=1; n<nfibre; ++n) {
           Tensor<T,ndim> N;
           copy_numpy(N,Nnp+3*n);
	   auto FN = matmul(F,N);
	   if (n==1) {
	      FN *= Gm;
	      auto innerFN = inner(FN,FN);
	      auto outerFN = outer(FN,FN);
	      auto FNFN_ijkl = einsum<Index<i,j>,Index<k,l>>(outerFN,outerFN);
	      // Passive component of muscles
	      T coeff1 = std::pow((innerFN-1.),2);
	      T coeff2 = std::exp(k2m*coeff1);
	      stress += 2.*k1m*GRnp[1]/J*(innerFN-1.)*coeff2*outerFN;
	      elasticity += 4.*k1m*GRnp[1]/J*(1.+2.*k2m*coeff1)*coeff2*FNFN_ijkl;
	      // Active component of muscles
	      T den0=1050., s_act=54000.*GRnp[1], stretch_m=1.4, stretch_a=1.,stretch_0=0.8;
	      T coeff3 = std::pow(((stretch_m-stretch_a)/(stretch_m-stretch_0)),2);
	      T coeff4 = std::pow(innerFN,2);
	      stress += (s_act/(den0*innerFN*J))*(1.-coeff3)*outerFN;
	      elasticity += -2.*(s_act/(den0*coeff4*J))*(1.-coeff3)*FNFN_ijkl;
	   }
	   else if (n>1) {
	      FN *= Gc;
	      auto innerFN = inner(FN,FN);
	      auto outerFN = outer(FN,FN);
	      auto FNFN_ijkl = einsum<Index<i,j>,Index<k,l>>(outerFN,outerFN);
	      // Component from the collagen
	      T coeff1 = std::pow((innerFN-1.),2);
	      T coeff2 = std::exp(k2c*coeff1);
	      stress += 2.*k1c*GRnp[n]/J*(innerFN-1.)*coeff2*outerFN;
	      elasticity += 4.*k1c*GRnp[n]/J*(1.+2.*k2c*coeff1)*coeff2*FNFN_ijkl;
	   }
	}

        auto hessian = voigt(elasticity);

        auto kinetics = std::make_tuple(stress,hessian);
        return kinetics;
    }



    template<typename T>
    void KineticMeasures(T *Snp, T* Hnp, int ndim, int ngauss, const T *Fnp, int nfibre,
		    const T *Nnp, const T *GRnp, const T *Gnp);

};

template<> template<>
void _ArterialWallMixture_<Real>::KineticMeasures<Real>(Real *Snp, Real* Hnp, int ndim, int ngauss, 
		const Real *Fnp, int nfibre, const Real *Nnp, const Real *GRnp, const Real *Gnp)
{

    if (ndim==3) {
        Tensor<Real,3> D;
        Tensor<Real,3,3> stress;
        Tensor<Real,6,6> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,3>(Fnp+9*g, Nnp, nfibre, GRnp, Gnp);
            copy_fastor(Snp,stress,g*9);
            copy_fastor(Hnp,hessian,g*36);
        }
    }
    else if (ndim==2) {
        Tensor<Real,2> D;
        Tensor<Real,2,2> stress;
        Tensor<Real,3,3> hessian;
        for (int g=0; g<ngauss; ++g) {
            std::tie(stress,hessian) =_KineticMeasures_<Real,2>(Fnp+4*g, Nnp, nfibre, GRnp, Gnp);
            copy_fastor(Snp,stress,g*4);
            copy_fastor(Hnp,hessian,g*9);
        }
    }
}
