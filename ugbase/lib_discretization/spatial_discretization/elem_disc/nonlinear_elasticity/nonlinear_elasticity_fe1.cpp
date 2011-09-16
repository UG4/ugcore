/*
 * fe1_nonlinear_elastictiy_fe1.cpp
 *
 *  Created on: 18.05.2011
 *      Author: raphaelprohl
 */

#include "fe1_nonlinear_elasticity.h"

#include "lib_discretization/spatial_discretization/disc_util/finite_element_geometry.h"
#include "common/util/provider.h"
#include "lib_discretization/local_finite_element/lagrange/lagrange.h"
#include "lib_discretization/quadrature/gauss_quad/gauss_quad.h"

#include "common/math/math_vector_matrix/math_matrix_functions.h"

////////////////////////////////////////////////////////////////////////////////
//	variables global to this file
////////////////////////////////////////////////////////////////////////////////

static number mu,kappa,K_0,K_inf,eta,Hard,omega,Exponent,poiss,emodul,MaxHardIter,HardAccur;


namespace ug{


////////////////////////////////////////////////////////////////////////////////
//	functions called at every ip
////////////////////////////////////////////////////////////////////////////////

static bool MaterialData()
{
	//properties for steel (in N/mm^2):

	emodul = 206900.00; //207e9 in N/m^2 = Pa
	poiss = 0.29;
	//lambda = //emodul*poiss/((1.0+poiss)*(1.0-2.0*poiss));
	//lambda = 110743.82;
	mu = emodul/(2.0 * (1.0 + poiss));
	//mu = 80193.80;
	kappa = 164206.00; //lambda + 2.0 * mu / 3.0;
	//for 4500bar=450e6 Pa steel begins to flow
	K_0 = 450.00; //Yields-stress (in N/mm^2)
	K_inf = 715.00;
	eta = 1.0;

	//set Hard = 0.0 & omega = 0.0 for perfect plasticity!
	//set Omega = 0.0 for linear hardening!
	Hard = 129.24; // (in N/mm^2); 129.24e6 Pa
	omega = 16.93;
	Exponent = 1.0;


	//accuracy-parameters
	MaxHardIter = 100; //max. steps in Newton-scheme for approximation of the nonlinear hardening law
	HardAccur = 1e-10; //for abs-truncation-criterion in nonlinear-hardening approximation (Newton-truncation)

	return true;
}

//"DeformationGradient": computes the deformationgradient at a ip, cf. Simo/Hughes 98 p.315
template<typename TDomain, typename TElem>
static bool
//DeformationGradient(size_t ip, const LocalVector& u, MathMatrix<TDomain::dim,TDomain::dim>& F)
DeformationGradient(size_t ip, const typename FE1NonlinearElasticityElemDisc<TDomain>::local_vector_type& u, MathMatrix<TDomain::dim,TDomain::dim>& F)
{
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	static const int dim = TDomain::dim;

	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
			= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	MatSet(F,0); //set F to Zero-Matrix

	//ip_loc = geo.local_ip(ip); //local-coordinates at ip

	for(size_t i = 0; i < geo.num_sh() ; ++i) // loop corners of element
	{
		MathVector<dim> Gradeps, Grad;
		Gradeps = geo.local_grad(ip,i); //local-gradient at ip
		Grad = geo.global_grad(ip,i); //global-gradient at ip

		for(size_t j = 0; j < (size_t)dim ; ++j)
		{
			for(size_t k = 0; k < (size_t)dim ; ++k)
			{
				F(j,k) += Grad[k] * u(i,j); //Grad(k)? u(i,j)? u(fct,dof)!
			}
		}

	}
	return false;
}

//"Hardening": nonlinear hardening law
static number Hardening(number alpha)
{
	return(K_0 + Hard * alpha + (K_inf - K_0) * (1.0 - exp(-omega * alpha)));
}

//"Hardening_d": derivative of nonlinear hardening law "Hardening"
static number Hardening_d(number alpha)
{
    return(Hard + (K_inf - K_0) * omega * exp(-omega * alpha));
}

//"SetCP": sets the plastic part of the symmetric right Cauchy_Green-strain tensor (3D-Version!)
template<typename TDomain>
static bool
SetCP(MathMatrix<TDomain::dim,TDomain::dim>& z, MathMatrix<TDomain::dim,TDomain::dim>& CPnew)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> IDENT;

	MatIdentity(IDENT);
	MatAdd(CPnew,z,IDENT); //CPnew = z + IDENT

	/*UG3-version:
	CPnew(0) = 1.0 + z(0);
	CPnew(1) = z(3);
	CPnew(2) = z(4);
	CPnew(3) = z(3);
	CPnew(4) = 1.0 + z(1);
	CPnew(5) = z(5);
	CPnew(6) = z(4);
	CPnew(7) = z(5);
	CPnew(8) = 1.0 + z(2);*/

	return false;
}

//"Inner": computes inner evolution of strains, hardening variables,...
//
//	input: dt: time increment
//			F: deformationgradient
//		CPold: old inner state parameters (quantity in convected config)
//
//	output: CPnew: new inner state parameters (quantity in convected config)
//			alpha: new inner hardening parameter
//
//	description: 	1. Transfer of CPold to spatial config by Push-Forward-Operation (trial elastic state)
//					2. Checking if yield-condition is violated in spatial config
//					3. if yield-condition is violated -> 	a) computing consistency-parameter gamma
//															b) Pull-Back in convected config
//															c) Update of CP (CPold to CPnew) by implicit euler step (return-mapping)
//					4. else CPnew = CPold

template<typename TDomain>
static bool
Inner(number dt, MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPold, MathMatrix<TDomain::dim,TDomain::dim>& CPnew, number& alpha)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> F1, FI, FT, FTI, CI, IDENT, help, btrial, devbtrial, mudevbtrial, normaltrial, Normaltrial, dCP;
	number detF, tracebtrial, snorm, flowcond, gamma, J2, J23;

	MatIdentity(IDENT);
	gamma = 0.0;

	MatAdd(F1,F,IDENT); //F1 = F + IDENT
	Transpose(FT,F1);
	Inverse(FI,F1);
	Transpose(FTI,FI);
	detF = Determinant(F1);
	J2 = detF * detF;
	J23 = 1.0 / pow(J2,1.0/dim);

	//TODO: anstatt MatMultiply TensorMultiply o.Ä. implementieren!
	//computing Inverse of right Cauchy-Green-Tensor
	MatMultiply(CI,FI,FTI);

	//Transfer of CPold to spatial config by Push-Forward-Operation (trial elastic state)
	SetCP(CPold,CPnew);

	//computing the finger tensor
	//btrial = F1 * CPnew * FT; //TODO: eventuell in math_matrix_function Multiplikation dreier Matrizen einfügen!
	MatMultiply(help,CPnew,FT);
	MatMultiply(btrial,F1,help);

	//compute btrialquer
	MatScale(btrial,J23,btrial);

	Trace(tracebtrial,btrial);

	//computing the deviatoric part of btrialquer
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			devbtrial(i,j) = btrial(i,j);
		}
		devbtrial(i,i) -= 1.0/3.0 * tracebtrial;
	}
	//TODO: Deviator(devbtrial,btrial);
	MatScale(mudevbtrial,mu,devbtrial);

	snorm = MatFrobeniusNorm(mudevbtrial);

	//checking yield-condition:
	flowcond = snorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	if (flowcond <= 0)
	{
		return true; //now: CPnew = CPold!
	}

	//return-mapping (corrector-step)

	number fdot, muquer, update;
	gamma = 0.0;
	muquer = mu * tracebtrial / 3.0;

	for(size_t i = 0; i < (size_t)MaxHardIter; ++i)
	{
		//fdot has to be equal to 0, due to the consistency condition
		fdot = snorm - 2.0 * gamma * muquer - sqrt(2.0/3.0) * Hardening(alpha + sqrt(2.0/3.0) * gamma);
		if (abs(fdot) < HardAccur * snorm)
		{
			//gamma found with "HardAccur * snorm" accuracy
			break;
		}
		gamma += fdot / (2.0 * muquer * (1.0 + Hardening_d(alpha + sqrt(2.0/3.0) * gamma) / (3.0 * muquer)));
	}

	MatScale(normaltrial,1.0/snorm,mudevbtrial);

	//computing Normaltrial for quantities in convected config
	MatMultiply(help,normaltrial,FTI);
	MatMultiply(Normaltrial,FI,help);

	//flow-rule modification:
	update = - 2.0 / 3.0 * gamma * tracebtrial;
	MatScale(dCP,update,Normaltrial);

	//CP-Update
	CPnew(0) = CPold(0) + dCP(0); //TODO: voller CP-Tensor?
	CPnew(1) = CPold(1) + dCP(1);
	CPnew(2) = CPold(2) + dCP(2);
	CPnew(3) = CPold(3) + dCP(3);
	CPnew(4) = CPold(4) + dCP(4);
	CPnew(5) = CPold(5) + dCP(5);
	CPnew(6) = CPold(6) + dCP(6);
	CPnew(7) = CPold(7) + dCP(7);
	CPnew(8) = CPold(8) + dCP(8);

	//Hardening-parameter Update
	alpha = alpha + sqrt(2.0/3.0) * gamma;

	//TODO: hier elementweisen data-Pointer füllen!
	/*UG3-Version:
	CPnew(0) = CPold(0) + dCP[0];
	CPnew(1) = CPold(1) + dCP[4];
	CPnew(2) = CPold(2) + dCP[8];
	CPnew(3) = CPold(3) + dCP[1];
	CPnew(4) = CPold(4) + dCP[2];
	CPnew(5) = CPold(5) + dCP(5);*/

	return false;
}

//"W_F_<..>": computes stress tensors after inner evolution
//
//	input: 	F: deformationgradient
//		CPnew: new inner state parameters (quantity in convected config)
//
//	output: tau: Kirchhoff stress tensor
//			PK1: first Piola-Kirchhoff stress tensor
//			PK2: second Piola-Kirchhoff stress tensor
//
//	description: 		1. CPnew is: either CPold or updated by return-mapping
//						2. Push-Forward in spatial config
//						3. Addition of elastic stress

template<typename TDomain>
static bool
W_F_tau(MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPnew, MathMatrix<TDomain::dim,TDomain::dim>& tau)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> F1, FI, FT, FTI, IDENT, help, b, devb, s, CP;
	number detF, J2, J23, traceb;

	MatIdentity(IDENT);

	MatAdd(F1,F,IDENT); //F1 = F + IDENT
	Transpose(FT,F1);
	Inverse(FI,F1);
	Transpose(FTI,FI);
	detF = Determinant(F1);
	J2 = detF * detF;
	J23 = 1.0 / pow(J2,1.0/dim);

	SetCP(CPnew,CP);

	//computing finger-tensor b in spatial config
	MatMultiply(help,CP,FT);
	MatMultiply(b,F1,help);

	//compute bquer
	MatScale(b,J23,b);

	Trace(traceb,b);

	//computing the deviatoric part of b
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			devb(i,j) = b(i,j);
		}
		devb(i,i) -= 1.0/3.0 * traceb;
	}

	MatScale(s,mu,devb);

	//computing the Kirchhoff stress tensor tau
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			tau(i,j) = s(i,j);
		}
		tau(i,i) += 0.5 * kappa * (J2 - 1.0);
	}

	return false;
}

template<typename TDomain>
static bool
W_F_PK1(MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPnew, MathMatrix<TDomain::dim,TDomain::dim>& PK1)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> F1, FI, FT, FTI, IDENT, help, b, devb, s, CP, tau;
	number detF, J2, J23, traceb;

	MatIdentity(IDENT);

	MatAdd(F1,F,IDENT); //F1 = F + IDENT
	Transpose(FT,F1);
	Inverse(FI,F1);
	Transpose(FTI,FI);
	detF = Determinant(F1);
	J2 = detF * detF;
	J23 = 1.0 / pow(J2,1.0/dim);

	SetCP(CPnew,CP);

	//computing finger-tensor b in spatial config
	MatMultiply(help,CP,FT);
	MatMultiply(b,F1,help);

	//compute bquer
	MatScale(b,J23,b);

	Trace(traceb,b);

	//computing the deviatoric part of b
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			devb(i,j) = b(i,j);
		}
		devb(i,i) -= 1.0/3.0 * traceb;
	}

	MatScale(s,mu,devb);

	//computing the Kirchhoff stress tensor tau
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			tau(i,j) = s(i,j);
		}
		tau(i,i) += 0.5 * kappa * (J2 - 1.0);
	}

	//computing first Piola Kirchhoff stress tensor PK1
	MatMultiply(PK1,tau,FTI); //PK1 = tau * F^-T

	return false;
}

template<typename TDomain>
static bool
W_F_PK2(MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPnew, MathMatrix<TDomain::dim,TDomain::dim>& PK2)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> F1, FI, FT, FTI, IDENT, help, b, devb, s, CP, tau, PK1;
	number detF, J2, J23, traceb;

	MatIdentity(IDENT);

	MatAdd(F1,F,IDENT); //F1 = F + IDENT
	Transpose(FT,F1);
	Inverse(FI,F1);
	Transpose(FTI,FI);
	detF = Determinant(F1);
	J2 = detF * detF;
	J23 = 1.0 / pow(J2,1.0/dim);

	SetCP(CPnew,CP);

	//computing finger-tensor b in spatial config
	MatMultiply(help,CP,FT);
	MatMultiply(b,F1,help);

	//compute bquer
	MatScale(b,J23,b);

	Trace(traceb,b);

	//computing the deviatoric part of b
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			devb(i,j) = b(i,j);
		}
		devb(i,i) -= 1.0/3.0 * traceb;
	}

	MatScale(s,mu,devb);

	//computing the Kirchhoff stress tensor tau
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			tau(i,j) = s(i,j);
		}
		tau(i,i) += 0.5 * kappa * (J2 - 1.0);
	}

	//computing second Piola Kirchhoff stress tensor PK2
	MatMultiply(PK1,tau,FTI); //PK1 = tau * F^-T
	MatMultiply(PK2,FI,PK1); //PK2 = F^-1 * PK1

	return false;
}

//"Stress_<..>": computes stress response at a ip
//
//	input: dt: time increment
//			F: deformationgradient
//		CPold: old inner state parameters (quantity in convected config)
//
//	output: tau: Kirchhoff stress tensor
//			PK1: first Piola Kirchhoff stress tensor
//			PK2: second Piola Kirchhoff stress tensor

template<typename TDomain>
static bool
Stress_tau(number dt, MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPold, MathMatrix<TDomain::dim,TDomain::dim>& tau)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim>& CPnew; //CPnew[2*dim+1];
	number alpha; //TODO: alpha mit CP vereinen für data-Pointer-Struktur!

	//compute inner evolution of strains, hardening variables,...
	Inner(dt,F,CPold,CPnew,alpha);
	//compute stored energy functional
	W_F_tau(F,CPnew,tau);

	return false;
}

template<typename TDomain>
static bool
Stress_PK1(number dt, MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPold, MathMatrix<TDomain::dim,TDomain::dim>& PK1)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim>& CPnew; //CPnew[2*dim+1];
	number alpha;

	//compute inner evolution of strains, hardening variables,...
	Inner(dt,F,CPold,CPnew,alpha);
	//compute stored energy functional
	W_F_PK1(F,CPnew,PK1);

	return false;
}

template<typename TDomain>
static bool
Stress_PK2(number dt, MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPold, MathMatrix<TDomain::dim,TDomain::dim>& PK2)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim>& CPnew; //CPnew[2*dim+1];
	number alpha;

	//compute inner evolution of strains, hardening variables,...
	Inner(dt,F,CPold,CPnew,alpha);
	//compute stored energy functional
	W_F_PK2(F,CPnew,PK2);

	return false;
}

template<typename TDomain>
static bool
TangentInnerParameter(number dt, MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPold, MathMatrix<TDomain::dim,TDomain::dim>& CPnew, number& alpha, number& muquer, number& snorm, MathMatrix<TDomain::dim,TDomain::dim>& normal, MathMatrix<TDomain::dim,TDomain::dim>& devnorm2, number& gamma, MathMatrix<TDomain::dim,TDomain::dim>& s)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> F1, FI, FT, FTI, CI, IDENT, help, btrial, devbtrial, mudevbtrial, normaltrial, Normaltrial, dCP;
	number detF, J2, J23, tracebtrial, flowcond;

	MatIdentity(IDENT);
	gamma = 0.0;

	MatAdd(F1,F,IDENT); //F1 = F + IDENT
	Transpose(FT,F1);
	Inverse(FI,F1);
	Transpose(FTI,FI);
	detF = Determinant(F1);
	J2 = detF * detF;
	J23 = 1.0/ pow(J2,1.0/dim);

	//computing Inverse of right Cauchy-Green-Tensor
	MatMultiply(CI,FI,FTI);

	//Transfer of CPold to spatial config by Push-Forward-Operation (trial elastic state)
	SetCP(CPold,CPnew);

	//computing the finger tensor
	//btrial = F1 * CPnew * FT;
	MatMultiply(help,CPnew,FT);
	MatMultiply(btrial,F1,help);

	//compute btrialquer
	MatScale(btrial,J23,btrial);

	Trace(tracebtrial,btrial);

	//computing the deviatoric part of btrial
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			devbtrial(i,j) = btrial(i,j);
		}
		devbtrial(i,i) -= 1.0/3.0 * tracebtrial;
	}

	MatScale(mudevbtrial,mu,devbtrial);

	snorm = MatFrobeniusNorm(mudevbtrial);

	//checking yield-condition:
	flowcond = snorm - sqrt(2.0 / 3.0) * Hardening(alpha);

	if (flowcond <= 0)
	{
		return true; //now: CPnew = CPold!
	}

	//return-mapping (corrector-step)

	number fdot, update;
	gamma = 0.0;
	muquer = mu * tracebtrial / 3.0;

	for(size_t i = 0; i < (size_t)MaxHardIter; ++i)
	{
		//fdot has to be equal to 0, due to the consistency condition
		fdot = snorm - 2.0 * gamma * muquer - sqrt(2.0/3.0) * Hardening(alpha + sqrt(2.0/3.0) * gamma);
		if (abs(fdot) < HardAccur * snorm)
		{
			//gamma found with "HardAccur * snorm" accuracy
			break;
		}
		gamma += fdot / (2.0 * muquer * (1.0 + Hardening_d(alpha + sqrt(2.0/3.0) * gamma) / (3.0 * muquer)));
	}

	MatScale(normaltrial,1.0/snorm,mudevbtrial);

	//computing Normaltrial for quantities in convected config
	MatMultiply(help,normaltrial,FTI);
	MatMultiply(Normaltrial,FI,help);

	//flow-rule modification:
	update = - 2.0 / 3.0 * gamma * tracebtrial;
	MatScale(dCP,update,Normaltrial);

	//CP-Update
	CPnew(0) = CPold(0) + dCP(0);
	CPnew(1) = CPold(1) + dCP(1);
	CPnew(2) = CPold(2) + dCP(2);
	CPnew(3) = CPold(3) + dCP(3);
	CPnew(4) = CPold(4) + dCP(4);
	CPnew(5) = CPold(5) + dCP(5);
	CPnew(6) = CPold(6) + dCP(6);
	CPnew(7) = CPold(7) + dCP(7);
	CPnew(8) = CPold(8) + dCP(8);

	//Hardening-parameter Update
	alpha = alpha + sqrt(2.0/3.0) * gamma;

	/*UG3-Version:
	CPnew(0) = CPold(0) + dCP[0];
	CPnew(1) = CPold(1) + dCP[4];
	CPnew(2) = CPold(2) + dCP[8];
	CPnew(3) = CPold(3) + dCP[1];
	CPnew(4) = CPold(4) + dCP[2];
	CPnew(5) = CPold(5) + dCP(5);*/

	return false;
}

//"TangentExact": computes Elasticity-tensor as the derivative of the stress tensor r.t. strain
// derived by exact linearization of the return-mapping-algorithm, cf. Simo/Hughes p. 321
// elasticity tensor C is defined in the spatial config
//
//	input: dt: time increment
//			F: deformationgradient
//		CPold: old inner parameter
//
//	output: C: elasticity tensor

template<typename TDomain>
static bool
TangentExact(number dt, MathMatrix<TDomain::dim,TDomain::dim>& F, MathMatrix<TDomain::dim,TDomain::dim>& CPold, number& alpha, MathTensor<4,TDomain::dim>& C)
{
	static const int dim = TDomain::dim;

	MathMatrix<dim,dim> F1, FI, FT, FTI, IDENT, normal, devnorm2, s, CPnew;
	number detF, J2, muquer, snorm, gamma, beta0, beta1, beta2, beta3, beta4;

	MatIdentity(IDENT);

	MatAdd(F1,F,IDENT); //F1 = F + IDENT
	Transpose(FT,F1);
	Inverse(FI,F1);
	Transpose(FTI,FI);
	detF = Determinant(F1);
	J2 = detF * detF;

	TangentInnerParameter(dt,F,CPold,CPnew,alpha,muquer,snorm,normal,devnorm2,gamma,s);

	beta0 = 1.0 + Hardening_d(alpha + sqrt(2.0/3.0) * gamma) / (3.0 * muquer);
	beta1 = 2.0 * muquer * gamma / snorm;
	beta2 = (1.0 - 1.0/beta0) * 2.0/3.0 * snorm/muquer * gamma;
	beta3 = 1.0/beta0 - beta1 + beta2;
	beta4 = (1.0/beta0 - beta1) * snorm/muquer;

	//filling the elasticity tensor by exact linearization of the return-mapping-algorithm in Simo/Hughes p. 321
	for(size_t i = 0; i < (size_t)dim; ++i) // loop dimension
	{
		for(size_t j = 0; j < (size_t)dim; ++j) // loop dimension
		{
			for(size_t k = 0; k < (size_t)dim; ++k) // loop dimension
			{
				for(size_t l = 0; l < (size_t)dim; ++l) // loop dimension
				{
					C(i,j,k,l) = 0.0;

					if ((i==j) && (k==l)) //1x1-term
					{
						C(i,j,k,l) += kappa * J2;
						C(i,j,k,l) -= 2.0 * muquer / 3.0;
						//beta1-term
						C(i,j,k,l) += beta1 * 2.0 * muquer / 3.0;
					}

					if ((i==k) && (j==l)) //first summand of I
					{
						C(i,j,k,l) -= 0.5 * kappa * (J2 - 1.0);
						C(i,j,k,l) += muquer;
						//beta1-term
						C(i,j,k,l) -= beta1 * muquer;
					}

					if ((i==l) && (j==k)) //second summand of I
					{
						C(i,j,k,l) -= 0.5 * kappa * (J2 - 1.0);
						C(i,j,k,l) += muquer;
						//beta1-term
						C(i,j,k,l) -= beta1 * muquer;
					}

					if (k==l) //nx1-term
					{
						C(i,j,k,l) -= 2.0 * snorm * normal(i,j) / 3.0;
						//beta1-term
						C(i,j,k,l) += beta1 * 2.0 * snorm * normal(i,j) / 3.0;
					}

					if (i==j) //1xn-term
					{
						C(i,j,k,l) -= 2.0 * snorm * normal(k,l) / 3.0;
						//beta1-term
						C(i,j,k,l) += beta1 * 2.0 * snorm * normal(k,l) / 3.0;
					}

					//beta3-term
					C(i,j,k,l) -= 2.0 * muquer * beta3 * normal(i,j) * normal(k,l);

					//beta4-term
					C(i,j,k,l) -= muquer * beta4 * normal(i,j) * devnorm2(k,l);
					C(i,j,k,l) -= muquer * beta4 * devnorm2(i,j) * normal(k,l);
				}
			}
		}
	}

	return false;
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
template<typename TElem>
bool
FE1NonlinearElasticityElemDisc<TDomain>::
prepare_element_loop()
{
	// all this will be performed outside of the loop over the elements.
	// Therefore it is not time critical.

//	set material data
	MaterialData();
// 	evaluate Elasticity Tensor
	m_ElasticityTensorFct(m_ElasticityTensor);
// 	evaluate Stress Tensor
	m_StressTensorFct(m_StressTensor);

	return true;
}

template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
prepare_element(TElem* elem, const local_vector_type& u)
{
//	get corners
	m_corners = this->template get_element_corners<TElem>(elem);

	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

// 	update Geometry for this element
	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
			= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	if(!geo.update(elem, m_corners, LFEID(LFEID::LAGRANGE,1), 2))
		{UG_LOG("FE1NonlinearElasticityElemDisc::prepare_element:"
				" Cannot update Finite Element Geometry.\n"); return false;}

	return true;
}

//assemble stiffness matrix
template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
			= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	//using Lagrange description of the linearized equations (cf. Bonet/Wood 1997 chapter 8/9)

	number detF;
	MathMatrix<dim,dim> F, F1, FI, FT, FTI, CI, IDENT, DEa, DEb, grad; //FI(i,j);
	//TODO: Schleife für ip rausziehen???

	MatIdentity(IDENT);

	for(size_t a = 0; a < geo.num_sh(); ++a) // loop corner
	{
		for(size_t i = 0; i < num_fct(); ++i) // loop component
		{
			for(size_t b = 0; b < geo.num_sh(); ++b) // loop corner
			{
				for(size_t j = 0; j < num_fct(); ++j) // loop component
				{
					// Compute entry A_{a, i, b, j}
					number integrand = 0;
					number integrandC = 0;
					number integrandS = 0;
					for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
					{
						//computing deformationgradient:
						DeformationGradient<TDomain,TElem>(ip,u,F);//mm: Anzahl der Komponenten des SolutionVectors auf Elementebene
						//-> F

						MatAdd(F1,F,IDENT); //F1 = F + IDENT
						Transpose(FT,F1);
						Inverse(FI,F1);
						Transpose(FTI,FI);
						detF = Determinant(F1);

						//computing Inverse of right Cauchy-Green-Tensor
						MatMultiply(CI,FI,FTI);

						//computing 2.PK-stress tensor: Stress_2PK-Aufruf!
						//-> T

						//computing elasticityMatrix: Funktionsaufruf für TangentExactSpatial
						//-> C

						/* for formulation in convected config:
						for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
						{
							for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
							{
								number Gsyma = 0.0;
								number Gsymb = 0.0;
								DEa(d1,d2) = 0.0;
								DEb(d1,d2) = 0.0;

								Gsyma += 0.5 * (FT(d1,i) * geo.global_grad(ip, a)[d2] + FT(d2,i) * geo.global_grad(ip, a)[d1]);
								DEa(d1,d2) += Gsyma;
								Gsymb += 0.5 * (FT(d1,j) * geo.global_grad(ip, b)[d2] + FT(d2,j) * geo.global_grad(ip, b)[d1]);
								DEb(d1,d2) += Gsymb;
							}
						}*/

						for(size_t k = 0; k < (size_t)dim; ++k) // loop dimension
						{
							for(size_t l = 0; l < (size_t)dim; ++l) // loop dimension
							{
								//computing the material contribution to the stiffness matrix

								/*for formulation in convected config:
								for(size_t C1 = 0; C1 < num_fct(); ++C1) // loop component
								{
									for(size_t C2 = 0; C2 < num_fct(); ++C2) // loop component
									{
										integrandC += DEa(C1,k) * m_ElasticityTensor[C1][k][C2][l] * DEb(C2,l); //geo.global_grad(ip, a)[k] * m_ElasticityTensor[c1][k][c2][l] * geo.global_grad(ip, b)[l];
									}
								}*/

								//for formulation in spatial config compute the "spatial grad"
								number grada = 0.0;
								number gradb = 0.0;

								for(size_t K = 0; K < (size_t)dim; ++K) // loop dimension
								{
									grada += geo.global_grad(ip, a)[K] * FI(K,k);
									gradb += geo.global_grad(ip, b)[K] * FI(K,l);
								}
								grad(a,k) = grada;
								grad(b,l) = gradb;

								integrandC += grad(a,k) * m_ElasticityTensor[i][k][j][l] * grad(b,l);

								//computing the geometrical contribution to the stiffness matrix
								integrandS += geo.global_grad(ip, a)[k] * m_StressTensor[k][l] * geo.global_grad(ip, b)[l];
							}
						}
						integrand = geo.weight(ip) * (integrandC + integrandS);
					}

					J(i, a, j, b) += integrand; //richtige Reihenfolge der Komponenten?
				}
			}
		}
	}

	return true;
}


//assemble mass matrix
template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
			= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
				// same value for all components
				number value = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);
				for(size_t c = 0; c < num_fct(); ++c)
				{
					J(c, i, c, j) += value;
				}
			}
		}
	}

	return true;
}


//assemble defect
template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

	// to be implemented
	FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> >& geo
			= Provider<FEGeometry<TElem, dim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 2> > >::get();

	number Vol, vol, detF;
	MathMatrix<dim,dim> F1, FI, F, FT, FTI, CI, IDENT, DE;

	MatIdentity(IDENT);

	//TODO: mean Dilatation Term fehlt noch!
	//TODO: Schleife für ip rausziehen???

	Vol = vol = 0.0;

	//computing element-volume in current and reference configuration
	/*for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
	{
		number FT[dim][dim], F[dim][dim];
		number detF;
		detF = 1.0;

		//DeformationGradient(ip,u,F);
		//computing deformationgradient F and det(F)
		for(size_t I = 0; I < (size_t)dim; ++I)
		{
			for(size_t j = 0; j < (size_t)dim; ++j)
			{
				FT[I][j] = F[j][I];
			}
			FT[I][I] += 1.0;
			detF *= FT[I][I];
		}
		vol += geo.weight(ip) * detF;
		Vol += geo.weight(ip);
	}*/

	for(size_t i = 0; i < geo.num_sh(); ++i) // loop corner
	{
		for(size_t c = 0; c < num_fct(); ++c) // loop component
		{
			number integrand = 0;
			for(size_t ip = 0; ip < geo.num_ip(); ++ip) // loop ip
			{
				//computing deformationgradient:
				DeformationGradient<TDomain,TElem>(ip,u,F);//mm: Anzahl der Komponenten des SolutionVectors auf Elementebene
				//-> F

				MatAdd(F1,F,IDENT); //F1 = F + IDENT
				Transpose(FT,F1);
				Inverse(FI,F1);
				Transpose(FTI,FI);
				detF = Determinant(F1);

				//Aufruf von Stress_1PK(ABS(time-time_old),F,data,T); //muss m_StressTensor füllen!

				for(size_t d1 = 0; d1 < (size_t)dim; ++d1) // loop dimension
				{
					integrand += geo.global_grad(ip, i)[d1] * m_StressTensor[c][d1];
					/*for(size_t d2 = 0; d2 < (size_t)dim; ++d2) // loop dimension
					{
						number Gsym = 0.0;

						DE(d1,d2) = 0.0;
						Gsym += 0.5 * (FT(d1,c) * geo.global_grad(ip, i)[d2] + FT(d2,c) * geo.global_grad(ip, i)[d1]);
						DE(d1,d2) += Gsym;

						integrand += DE(d1,d2) * m_StressTensor[d1][d2];
					}*/

				}
				integrand *= geo.weight(ip); //in weight(ip) steckt schon die Multiplikation mit Jdet!
			}

			d(i,c) -= integrand; //+= oder -=?
		}
	}
	return true;
}


template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// Not implemented
	return false;
}


template<typename TDomain>
template<typename TElem>
inline
bool
FE1NonlinearElasticityElemDisc<TDomain>::
assemble_f(local_vector_type& d)
{
	// Not implemented
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FE1NonlinearElasticityElemDisc<TDomain>::
FE1NonlinearElasticityElemDisc()
{
	register_all_fe1_funcs();
};


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FE1NonlinearElasticityElemDisc<TDomain>::
register_all_fe1_funcs()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>( RegisterFE1(this) );
}

template<typename TDomain>
template<typename TElem>
void
FE1NonlinearElasticityElemDisc<TDomain>::
register_fe1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem>);
	set_prep_elem_fct(	 id, &T::template prepare_element<TElem>);
	set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem>);
	set_ass_JA_elem_fct(		 id, &T::template assemble_JA<TElem>);
	set_ass_JM_elem_fct(		 id, &T::template assemble_JM<TElem>);
	set_ass_dA_elem_fct(		 id, &T::template assemble_A<TElem>);
	set_ass_dM_elem_fct(		 id, &T::template assemble_M<TElem>);
	set_ass_rhs_elem_fct(	 id, &T::template assemble_f<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FE1NonlinearElasticityElemDisc<Domain1d>;
template class FE1NonlinearElasticityElemDisc<Domain2d>;
template class FE1NonlinearElasticityElemDisc<Domain3d>;

} // namespace ug
