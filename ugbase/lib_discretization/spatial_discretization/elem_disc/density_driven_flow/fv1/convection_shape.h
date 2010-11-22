/*
 * convection_shape.h
 *
 *  Created on: 29.10.2010
 *      Author: andreasvogel
 *
 *
 *  Note: 	This file is adapted from 'upwind.c' in the d3f-package.
 *  		The original implementation is by Dr. D. Logashenko
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONVECTION_SHAPE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONVECTION_SHAPE__

// other ug4 modules
#include "common/common.h"
#include "common/math/ugmath.h"

namespace ug{

/**
 * \param[in]		fvg				Finite Volume Geometry
 * \param[in]		scvf			SubControlVolumeFace to consider (ip)
 * \param[in]		vel				Convection Velocity at ip
 * \param[in]		Diffusion		Diffusion Tensor at ip
 * \param[out]		shape			Convection shape for each dof
 * \param[out]		D_shape_vel		Derivative of the shape w.r.t. the Velocity for each dof
 * \param[out]		D_shape_D		Derivative of the shape w.r.t. the Diffusion for each dof
 * \param[out]		bNonZeroDerivD	nonzero iff the deriv. w.r.t. the diffusion tensor is non-zero
 *
 */
template <typename TFVGeometry>
bool ConvectionShapes_NoUpwind
(
	const TFVGeometry& fvg,
    const typename TFVGeometry::SCVF& scvf,
    const MathVector<TFVGeometry::world_dim>& Vel,
    const MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>& Diffusion,
    number* shape,
    MathVector<TFVGeometry::world_dim>* D_shape_vel,
    MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>* D_shape_D,
    bool* bNonZeroDerivD
)
{
//	Currently only first order
	static const size_t ip = 0;
	UG_ASSERT(scvf.num_ip() == 1, "Only first order implemented.");

//	Compute flux
    const number flux = VecDot(scvf.normal(), Vel);

//	Write Shapes
    for(size_t sh = 0; sh < scvf.num_sh(); sh++)
    {
        shape[sh] = flux * scvf.shape(sh, ip);
    }

//	Write Derivatives if wanted
    if (D_shape_vel != NULL)
    {
    	for (size_t sh = 0; sh < scvf.num_sh(); sh++)
        {
            VecScale(D_shape_vel[sh], scvf.normal(), scvf.shape(sh, ip));
        }
    }

//	The shapes do not depend of the diffusion tensor
    if(bNonZeroDerivD != NULL) bNonZeroDerivD = false;

//	we're done
    return true;
}

/**
 * \param[in]		fvg				Finite Volume Geometry
 * \param[in]		scvf			SubControlVolumeFace to consider (ip)
 * \param[in]		vel				Convection Velocity at ip
 * \param[in]		Diffusion		Diffusion Tensor at ip
 * \param[out]		shape			Convection shape for each dof
 * \param[out]		D_shape_vel		Derivative of the shape w.r.t. the Velocity for each dof
 * \param[out]		D_shape_D		Derivative of the shape w.r.t. the Diffusion for each dof
 * \param[out]		bNonZeroDerivD	nonzero iff the deriv. w.r.t. the diffusion tensor is non-zero
 *
 */
template <typename TFVGeometry>
bool ConvectionShapes_FullUpwind
(
	const TFVGeometry& fvg,
    const typename TFVGeometry::SCVF& scvf,
    const MathVector<TFVGeometry::world_dim>& Vel,
    const MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>& Diffusion,
    number* shape,
    MathVector<TFVGeometry::world_dim>* D_shape_vel,
    MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>* D_shape_D,
    bool* bNonZeroDerivD
)
{
//	Currently only first order
	UG_ASSERT(scvf.num_ip() == 1, "Only first order implemented.");

//	Compute flux
	const number flux = VecDot(scvf.normal(), Vel);

//	Choose Upwind corner
	size_t up = (flux >= 0) ? scvf.from() : scvf.to();

//	Write Shapes
    memset(shape, 0, scvf.num_sh() * sizeof (number));
    shape[up] = flux;

//	Write Derivatives if wanted
	if (D_shape_vel != NULL)
	{
		memset(D_shape_vel, 0, scvf.num_sh() * sizeof ( MathVector<TFVGeometry::world_dim> ));
		D_shape_vel [up] = scvf.normal();
    }

//	The shapes do not depend of the diffusion tensor
    if(bNonZeroDerivD != NULL) bNonZeroDerivD = false;

//	we're done
	return true;
}

/*** Convection shapes for "partial upwind": ***/

/* GetEffDiffusionH - computes the effective scalar diffusion multiplyed
 * by the distance between the nodes for the computation of the Peclet
 * number in the partial and the exponential upwind methods. Additionally,
 * the derivatives of this coefficient w.r.t. the diffusion tensor are
 * computed. Note that the derivatives are computed only if the coefficient
 * is > 0.
 * Remark: In the papers by Frolkovic, this coefficient would be
 * denoted by $\lambda^{e}_{ij} \cdot | \Gamma^{e}_{ij} | / h_{ij}$.
 * The function returns 0 if OK, non-zero on an error:
 */
template <typename TFVGeometry>
bool GetEffDiffusion
(
	const TFVGeometry& fvg,
	const typename TFVGeometry::SCVF& scvf,
	const MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>& Diffusion,
	number& lambda, /* to save the coefficient */
	MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>& deriv /* the derivatives */
)
{
//	Currently only first order
	static const size_t ip = 0;
	UG_ASSERT(scvf.num_ip() == 1, "Only first order implemented.");

//	Dimension
	static const int dim = TFVGeometry::world_dim;

//	Compute Volume of Element
	number Volume = ElementSize<typename TFVGeometry::ref_elem_type, dim>(fvg.corners());

//  Get Gradients
	MathVector<dim> DiffGrad;
	const MathVector<dim>& gradTo = scvf.global_grad(scvf.to(), ip);
	const MathVector<dim>& gradFrom = scvf.global_grad(scvf.from(), ip);

//	Compute DiffGrad = D * Grad Phi_to
	MatVecMult(DiffGrad, Diffusion, gradTo);

//	Compute GradDiffGrad = < Grad Phi_from, DiffGrad >
	number GradDiffGrad = VecDot(DiffGrad,  gradFrom);

//	Set lambda
	if(! ((lambda = - GradDiffGrad * Volume) > 0) )
		return true;

//	Compute derivatives
	for (size_t i = 0; i < (size_t)dim; i++)
		for (size_t j = 0; j < (size_t)dim; j++)
			deriv(i, j) = - gradFrom[i] * gradTo[j] * Volume;

//	we're done
	return true;
}

/**
 * \param[in]		fvg				Finite Volume Geometry
 * \param[in]		scvf			SubControlVolumeFace to consider (ip)
 * \param[in]		vel				Convection Velocity at ip
 * \param[in]		Diffusion		Diffusion Tensor at ip
 * \param[out]		shape			Convection shape for each dof
 * \param[out]		D_shape_vel		Derivative of the shape w.r.t. the Velocity for each dof
 * \param[out]		D_shape_D		Derivative of the shape w.r.t. the Diffusion for each dof
 * \param[out]		bNonZeroDerivD	nonzero iff the deriv. w.r.t. the diffusion tensor is non-zero
 *
 */
template <typename TFVGeometry>
bool ConvectionShapes_PartialUpwind
(
	const TFVGeometry& fvg,
    const typename TFVGeometry::SCVF& scvf,
    const MathVector<TFVGeometry::world_dim>& Vel,
    const MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>& Diffusion,
    number* shape,
    MathVector<TFVGeometry::world_dim>* D_shape_vel,
    MathMatrix<TFVGeometry::world_dim, TFVGeometry::world_dim>* D_shape_D,
    bool* bNonZeroDerivD
)
{
//	Currently only first order
	UG_ASSERT(scvf.num_ip() == 1, "Only first order implemented.");
	static const int dim = TFVGeometry::world_dim;

//	Get Effective Diffusion
	MathMatrix<dim, dim> D_lambda;
	number lambda;
	if (!GetEffDiffusion (fvg, scvf, Diffusion, lambda, D_lambda))
		return false;

//	The case of the "non-positive diffusive flux" (lambda <= 0)
	if (lambda <= 0)
		return ConvectionShapes_FullUpwind (fvg, scvf, Vel, Diffusion, shape,
											D_shape_vel, D_shape_D, bNonZeroDerivD);

//	Compute Convective Flux
	number convFlux = VecDot(scvf.normal(), Vel);

//	get corners
    const size_t co_from = scvf.from();
    const size_t co_to = scvf.to();

    memset (shape, 0, scvf.num_sh() * sizeof (number));

//	The case of the diffusion dominance (central differences)
	if (2 * lambda > fabs(convFlux))
	{
		if (bNonZeroDerivD != NULL) *bNonZeroDerivD = false;

		shape[co_from] = convFlux / 2;
		shape[co_to] = convFlux / 2;

		if (D_shape_vel != NULL)
		{
			memset (D_shape_vel, 0, scvf.num_sh() * sizeof (MathVector<dim>));

			VecScale(D_shape_vel[co_from], scvf.normal(), 1.0/2.0);
			VecScale(D_shape_vel[co_to], scvf.normal(), 1.0/2.0);
		}

		return true;
	}

//	The cases of the convection dominance
	if (bNonZeroDerivD != NULL) *bNonZeroDerivD = true;
	if (convFlux >= 0)
	{
		shape[co_from] = convFlux - lambda;
		shape[co_to] = lambda;

		if (D_shape_vel != NULL)
		{
			memset (D_shape_vel, 0, scvf.num_sh() * sizeof (MathVector<dim>));
			D_shape_vel[co_from] = scvf.normal();
		}
	}
	else
	{
		shape[co_from] = - lambda;
		shape[co_to] = convFlux + lambda;

		if (D_shape_vel != NULL)
		{
			memset (D_shape_vel, 0, scvf.num_sh() * sizeof (MathVector<dim>));
			D_shape_vel[co_to] = scvf.normal();
		}
	}

	if (D_shape_D != NULL)
	{
		memset (D_shape_D, 0, scvf.num_sh() * sizeof (MathMatrix<dim,dim>));
		for (size_t i = 0; i < (size_t)dim; i++)
			for (size_t j = 0; j < (size_t)dim; j++)
			{
				D_shape_D[co_from](i,j) = - D_lambda(i,j);
				D_shape_D[co_to](i,j) = D_lambda(i,j);
			}
	}

	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONVECTION_SHAPE__ */
