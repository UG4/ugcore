/*
 * conv_shape.h
 *
 *  Created on: 08.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__

#include "conv_shape_interface.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim, typename TAlgebra>
class ConvectionShapesNoUpwind
	: public IConvectionShapes<TDim, TAlgebra>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim, TAlgebra> base_type;

	///	This class
		typedef ConvectionShapesNoUpwind<TDim, TAlgebra> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::conv_shape_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesNoUpwind()
		{
		//	the shapes do not depend on the Diffusion. Thus, we can set the
		//	derivative to be always zero w.r.t. the Diffusion for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo,
		            const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
		            const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
		            bool computeDeriv);

	private:
		void register_func(Int2Type<1>)
		{register_func<Edge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func<Triangle>();
			register_func<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(const TGeom* geo,
					const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
					const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
					bool computeDeriv);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
ConvectionShapesNoUpwind<TDim, TAlgebra>::
update(const FV1Geometry<TElem, dim>* geo,
       const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
       const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
       bool computeDeriv)
{
//	loop subcontrol volume faces
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get subcontrol volume face
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//	Currently only first order
		static const size_t ip = 0;
		UG_ASSERT(scvf.num_ip() == 1, "Only first order implemented.");

	//	Compute flux
		const number flux = VecDot(scvf.normal(), DarcyVelocity[i]);

	//	Write Shapes
		for(size_t sh = 0; sh < scvf.num_sh(); sh++)
			conv_shape(i, sh) = flux * scvf.shape(sh, ip);

	//	Write Derivatives if wanted
		if(computeDeriv)
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				VecScale(conv_shape_vel(i, sh),
				         scvf.normal(), scvf.shape(sh, ip));

	//	The shapes do not depend of the diffusion tensor
	}

//	we're done
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim, typename TAlgebra>
class ConvectionShapesFullUpwind
	: public IConvectionShapes<TDim, TAlgebra>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim, TAlgebra> base_type;

	///	This class
		typedef ConvectionShapesFullUpwind<TDim, TAlgebra> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::conv_shape_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesFullUpwind()
		{
		//	the shapes do not depend on the Diffusion. Thus, we can set the
		//	derivative to be always zero w.r.t. the Diffusion for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo,
		            const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
		            const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
		            bool computeDeriv);

	private:
		void register_func(Int2Type<1>)
		{register_func<Edge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func<Triangle>();
			register_func<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(const TGeom* geo,
					const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
					const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
					bool computeDeriv);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
ConvectionShapesFullUpwind<TDim, TAlgebra>::
update(const FV1Geometry<TElem, dim>* geo,
       const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
       const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
       bool computeDeriv)
{
//	loop subcontrol volume faces
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get subcontrol volume face
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//	Compute flux
		const number flux = VecDot(scvf.normal(), DarcyVelocity[i]);

	//	Choose Upwind corner
		const size_t up = (flux >= 0) ? scvf.from() : scvf.to();

	//	Write Shapes
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh) conv_shape(i, sh) = 0.0;
		conv_shape(i, up) = flux;

	//	Write Derivatives if wanted
		if(computeDeriv)
		{
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecSet(conv_shape_vel(i, sh), 0.0);
			conv_shape_vel(i, up) = scvf.normal();
		}

	//	The shapes do not depend of the diffusion tensor
	}

//	we're done
	return true;
}


/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim, typename TAlgebra>
class ConvectionShapesPartialUpwind
	: public IConvectionShapes<TDim, TAlgebra>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim, TAlgebra> base_type;

	///	This class
		typedef ConvectionShapesPartialUpwind<TDim, TAlgebra> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::conv_shape_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesPartialUpwind()
		{
		//	the shapes do not depend on the Diffusion. Thus, we can set the
		//	derivative to be always zero w.r.t. the Diffusion for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo,
		            const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
		            const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
		            bool computeDeriv);

	private:
		void register_func(Int2Type<1>)
		{register_func<Edge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func<Triangle>();
			register_func<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(const TGeom* geo,
					const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
					const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
					bool computeDeriv);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

template <int TDim, typename TAlgebra>
template <typename TElem>
bool
ConvectionShapesPartialUpwind<TDim, TAlgebra>::
update(const FV1Geometry<TElem, dim>* geo,
       const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
       const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
       bool computeDeriv)
{
//	loop subcontrol volume faces
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get subcontrol volume face
		const typename FV1Geometry<TElem, dim>::SCVF& scvf = geo->scvf(i);

	//	Currently only first order
		static const size_t ip = 0;
		UG_ASSERT(scvf.num_ip() == 1, "Only first order implemented.");

	//	Compute Volume of Element
		number Volume = ElementSize<typename FV1Geometry<TElem, dim>::ref_elem_type, dim>(geo->corners());

	//  Get Gradients
		MathVector<dim> DiffGrad;
		const MathVector<dim>& gradTo = scvf.global_grad(scvf.to(), ip);
		const MathVector<dim>& gradFrom = scvf.global_grad(scvf.from(), ip);

	//	Compute DiffGrad = D * Grad Phi_to
		MatVecMult(DiffGrad, Diffusion[i], gradTo);

	//	Compute GradDiffGrad = < Grad Phi_from, DiffGrad >
		const number GradDiffGrad = VecDot(DiffGrad,  gradFrom);

	//	Set lambda
		const number lambda = - GradDiffGrad * Volume;

	//	Compute flux
		const number flux = VecDot(scvf.normal(), DarcyVelocity[i]);

	//	switch if full upwind has to be used
		if(lambda <= 0)
		{
		//	Choose Upwind corner
			const size_t up = (flux >= 0) ? scvf.from() : scvf.to();

		//	Write Shapes
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh) conv_shape(i, sh) = 0.0;
			conv_shape(i, up) = flux;

		//	Write Derivatives if wanted
			if(computeDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecSet(conv_shape_vel(i, sh), 0.0);
				conv_shape_vel(i, up) = scvf.normal();

			//	does not depend on diffusion
				set_non_zero_deriv_diffusion_flag(false);
			}

		//	everything done
			return true;
		}

	//	Compute derivatives
		MathMatrix<dim, dim> D_lambda;
		for (size_t k = 0; k < (size_t)dim; k++)
			for (size_t l = 0; l < (size_t)dim; l++)
				D_lambda(k, l) = - gradFrom[k] * gradTo[l] * Volume;

	//	get corners
		const size_t co_from = scvf.from();
		const size_t co_to = scvf.to();


	//	reset values
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh) conv_shape(i, sh) = 0.0;

	//	The case of the diffusion dominance (central differences)
		if (2 * lambda > fabs(flux))
		{
			conv_shape(i, co_from) = flux / 2;
			conv_shape(i, co_to) = flux / 2;

			if(computeDeriv)
			{
				set_non_zero_deriv_diffusion_flag(false);

				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecSet(conv_shape_vel(i, sh), 0.0);

				VecScale(conv_shape_vel(i,co_from), scvf.normal(), 1.0/2.0);
				VecScale(conv_shape_vel(i, co_to), scvf.normal(), 1.0/2.0);
			}

		//	everything done
			return true;
		}

	//	The cases of the convection dominance
		set_non_zero_deriv_diffusion_flag(true);
		if (flux >= 0)
		{
			conv_shape(i, co_from) = flux - lambda;
			conv_shape(i, co_to) = lambda;

			if(computeDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecSet(conv_shape_vel(i, sh), 0.0);

				conv_shape_vel(i,co_from) = scvf.normal();
			}
		}
		else
		{
			conv_shape(i, co_from) = - lambda;
			conv_shape(i, co_to) = flux + lambda;

			if(computeDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecSet(conv_shape_vel(i, sh), 0.0);

				conv_shape_vel(i,co_to) = scvf.normal();
			}
		}

		if (computeDeriv)
		{
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatSet(conv_shape_diffusion(i, sh), 0.0);

			for (size_t k = 0; k < (size_t)dim; k++)
				for (size_t l = 0; l < (size_t)dim; l++)
				{
					conv_shape_diffusion(i, co_from)(k,l) = - D_lambda(k,l);
					conv_shape_diffusion(i, co_to)(k,l) = D_lambda(k,l);
				}
		}
	}

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__ */
