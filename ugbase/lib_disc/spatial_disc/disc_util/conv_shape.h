
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__

#include "conv_shape_interface.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class ConvectionShapesNoUpwind
	: public IConvectionShapes<TDim>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim> base_type;

	///	This class
		typedef ConvectionShapesNoUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::D_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesNoUpwind()
		{
		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv);

	private:
		void register_func(Int2Type<1>)
		{	register_func_for_refDim<1>();
			register_func_for_elem<RegularEdge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func_for_refDim<2>();
			register_func_for_elem<Triangle>();
			register_func_for_elem<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func_for_refDim<3>();
			register_func_for_elem<Tetrahedron>();
			register_func_for_elem<Pyramid>();
			register_func_for_elem<Prism>();
			register_func_for_elem<Hexahedron>();
			register_func_for_elem<Octahedron>();}

		template <typename TElem>
		void register_func_for_elem()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);

			typedef HFV1Geometry<TElem, dim> THGeom;
			typedef bool (this_type::*THFunc)
					(  const THGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<THGeom, THFunc>(&this_type::template update<THGeom>);
		}
		
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
		}
};

template <int TDim>
template <typename TFVGeom>
bool
ConvectionShapesNoUpwind<TDim>::
update(const TFVGeom* geo,
       const MathVector<dim>* Velocity,
       const MathMatrix<dim, dim>* DiffDisp,
       bool computeDeriv)
{
	UG_ASSERT(geo != NULL, "Null pointer");
	UG_ASSERT(Velocity != NULL, "Null pointer");

//	\todo: think about: this should be something like scvf.num_sh()
	const size_t numSH = geo->num_sh();

//	loop subcontrol volume faces
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get subcontrol volume face
		const typename TFVGeom::SCVF& scvf = geo->scvf(ip);

	//	Compute flux
		const number flux = VecDot(scvf.normal(), Velocity[ip]);

	//	Write Shapes
		for(size_t sh = 0; sh < scvf.num_sh(); sh++)
			conv_shape(ip, sh) = flux * scvf.shape(sh);

	//	this is introduced here, hopefully temporarily: The problem is, that
	//	for hanging nodes the number of shape function is not the number of
	//	corners, but scvf.num_sh() currently returns the number of corners.
	//	this is actually enough to interpolate the function, but still we
	//	should reset the interpolation adding for hanging dofs to zero
		for(size_t sh = scvf.num_sh(); sh < numSH; sh++)
			conv_shape(ip, sh) = 0.0;

	//	Write Derivatives if wanted
		if(computeDeriv){
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				VecScale(D_vel(ip, sh), scvf.normal(), scvf.shape(sh));

			// temporary, see above
			for(size_t sh = scvf.num_sh(); sh < numSH; sh++)
				VecSet(D_vel(ip, sh), 0.0);
		}

	//	The shapes do not depend of the diffusion tensor
	}

//	we're done
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class ConvectionShapesFullUpwind
	: public IConvectionShapes<TDim>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim> base_type;

	///	This class
		typedef ConvectionShapesFullUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::D_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesFullUpwind()
		{
		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv);

	private:
		void register_func(Int2Type<1>)
		{	register_func_for_refDim<1>();
			register_func_for_elem<RegularEdge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func_for_refDim<2>();
			register_func_for_elem<Triangle>();
			register_func_for_elem<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func_for_refDim<3>();
			register_func_for_elem<Tetrahedron>();
			register_func_for_elem<Pyramid>();
			register_func_for_elem<Prism>();
			register_func_for_elem<Hexahedron>();
			register_func_for_elem<Octahedron>();}

		template <typename TElem>
		void register_func_for_elem()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);

			typedef HFV1Geometry<TElem, dim> THGeom;
			typedef bool (this_type::*THFunc)
					(  const THGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<THGeom, THFunc>(&this_type::template update<THGeom>);
		}
		
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
		}
};

template <int TDim>
template <typename TFVGeom>
bool ConvectionShapesFullUpwind<TDim>::
update(const TFVGeom* geo,
       const MathVector<dim>* Velocity,
       const MathMatrix<dim, dim>* DiffDisp,
       bool computeDeriv)
{
	UG_ASSERT(geo != NULL, "Null pointer");
	UG_ASSERT(Velocity != NULL, "Null pointer");

//	\todo: think about: this should be something like scvf.num_sh()
	const size_t numSH = geo->num_sh();

//	loop subcontrol volume faces
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get subcontrol volume face
		const typename TFVGeom::SCVF& scvf = geo->scvf(ip);

	//	Compute flux
		const number flux = VecDot(scvf.normal(), Velocity[ip]);

		size_t from = scvf.from();
		size_t to = scvf.to();

	//	if e.g. hanging nodes are involved, no upwind can be performed...		
		if((from >= scvf.num_sh()) || (to >= scvf.num_sh())){
		//	No upwind...
		//	Write Shapes
			for(size_t sh = 0; sh < scvf.num_sh(); sh++)
				conv_shape(ip, sh) = flux * scvf.shape(sh);

			UG_ASSERT(scvf.num_sh() == numSH, "sh's have to match!");

		//	Write Derivatives if wanted
			if(computeDeriv){
				for (size_t sh = 0; sh < scvf.num_sh(); sh++)
					VecScale(D_vel(ip, sh), scvf.normal(), scvf.shape(sh));
			}
			
			continue;
		}
		
	//	Full upwind below...
	//	Choose Upwind corner
		size_t up = (flux >= 0) ? from : to;

	//	Write Shapes
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh) conv_shape(ip, sh) = 0.0;

		conv_shape(ip, up) = flux;

	//	Write Derivatives if wanted
		if(computeDeriv)
		{
			for(size_t sh = 0; sh < numSH; ++sh) VecSet(D_vel(ip, sh), 0.0);
			D_vel(ip, up) = scvf.normal();
		}

	//	The shapes do not depend of the diffusion tensor
	}

//	we're done
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Weighted Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class ConvectionShapesWeightedUpwind
	: public IConvectionShapes<TDim>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim> base_type;

	///	This class
		typedef ConvectionShapesWeightedUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::D_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesWeightedUpwind() : m_weight(0.5)
		{
		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	constructor
		ConvectionShapesWeightedUpwind(number weight)
		{
			set_weight(weight);

		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	set weighting between full upwind (1.0) and no upwind (0.0)
		void set_weight(number weight) {m_weight = weight;}

	///	update of values for FV1Geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv);

	private:
	//	weight between no and full upwind (1.0 -> full upwind, 0.0 -> no upwind)
		number m_weight;

		void register_func(Int2Type<1>)
		{	register_func_for_refDim<1>();
			register_func_for_elem<RegularEdge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func_for_refDim<2>();
			register_func_for_elem<Triangle>();
			register_func_for_elem<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func_for_refDim<3>();
			register_func_for_elem<Tetrahedron>();
			register_func_for_elem<Pyramid>();
			register_func_for_elem<Prism>();
			register_func_for_elem<Hexahedron>();
			register_func_for_elem<Octahedron>();}

		template <typename TElem>
		void register_func_for_elem()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);

			typedef HFV1Geometry<TElem, dim> THGeom;
			typedef bool (this_type::*THFunc)
					(  const THGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<THGeom, THFunc>(&this_type::template update<THGeom>);
		}
		
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
		}
};

template <int TDim>
template <typename TFVGeom>
bool ConvectionShapesWeightedUpwind<TDim>::
update(const TFVGeom* geo,
       const MathVector<dim>* Velocity,
       const MathMatrix<dim, dim>* DiffDisp,
       bool computeDeriv)
{
	UG_ASSERT(geo != NULL, "Null pointer");
	UG_ASSERT(Velocity != NULL, "Null pointer");

//	\todo: think about: this should be something like scvf.num_sh()
	const size_t numSH = geo->num_sh();

//	loop subcontrol volume faces
	for(size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
	//	get subcontrol volume face
		const typename TFVGeom::SCVF& scvf = geo->scvf(ip);

	//	Compute flux
		const number flux = VecDot(scvf.normal(), Velocity[ip]);

		size_t from = scvf.from();
		size_t to = scvf.to();

	//	if e.g. hanging nodes are involved, no upwind can be performed...		
		if((from >= scvf.num_sh()) || (to >= scvf.num_sh())){
		//	No upwind...
		//	Write Shapes
			for(size_t sh = 0; sh < scvf.num_sh(); sh++)
				conv_shape(ip, sh) = flux * scvf.shape(sh);

			UG_ASSERT(scvf.num_sh() == numSH, "sh's have to match!");

		//	Write Derivatives if wanted
			if(computeDeriv){
				for (size_t sh = 0; sh < scvf.num_sh(); sh++)
					VecScale(D_vel(ip, sh), scvf.normal(), scvf.shape(sh));
			}
			
			continue;
		}
		
	//	Choose Upwind corner
		size_t up = (flux >= 0) ? from : to;
		
	//	write no upwind part of shapes
		const number noUpFlux = (1.-m_weight)*flux;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			conv_shape(ip, sh) = noUpFlux * scvf.shape(sh);

	//	add full upwind part of shapes
		conv_shape(ip, up) += m_weight * flux;

	//	Write Derivatives if wanted
		if(computeDeriv)
		{
		//	write no upwind part of derivatives
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				VecScale(D_vel(ip, sh), scvf.normal(),
				         	 	 	 	 	 (1.-m_weight)*scvf.shape(sh));
		//	see comment above
			for(size_t sh = scvf.num_sh(); sh < numSH; sh++)
				VecSet(D_vel(ip, sh), 0.0);

		//	add full upwind part of derivatives
			VecScaleAppend(D_vel(ip, up), m_weight, scvf.normal());
		}

	//	The shapes do not depend of the diffusion tensor
	}

//	we're done
	return true;
}


/////////////////////////////////////////////////////////////////////////////
// Partial Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class ConvectionShapesPartialUpwind
	: public IConvectionShapes<TDim>
{
	public:
	///	Base class
		typedef IConvectionShapes<TDim> base_type;

	///	This class
		typedef ConvectionShapesPartialUpwind<TDim> this_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_non_zero_deriv_diffusion_flag;
		using base_type::conv_shape;
		using base_type::D_vel;
		using base_type::conv_shape_diffusion;
		using base_type::non_zero_deriv_diffusion;
		using base_type::register_update_func;

	public:
	///	constructor
		ConvectionShapesPartialUpwind()
		{
		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
					bool computeDeriv);

	private:
		void register_func(Int2Type<1>)
		{	register_func_for_refDim<1>();
			register_func_for_elem<RegularEdge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func_for_refDim<2>();
			register_func_for_elem<Triangle>();
			register_func_for_elem<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func_for_refDim<3>();
			register_func_for_elem<Tetrahedron>();
			register_func_for_elem<Pyramid>();
			register_func_for_elem<Prism>();
			register_func_for_elem<Hexahedron>();
			register_func_for_elem<Octahedron>();}

		template <typename TElem>
		void register_func_for_elem()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);

		}

		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc)
					(  const TGeom* geo,
					   const MathVector<dim>* Velocity,
					   const MathMatrix<dim, dim>* DiffDisp,
					   bool computeDeriv);

			base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
		}
};

template <int TDim>
template <typename TFVGeom>
bool
ConvectionShapesPartialUpwind<TDim>::
update(const TFVGeom* geo,
       const MathVector<dim>* Velocity,
       const MathMatrix<dim, dim>* DiffDisp,
       bool computeDeriv)
{
	UG_ASSERT(geo != NULL, "Null pointer");
	UG_ASSERT(Velocity != NULL, "Null pointer");
//	UG_ASSERT(DiffDisp != NULL, "Null pointer");

//	Compute Volume of Element
//	typedef typename TFVGeom::ref_elem_type ref_elem_type;
	const number vol = ElementSize<dim>(geo->roid(), geo->corners());

//	loop subcontrol volume faces
	for(size_t i = 0; i < geo->num_scvf(); ++i)
	{
	//	get subcontrol volume face
		const typename TFVGeom::SCVF& scvf = geo->scvf(i);

	//	get corners
		const size_t from = scvf.from();
		const size_t to = scvf.to();

	//	get gradients
		const MathVector<dim>& gradTo = scvf.global_grad(to);
		const MathVector<dim>& gradFrom = scvf.global_grad(from);

	//	set lambda negative as default
		number lambda = -1;

	//	if DiffDisp-Tensor passed, compute lambda
		if(DiffDisp != NULL)
		{
		//  Get Gradients
			MathVector<dim> DiffGrad;

		//	Compute DiffGrad = D * Grad Phi_to
			MatVecMult(DiffGrad, DiffDisp[i], gradTo);

		//	Compute GradDiffGrad = < Grad Phi_from, DiffGrad >
			const number GradDiffGrad = VecDot(DiffGrad,  gradFrom);

		//	Set lambda
			lambda = - GradDiffGrad * vol;
		}

	//	Compute flux
		const number flux = VecDot(scvf.normal(), Velocity[i]);

	//	reset values
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			conv_shape(i, sh) = 0.0;
		if(computeDeriv)
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecSet(D_vel(i, sh), 0.0);


	//	check: Currently hanging nodes not supported.
	//	/todo: add support for hanging nodes
		if(from >= scvf.num_sh() || to >= scvf.num_sh())
			UG_THROW("PartialUpwind: Currently not implemented for hanging nodes.")

	///////////////////////////////////////////////////////////////////
	//	Case 1:
	//	full upwind is used
	///////////////////////////////////////////////////////////////////
		if(lambda <= 0 || DiffDisp == NULL)
		{
		//	Choose Upwind corner
			const size_t up = (flux >= 0) ? scvf.from() : scvf.to();

		//	Write Shapes
			conv_shape(i, up) = flux;

		//	Write Derivatives if wanted
			if(computeDeriv)
			{
			//	set derivative
				D_vel(i, up) = scvf.normal();

			//	does not depend on diffusion
				set_non_zero_deriv_diffusion_flag(false);
			}

		//	everything done
			continue;
		}

	///////////////////////////////////////////////////////////////////
	//	Case 2:
	//	The case of the diffusion dominance (central differences)
	///////////////////////////////////////////////////////////////////
		if (2 * lambda > fabs(flux))
		{
			conv_shape(i, from) = flux / 2.0;
			conv_shape(i, to) = flux / 2.0;

			if(computeDeriv)
			{
				set_non_zero_deriv_diffusion_flag(false);

				VecScale(D_vel(i,from), scvf.normal(), 1.0/2.0);
				VecScale(D_vel(i, to), scvf.normal(), 1.0/2.0);
			}

		//	everything done
			continue;
		}

	///////////////////////////////////////////////////////////////////
	//	Case 3:
	//	The cases of the convection dominance
	///////////////////////////////////////////////////////////////////
		set_non_zero_deriv_diffusion_flag(true);
		if (flux >= 0)
		{
			conv_shape(i, from) = flux - lambda;
			conv_shape(i, to) = lambda;

			if(computeDeriv)
				D_vel(i,from) = scvf.normal();
		}
		else
		{
			conv_shape(i, from) = - lambda;
			conv_shape(i, to) = flux + lambda;

			if(computeDeriv)
				D_vel(i,to) = scvf.normal();
		}

		if (computeDeriv)
		{
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatSet(conv_shape_diffusion(i, sh), 0.0);

			for (size_t k = 0; k < (size_t)dim; k++)
				for (size_t l = 0; l < (size_t)dim; l++)
				{
					conv_shape_diffusion(i, from)(k,l) = gradFrom[k]*gradTo[l]*vol;
					conv_shape_diffusion(i, to)(k,l) = - gradFrom[k]*gradTo[l]*vol;
				}
		}
	}

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__ */
