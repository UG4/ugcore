/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__

#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

// ug4 headers
#include "conv_shape_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
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
			boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
			boost::mpl::for_each< boost::mpl::range_c<int, 1, dim+1> >(RegisterRefDimFunc(this));
		}

	///	update of values for FV geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv);

	private:
		
	///	functor for registering the shapes for the element-templated FV geometries
		struct RegisterElemFunc
		{
			this_type* m_pThis;
			RegisterElemFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TElem> void operator() (TElem &)
			{
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1Geometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1CondensedGeometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, HFV1Geometry<TElem, dim> >();
			}
		};
	
	/// registers the update function for an element type and a FV geometry
		template <typename TElem, typename TFVGeom>
		void register_func_for_elem_fvgeom()
		{
			typedef bool (this_type::*TFunc) (const TFVGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
			base_type::template register_update_func<TFVGeom, TFunc>(&this_type::template update<TFVGeom>);
		}
		
	///	functor for registering the shapes for the reference-dimension-templated FV geometries
		struct RegisterRefDimFunc
		{
			this_type* m_pThis;
			RegisterRefDimFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TRefDim> void operator() (TRefDim &) {m_pThis->register_func_for_refDim<TRefDim::value> ();}
		};

	/// registers the update function for a reference dimension
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc) (const TGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
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
	UG_ASSERT(geo != nullptr, "Null pointer");
	UG_ASSERT(Velocity != nullptr, "Null pointer");

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
			boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
			boost::mpl::for_each< boost::mpl::range_c<int, 1, dim+1> >(RegisterRefDimFunc(this));
		}

	///	update of values for FV1Geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv);

	private:
		
	///	functor for registering the shapes for the element-templated FV geometries
		struct RegisterElemFunc
		{
			this_type* m_pThis;
			RegisterElemFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TElem> void operator() (TElem &)
			{
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1Geometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1CondensedGeometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, HFV1Geometry<TElem, dim> >();
			}
		};
	
	/// registers the update function for an element type and a FV geometry
		template <typename TElem, typename TFVGeom>
		void register_func_for_elem_fvgeom()
		{
			typedef bool (this_type::*TFunc) (const TFVGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
			base_type::template register_update_func<TFVGeom, TFunc>(&this_type::template update<TFVGeom>);
		}
		
	///	functor for registering the shapes for the reference-dimension-templated FV geometries
		struct RegisterRefDimFunc
		{
			this_type* m_pThis;
			RegisterRefDimFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TRefDim> void operator() (TRefDim &) {m_pThis->register_func_for_refDim<TRefDim::value> ();}
		};

	/// registers the update function for a reference dimension
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc) (const TGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
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
	UG_ASSERT(geo != nullptr, "Null pointer");
	UG_ASSERT(Velocity != nullptr, "Null pointer");

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
			boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
			boost::mpl::for_each< boost::mpl::range_c<int, 1, dim+1> >(RegisterRefDimFunc(this));
		}

	///	constructor
		ConvectionShapesWeightedUpwind(number weight)
		{
			set_weight(weight);

		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
			set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
			boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
			boost::mpl::for_each< boost::mpl::range_c<int, 1, dim+1> >(RegisterRefDimFunc(this));
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
		
	private:
	
	///	functor for registering the shapes for the element-templated FV geometries
		struct RegisterElemFunc
		{
			this_type* m_pThis;
			RegisterElemFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TElem> void operator() (TElem &)
			{
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1Geometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1CondensedGeometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, HFV1Geometry<TElem, dim> >();
			}
		};
	
	/// registers the update function for an element type and a FV geometry
		template <typename TElem, typename TFVGeom>
		void register_func_for_elem_fvgeom()
		{
			typedef bool (this_type::*TFunc) (const TFVGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
			base_type::template register_update_func<TFVGeom, TFunc>(&this_type::template update<TFVGeom>);
		}
		
	///	functor for registering the shapes for the reference-dimension-templated FV geometries
		struct RegisterRefDimFunc
		{
			this_type* m_pThis;
			RegisterRefDimFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TRefDim> void operator() (TRefDim &) {m_pThis->register_func_for_refDim<TRefDim::value> ();}
		};

	/// registers the update function for a reference dimension
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc) (const TGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
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
	UG_ASSERT(geo != nullptr, "Null pointer");
	UG_ASSERT(Velocity != nullptr, "Null pointer");

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
			boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
			boost::mpl::for_each< boost::mpl::range_c<int, 1, dim+1> >(RegisterRefDimFunc(this));
		}

	///	update of values for FV1Geometry
		template <typename TFVGeom>
		bool update(const TFVGeom* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
					bool computeDeriv);

	private:
		
	///	functor for registering the shapes for the element-templated FV geometries
		struct RegisterElemFunc
		{
			this_type* m_pThis;
			RegisterElemFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TElem> void operator() (TElem &)
			{
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1Geometry<TElem, dim> >();
				m_pThis->template register_func_for_elem_fvgeom< TElem, FV1CondensedGeometry<TElem, dim> >();
			}
		};
	
	/// registers the update function for an element type and a FV geometry
		template <typename TElem, typename TFVGeom>
		void register_func_for_elem_fvgeom()
		{
			typedef bool (this_type::*TFunc) (const TFVGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
			base_type::template register_update_func<TFVGeom, TFunc>(&this_type::template update<TFVGeom>);
		}
		
	///	functor for registering the shapes for the reference-dimension-templated FV geometries
		struct RegisterRefDimFunc
		{
			this_type* m_pThis;
			RegisterRefDimFunc(this_type * pThis) : m_pThis(pThis) {}
			template<typename TRefDim> void operator() (TRefDim &) {m_pThis->register_func_for_refDim<TRefDim::value> ();}
		};

	/// registers the update function for a reference dimension
		template <int refDim>
		void register_func_for_refDim()
		{
			typedef DimFV1Geometry<refDim, dim> TGeom;
			typedef bool (this_type::*TFunc) (const TGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
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
	UG_ASSERT(geo != nullptr, "Null pointer");
	UG_ASSERT(Velocity != nullptr, "Null pointer");
//	UG_ASSERT(DiffDisp != nullptr, "Null pointer");

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
		if(DiffDisp != nullptr)
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
		if(lambda <= 0 || DiffDisp == nullptr)
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

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class ConvectionShapesSkewedUpwind
	: public IConvectionShapes<TDim>
{
public:
	///	Base class
	typedef IConvectionShapes<TDim> base_type;

	///	This class
	typedef ConvectionShapesSkewedUpwind<TDim> this_type;

	///	Dimension
	static const int dim = TDim;

protected:
	//	explicitly forward some function
	using base_type::conv_shape;
	using base_type::conv_shape_diffusion;
	using base_type::D_vel;
	using base_type::non_zero_deriv_diffusion;
	using base_type::register_update_func;
	using base_type::set_non_zero_deriv_diffusion_flag;

public:
	///	constructor
	ConvectionShapesSkewedUpwind()
	{
		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
		set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
		boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
		boost::mpl::for_each<boost::mpl::range_c<int, 1, dim + 1>>(RegisterRefDimFunc(this));
	}

	///	update of values for FV geometry
	template <typename TFVGeom>
	bool update(const TFVGeom *geo,
				const MathVector<dim> *Velocity,
				const MathMatrix<dim, dim> *DiffDisp,
				bool computeDeriv);

private:
	///	functor for registering the shapes for the element-templated FV geometries
	struct RegisterElemFunc
	{
		this_type *m_pThis;
		RegisterElemFunc(this_type *pThis) : m_pThis(pThis) {}
		template <typename TElem>
		void operator()(TElem &)
		{
			m_pThis->template register_func_for_elem_fvgeom<TElem, FV1Geometry<TElem, dim>>();
			//m_pThis->template register_func_for_elem_fvgeom< TElem, FV1CondensedGeometry<TElem, dim> >();
			//m_pThis->template register_func_for_elem_fvgeom< TElem, HFV1Geometry<TElem, dim> >();
		}
	};

	/// registers the update function for an element type and a FV geometry
	template <typename TElem, typename TFVGeom>
	void register_func_for_elem_fvgeom()
	{
		typedef bool (this_type::*TFunc)(const TFVGeom *, const MathVector<dim> *, const MathMatrix<dim, dim> *, bool);
		base_type::template register_update_func<TFVGeom, TFunc>(&this_type::template update<TFVGeom>);
	}

	///	functor for registering the shapes for the reference-dimension-templated FV geometries
	struct RegisterRefDimFunc
	{
		this_type *m_pThis;
		RegisterRefDimFunc(this_type *pThis) : m_pThis(pThis) {}
		template <typename TRefDim>
		void operator()(TRefDim &) { m_pThis->register_func_for_refDim<TRefDim::value>(); }
	};

	/// registers the update function for a reference dimension
	template <int refDim>
	void register_func_for_refDim()
	{
		//typedef DimFV1Geometry<refDim, dim> TGeom;
		//typedef bool (this_type::*TFunc) (const TGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
		//base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
	}
};

/// computes the closest node to a elem side ray intersection
template <typename TRefElem, int TWorldDim>
void GetNodeNextToCut(size_t &coOut,
					  const MathVector<TWorldDim> &IP,
					  const MathVector<TWorldDim> &IPVel,
					  const MathVector<TWorldDim> *vCornerCoords)
{
	//	help variables
	size_t side = 0;
	MathVector<TWorldDim> globalIntersection;
	MathVector<TRefElem::dim> localIntersection;

	//	compute intersection of ray in direction of ip velocity with elem side
	//	we search the ray only in upwind direction
	if (!ElementSideRayIntersection<TRefElem, TWorldDim>(side, globalIntersection, localIntersection,
														 IP, IPVel, false /* i.e. search upwind */, vCornerCoords))
		UG_THROW("GetNodeNextToCut: ElementSideRayIntersection failed");

	//	get reference element
	static const TRefElem &rRefElem = Provider<TRefElem>::get();
	const int dim = TRefElem::dim;

	// 	reset minimum
	number min = std::numeric_limits<number>::max();

	// 	loop corners of side
	for (size_t i = 0; i < rRefElem.num(dim - 1, side, 0); ++i)
	{
		// 	get corner
		const size_t co = rRefElem.id(dim - 1, side, 0, i);

		// 	Compute Distance to intersection
		number dist = VecDistanceSq(globalIntersection, vCornerCoords[co]);

		// 	if distance is smaller, choose this node
		if (dist < min)
		{
			min = dist;
			coOut = co;
		}
	}
}

template <int TDim>
template <typename TFVGeom>
bool ConvectionShapesSkewedUpwind<TDim>::
	update(const TFVGeom *geo,
		   const MathVector<dim> *Velocity,
		   const MathMatrix<dim, dim> *DiffDisp,
		   bool computeDeriv)
{
	UG_ASSERT(geo != nullptr, "Null pointer");
	UG_ASSERT(Velocity != nullptr, "Null pointer");

	//	\todo: think about: this should be something like scvf.num_sh()
	const size_t numSH = geo->num_sh();

	//	loop subcontrol volume faces
	for (size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
		//	get subcontrol volume face
		const typename TFVGeom::SCVF &scvf = geo->scvf(ip);

		//	Compute flux
		const number flux = VecDot(scvf.normal(), Velocity[ip]);

		// Switch to no upwind in case of small velocity
		if (VecTwoNorm(Velocity[ip]) < 1e-14)
		{
			// No upwind!
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				conv_shape(ip, sh) = flux * scvf.shape(sh);
			for (size_t sh = scvf.num_sh(); sh < numSH; sh++)
				conv_shape(ip, sh) = 0.0;
			if (computeDeriv)
			{
				for (size_t sh = 0; sh < scvf.num_sh(); sh++)
					VecScale(D_vel(ip, sh), scvf.normal(), scvf.shape(sh));
				// temporary: case of hanging nodes
				for (size_t sh = scvf.num_sh(); sh < numSH; sh++)
					VecSet(D_vel(ip, sh), 0.0);
			}
			continue;
		}

		//	initialize shapes with zeroes
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			conv_shape(ip, sh) = 0.0;

		//  Hanging nodes
		//	this is introduced here, hopefully temporarily: The problem is, that
		//	for hanging nodes the number of shape function is not the number of
		//	corners, but scvf.num_sh() currently returns the number of corners.
		//	this is actually enough to interpolate the function, but still we
		//	should reset the interpolation adding for hanging dofs to zero
		for (size_t sh = scvf.num_sh(); sh < numSH; sh++)
			conv_shape(ip, sh) = 0.0;

		const MathVector<dim> *vCornerCoords = geo->corners();
		size_t up = 0;
		try
		{
			GetNodeNextToCut<typename TFVGeom::ref_elem_type, dim>(up, scvf.global_ip(), Velocity[ip], vCornerCoords);
		}
		UG_CATCH_THROW("Skewed upwind: Cannot find node next to cut");

		conv_shape(ip, up) = flux;

		//	Write Derivatives if wanted
		if (computeDeriv)
		{
			for (size_t sh = 0; sh < numSH; ++sh)
				VecSet(D_vel(ip, sh), 0.0);
			D_vel(ip, up) = scvf.normal();
		}

		//	The shapes do not depend of the diffusion tensor
	}

	//	we're done
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class ConvectionShapesLinearProfileSkewedUpwind
	: public IConvectionShapes<TDim>
{
public:
	///	Base class
	typedef IConvectionShapes<TDim> base_type;

	///	This class
	typedef ConvectionShapesLinearProfileSkewedUpwind<TDim> this_type;

	///	Dimension
	static const int dim = TDim;

protected:
	//	explicitly forward some function
	using base_type::conv_shape;
	using base_type::conv_shape_diffusion;
	using base_type::D_vel;
	using base_type::non_zero_deriv_diffusion;
	using base_type::register_update_func;
	using base_type::set_non_zero_deriv_diffusion_flag;

public:
	///	constructor
	ConvectionShapesLinearProfileSkewedUpwind()
	{
		//	the shapes do not depend on the DiffDisp. Thus, we can set the
		//	derivative to be always zero w.r.t. the DiffDisp for all shapes
		set_non_zero_deriv_diffusion_flag(false);

		//	register evaluation function
		boost::mpl::for_each<typename domain_traits<dim>::AllElemList>(RegisterElemFunc(this));
		boost::mpl::for_each<boost::mpl::range_c<int, 1, dim + 1>>(RegisterRefDimFunc(this));
	}

	///	update of values for FV geometry
	template <typename TFVGeom>
	bool update(const TFVGeom *geo,
				const MathVector<dim> *Velocity,
				const MathMatrix<dim, dim> *DiffDisp,
				bool computeDeriv);

private:
	///	functor for registering the shapes for the element-templated FV geometries
	struct RegisterElemFunc
	{
		this_type *m_pThis;
		RegisterElemFunc(this_type *pThis) : m_pThis(pThis) {}
		template <typename TElem>
		void operator()(TElem &)
		{
			m_pThis->template register_func_for_elem_fvgeom<TElem, FV1Geometry<TElem, dim>>();
			//m_pThis->template register_func_for_elem_fvgeom< TElem, FV1CondensedGeometry<TElem, dim> >();
			//m_pThis->template register_func_for_elem_fvgeom< TElem, HFV1Geometry<TElem, dim> >();
		}
	};

	/// registers the update function for an element type and a FV geometry
	template <typename TElem, typename TFVGeom>
	void register_func_for_elem_fvgeom()
	{
		typedef bool (this_type::*TFunc)(const TFVGeom *, const MathVector<dim> *, const MathMatrix<dim, dim> *, bool);
		base_type::template register_update_func<TFVGeom, TFunc>(&this_type::template update<TFVGeom>);
	}

	///	functor for registering the shapes for the reference-dimension-templated FV geometries
	struct RegisterRefDimFunc
	{
		this_type *m_pThis;
		RegisterRefDimFunc(this_type *pThis) : m_pThis(pThis) {}
		template <typename TRefDim>
		void operator()(TRefDim &) { m_pThis->register_func_for_refDim<TRefDim::value>(); }
	};

	/// registers the update function for a reference dimension
	template <int refDim>
	void register_func_for_refDim()
	{
		//typedef DimFV1Geometry<refDim, dim> TGeom;
		//typedef bool (this_type::*TFunc) (const TGeom*, const MathVector<dim>*, const MathMatrix<dim, dim>*, bool);
		//base_type::template register_update_func<TGeom, TFunc>(&this_type::template update<TGeom>);
	}
};

template <int TDim>
template <typename TFVGeom>
bool ConvectionShapesLinearProfileSkewedUpwind<TDim>::
	update(const TFVGeom *geo,
		   const MathVector<dim> *Velocity,
		   const MathMatrix<dim, dim> *DiffDisp,
		   bool computeDeriv)
{
	UG_ASSERT(geo != nullptr, "Null pointer");
	UG_ASSERT(Velocity != nullptr, "Null pointer");

	//	\todo: think about: this should be something like scvf.num_sh()
	const size_t numSH = geo->num_sh();

	//	loop subcontrol volume faces
	for (size_t ip = 0; ip < geo->num_scvf(); ++ip)
	{
		//	get subcontrol volume face
		const typename TFVGeom::SCVF &scvf = geo->scvf(ip);

		//	Compute flux
		const number flux = VecDot(scvf.normal(), Velocity[ip]);

		// Switch to no upwind in case of small velocity
		if (VecTwoNorm(Velocity[ip]) < 1e-14)
		{
			// No upwind!
			for (size_t sh = 0; sh < scvf.num_sh(); sh++)
				conv_shape(ip, sh) = flux * scvf.shape(sh);
			for (size_t sh = scvf.num_sh(); sh < numSH; sh++)
				conv_shape(ip, sh) = 0.0;
			if (computeDeriv)
			{
				for (size_t sh = 0; sh < scvf.num_sh(); sh++)
					VecScale(D_vel(ip, sh), scvf.normal(), scvf.shape(sh));
				// temporary, see above
				for (size_t sh = scvf.num_sh(); sh < numSH; sh++)
					VecSet(D_vel(ip, sh), 0.0);
			}
			continue;
		}

		//	initialize shapes with zeroes
		for (size_t sh = 0; sh < scvf.num_sh(); sh++)
			conv_shape(ip, sh) = 0.0;

		//  Hanging nodes
		//	this is introduced here, hopefully temporarily: The problem is, that
		//	for hanging nodes the number of shape function is not the number of
		//	corners, but scvf.num_sh() currently returns the number of corners.
		//	this is actually enough to interpolate the function, but still we
		//	should reset the interpolation adding for hanging dofs to zero
		for (size_t sh = scvf.num_sh(); sh < numSH; sh++)
			conv_shape(ip, sh) = 0.0;

		const MathVector<dim> *vCornerCoords = geo->corners();
		// Reference element type
		typedef typename TFVGeom::ref_elem_type TRefElem;
		size_t side = 0;
		MathVector<dim> globalIntersection;
		MathVector<TRefElem::dim> localIntersection;

		try
		{
			ElementSideRayIntersection<TRefElem, dim>(side, globalIntersection, localIntersection,
													  scvf.global_ip(), Velocity[ip], false /* search upwind */, vCornerCoords);
		}
		UG_CATCH_THROW("GetLinearProfileSkewedUpwindShapes: Cannot find cut side.");

		// 	get linear trial space
		static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
		const LocalShapeFunctionSet<TRefElem::dim> &TrialSpace =
			LocalFiniteElementProvider::get<TRefElem::dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

		// 	get Reference Element
		static const TRefElem &rRefElem = Provider<TRefElem>::get();

        // Shape function values 
		std::vector<number> vShape;
		// Values of shape function gradients
		std::vector<MathVector<TRefElem::dim>> vShapeGrad;
		// Get values at the intersection point
		TrialSpace.shapes(vShape, localIntersection);
		TrialSpace.grads(vShapeGrad, localIntersection);

		size_t num_corners_inters_side = rRefElem.num(dim - 1, side, 0);

		// 	loop corners of side
		for (size_t j = 0; j < num_corners_inters_side; ++j)
		{
			// 	get corner
			const size_t co = rRefElem.id(dim-1, side, 0, j);

			//	write shape
			conv_shape(ip, co) = flux * vShape[co];
		}

		//	Write Derivatives if wanted
		if (computeDeriv)
		{
			// No derivatives yet!
		}

		//	The shapes do not depend of the diffusion tensor
	}

	//	we're done
	return true;
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE__ */
