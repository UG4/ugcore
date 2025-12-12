/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include "local_finite_element_provider.h"

#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"

// include spaces
#include "lagrange/lagrangep1.h"
#include "lagrange/lagrange.h"
#include "crouzeix-raviart/crouzeix_raviart.h"
#include "piecewise_constant/piecewise_constant.h"
#include "mini/mini.h"
#include "nedelec/nedelec.h"


namespace ug {

DebugID DID_LOCAL_FINITE_ELEMENT_PROVIDER("LOCAL_FINITE_ELEMENT_PROVIDER");

////////////////////////////////////////////////////////////////////////////////
//	Wrapper for LocalShapeFunctionSets
////////////////////////////////////////////////////////////////////////////////

/// wrapper class implementing the LocalShapeFunctionSet interface
/**
 * This class wrappes a class passed by the template argument into the
 * virtual ILocalShapeFunctionSet interface and makes it thus usable in that
 * context on the price of virtual functions.
 *
 * \tparam 	TImpl		Implementation of a Local Shape Function Set
 */
template <typename TImpl>
class LocalShapeFunctionSetWrapper
	: public LocalShapeFunctionSet<TImpl::dim,
	  	  	  	  	  	  	  	   typename TImpl::shape_type,
	  	  	  	  	  	  	  	   typename TImpl::grad_type>,
	  public TImpl
{
	/// Implementation
		using ImplType = TImpl;

	public:
	///	reference element dimension
		static constexpr int dim = TImpl::dim;

	///	Shape type
		using shape_type = typename ImplType::shape_type;

	///	Gradient type
		using grad_type = typename ImplType::grad_type;

	public:
	///	constructor
		LocalShapeFunctionSetWrapper(){}

	public:
	///	\copydoc ug::LocalDoFSet::type()
		virtual ReferenceObjectID roid() const {return ImplType::roid();}

	///	\copydoc ug::LocalDoFSet::num_sh()
		virtual size_t num_sh() const {return ImplType::num_sh();}

	///	\copydoc ug::LocalDoFSet::num_dof()
		virtual size_t num_dof(ReferenceObjectID roid) const {return ImplType::num_dof(roid);}

	///	\copydoc ug::LocalDoFSet::local_dof()
		virtual const LocalDoF& local_dof(size_t dof) const {return ImplType::local_dof(dof);}

	///	\copydoc ug::LocalDoFSet::exact_position_available()
		virtual bool exact_position_available() const {return ImplType::exact_position_available();}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		virtual bool position(size_t i, MathVector<dim>& pos) const
		{
			return ImplType::position(i, pos);
		}

	public:
	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		virtual bool continuous() const {return ImplType::continuous();}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return ImplType::shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual void shape(shape_type& s, size_t i, const MathVector<dim>& x) const
		{
			using baseType = BaseLSFS<TImpl, TImpl::dim,
				typename TImpl::shape_type,
				typename TImpl::grad_type>;

				  	  	  	  	  	  	  	   baseType::shape(s, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* vShape, const MathVector<dim>& x) const
		{
			ImplType::shapes(vShape, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(std::vector<shape_type>& vShape, const MathVector<dim>& x) const
		{
			ImplType::shapes(vShape, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(std::vector<std::vector<shape_type> >& vvShape,
							const std::vector<MathVector<dim> >& vLocPos) const
		{
			ImplType::shapes(vvShape, vLocPos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
			ImplType::grad(g, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* vGrad, const MathVector<dim>& x) const
		{
			ImplType::grads(vGrad, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(std::vector<grad_type>& vGrad, const MathVector<dim>& x) const
		{
			ImplType::grads(vGrad, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(std::vector<std::vector<grad_type> >& vvGrad,
						   const std::vector<MathVector<dim> >& vLocPos) const
		{
			ImplType::grads(vvGrad, vLocPos);
		}
};

////////////////////////////////////////////////////////////////////////////////
// 	SubLocalDoFSet
////////////////////////////////////////////////////////////////////////////////

/**
 * Intersection of local dofs
 */
template <int TDim>
class SubLocalDoFSet : public DimLocalDoFSet<TDim>
{
	static constexpr int dim = TDim;

	public:
		template <int setDim>
		SubLocalDoFSet(const ReferenceObjectID roid,
					   const DimLocalDoFSet<setDim>& set)
		   : m_roid(roid),
			 m_bExactPos(set.exact_position_available()),
			 m_bInit(false),
			 m_vNumDoF(NUM_REFERENCE_OBJECTS, 0)
		{
			if(ReferenceElementDimension(roid) != dim)
				UG_THROW("SubLocalDoFSet: templated Reference Dimension "<<dim<<
				         " does not match Reference Element dimension of "<<roid);

			const DimReferenceElement<setDim>& setRefElem =
					ReferenceElementProvider::get<setDim>(set.roid());

		// loop all subs of roid type
			for(size_t id = 0; id < setRefElem.num(dim); ++id){
				if(setRefElem.roid(dim, id) != m_roid) continue;

			// 	get mapping to sub
				std::vector<MathVector<setDim> > vCorner;
				for(size_t co = 0; co < setRefElem.num(dim, id, 0); ++co){
					vCorner.push_back(setRefElem.corner(setRefElem.id(dim, id, 0, co)));
				}
				const DimReferenceMapping<dim, setDim>& map =
							ReferenceMappingProvider::get<dim, setDim>(m_roid, vCorner);

				std::vector<size_t> vNumDoF(NUM_REFERENCE_OBJECTS, 0);
				std::vector<LocalDoF> vLocalDoF;
				std::vector<MathVector<dim> > vLocalPos;

			//	loop subselements of sub
				for(int subDim = 0; subDim <= dim; ++subDim){

					for(size_t i = 0; i < setRefElem.num(dim, id, subDim); ++i){
						const size_t subID = setRefElem.id(dim, id, subDim, i);

						const ReferenceObjectID sroid =  setRefElem.roid(subDim, subID);
						vNumDoF[sroid] = set.num_dof(sroid);

					//	get local DoFs on sub
						std::vector<size_t> vLocalDoFID;
						for(size_t dof = 0; dof < set.num_dof(); ++dof){
							const LocalDoF& localDoF = set.local_dof(dof);
							if(localDoF.dim() != subDim) continue;
							if(localDoF.id() != subID) continue;

							vLocalDoFID.push_back(dof);
						}

					//	loop sub DoFs in correct order
						for(size_t offset = 0; offset < vLocalDoFID.size(); ++offset){

							const LocalDoF* pLocalDoF = nullptr;
							size_t localDoFID = (size_t)-1;
							for(size_t dof = 0; dof < vLocalDoFID.size(); ++dof){
								if(set.local_dof(vLocalDoFID[dof]).offset() == offset){
									pLocalDoF = &set.local_dof(vLocalDoFID[dof]);
									localDoFID = vLocalDoFID[dof];
								}
							}

							if(pLocalDoF == nullptr)
								UG_THROW("SubLocalDoFSet: Cannot find local dof "
										"with offset "<<offset<<", but must "
										"exist. Check Implementation of LocalDoFSet ");

						//	map local dof
							MathVector<dim> locPos(0.0);
							MathVector<setDim> globPos;
							set.position(localDoFID, globPos);

							try{
								map.global_to_local(locPos, globPos);
							}
							catch(UGError& err){
								std::stringstream ss;
								ss<<"SubLocalDoFSet: Cannot find local position for "
								"global Position "<<globPos<<" of the "<<localDoFID<<
								"'th LocalDoF on "<<set.roid()<<
								" for the "<<id<<"'th "<<roid<<" and its "
								"Sub-"<<sroid<<" with id "<<subID<<" with "
								"corners:\n";
								for(size_t i = 0; i < vCorner.size(); ++i)
									ss << " "<<i<<": "<<vCorner[i] << "\n";

								err.push_msg(ss.str(),__FILE__,__LINE__);
								throw(err);
							}

						//	add
							vLocalDoF.push_back(LocalDoF(subDim, i, offset));
							vLocalPos.push_back(locPos);
						}
					}
				}

				this->set(vNumDoF, vLocalDoF, vLocalPos);
			}
		}

	public:
		virtual ReferenceObjectID roid() const {return m_roid;};
		virtual size_t num_dof(ReferenceObjectID roid) const {return m_vNumDoF[roid];}
		virtual size_t num_sh() const {return m_vLocalDoF.size();}
		virtual const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}
		virtual bool exact_position_available() const {return m_bExactPos;}
		virtual bool position(size_t i, MathVector<TDim>& pos) const
		{
			pos = m_vLocalPos[i];
			return exact_position_available();
		}

	protected:
		void set(const std::vector<size_t>& vNumDoF,
				 const std::vector<LocalDoF>& vLocalDoF,
				 const std::vector<MathVector<TDim> >& vLocalPos)
		{
			if(m_bInit){
				if(vNumDoF != m_vNumDoF) UG_THROW("NumDoF mismatch");
				if(vLocalDoF != m_vLocalDoF) UG_THROW("vLocalDoF mismatch");
				if(vLocalPos != m_vLocalPos) UG_THROW("vLocalPos mismatch");
			}
			else{
				m_vNumDoF = vNumDoF;
				m_vLocalDoF = vLocalDoF;
				m_vLocalPos = vLocalPos;
			}
		}


	protected:
		const ReferenceObjectID m_roid; ///< Reference ID this DoF Set is for
		const bool m_bExactPos;

		bool m_bInit;
		std::vector<size_t> m_vNumDoF;
		std::vector<LocalDoF> m_vLocalDoF; ///< Local DoFs of this set
		std::vector<MathVector<TDim> > m_vLocalPos; ///< Local Positions of DoFs
};

////////////////////////////////////////////////////////////////////////////////
//	Provider for LocalShapeFunctionSets
////////////////////////////////////////////////////////////////////////////////

template <typename TRefElem>
void LocalFiniteElementProvider::create_lagrange_set(const LFEID& id)
{
//	reference object id
	static constexpr ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static constexpr int dim = TRefElem::dim;

//	only order >= 1 available
	if(id.order() < 1) return;

//	if refdim == dim create the space
	if(dim == id.dim()){
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set;
		switch(id.order()){
			case 1: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem> >);
					break;
			case 2: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem,2> >);
					break;
			break;
			case 3: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem,3> >);
					break;
			break;
			case 4: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						  (new LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem,4> >);
					break;
			break;

			default:{
				SmartPtr<LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> > >
					sSetFlexLagrange(new LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >);
				sSetFlexLagrange->set_order(id.order());
				set = sSetFlexLagrange;
			}
		}
		register_set(id, set);
		return;
	}
//	if refdim < dim, the restriction of to the subelement is the lagrange space
//	of the refdim
	else if (dim < id.dim()){
	//  get (and maybe create) lagrange space for dim == refdim (always exist)
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set =
			getptr<dim>(roid, LFEID(id.type(), dim, id.order()));

	//	register this space as space on subelement
		register_set(id, set);
	}

	// else: for refdim > dim there exist no restrictions
}

void LocalFiniteElementProvider::
create_lagrange_set(ReferenceObjectID roid, const LFEID& id)
{
	UG_DLOG(DID_LOCAL_FINITE_ELEMENT_PROVIDER, 2, ">>OCT_DISC_DEBUG: " << "local_finite_element_provider.cpp: " << "create_lagrange_set(): " << roid << std::endl);
	switch(roid){
		case ROID_VERTEX:		create_lagrange_set<ReferenceVertex>(id); return;
		case ROID_EDGE:			create_lagrange_set<ReferenceEdge>(id); return;
		case ROID_TRIANGLE:		create_lagrange_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_lagrange_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_lagrange_set<ReferenceTetrahedron>(id); return;
		case ROID_PRISM:		create_lagrange_set<ReferencePrism>(id); return;
		case ROID_HEXAHEDRON:	create_lagrange_set<ReferenceHexahedron>(id); return;
		case ROID_PYRAMID:
			// only space available for order 1
			if(id.order() != 1)
				return;
			else
			{
				register_set(LFEID(LFEID::LAGRANGE, 3, 1),
							 ConstSmartPtr<LocalShapeFunctionSet<3> >
							(new LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid> >));
			}
			return;
		case ROID_OCTAHEDRON:
			// only space available for order 1
			if(id.order() != 1)
				return;
			else
			{
				register_set(LFEID(LFEID::LAGRANGE, 3, 1),
							 ConstSmartPtr<LocalShapeFunctionSet<3> >
							(new LocalShapeFunctionSetWrapper<LagrangeP1<ReferenceOctahedron> >));
			}
			return;
		default: return;
	}
}

template <typename TRefElem>
void LocalFiniteElementProvider::create_mini_bubble_set(const LFEID& id)
{
//	reference object id
	static constexpr ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static constexpr int dim = TRefElem::dim;

//	if refdim == dim create the space
	if(dim == id.dim()){

		if(id == LFEID(LFEID::MINI, dim, 1))
			register_set(id,  ConstSmartPtr<LocalShapeFunctionSet<dim> >(
						 new LocalShapeFunctionSetWrapper<MiniBubbleLSFS<TRefElem> >));
		return;

		/*
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set;

		switch(id.order()){
			case 1: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						(new LocalShapeFunctionSetWrapper<MiniBubbleLSFS<TRefElem> >);
					break;
			case 2: set = ConstSmartPtr<LocalShapeFunctionSet<dim> >
						(new LocalShapeFunctionSetWrapper<MiniBubbleLSFS<TRefElem,2> >);
					break;

		}*/

	}
//	if refdim < dim, the restriction of to the subelement is the lagrange space
//	of the refdim
	else if (dim < id.dim()){
	//  get (and maybe create) lagrange space for dim == refdim (always exist)
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set =
			getptr<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

	//	register this space as space on subelement
		register_set(id, set);
	}

	// else: for refdim > dim there exist no restrictions
}

void LocalFiniteElementProvider::
create_mini_bubble_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_EDGE:			create_mini_bubble_set<ReferenceEdge>(id); return;
		case ROID_TRIANGLE:		create_mini_bubble_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_mini_bubble_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_mini_bubble_set<ReferenceTetrahedron>(id); return;
		// no spaces implemented for other 3d elems
		default: return;
	}
}

template <typename TRefElem>
void LocalFiniteElementProvider::create_nedelec_set(const LFEID& id)
{
	static constexpr int dim = TRefElem::dim;
	if(id == LFEID(LFEID::NEDELEC, dim, 1))
		register_set(id, ConstSmartPtr<LocalShapeFunctionSet<dim, MathVector<dim>, MathMatrix<dim,dim> > >(
					 new LocalShapeFunctionSetWrapper<NedelecLSFS<TRefElem> >));
}

void LocalFiniteElementProvider::
create_nedelec_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_TRIANGLE:		create_nedelec_set<ReferenceTriangle>(id); return;
		case ROID_TETRAHEDRON:	create_nedelec_set<ReferenceTetrahedron>(id); return;
		default: return;
	}
}

template <typename TRefElem>
void LocalFiniteElementProvider::create_piecewise_constant_set(const LFEID& id)
{
	static constexpr int dim = TRefElem::dim;
	if(id == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
		register_set(id, ConstSmartPtr<LocalShapeFunctionSet<dim> >(
					 new LocalShapeFunctionSetWrapper<PiecewiseConstantLSFS<TRefElem> >));
}

void LocalFiniteElementProvider::
create_piecewise_constant_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		case ROID_EDGE:			create_piecewise_constant_set<ReferenceEdge>(id); return;
		case ROID_TRIANGLE:		create_piecewise_constant_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_piecewise_constant_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_piecewise_constant_set<ReferenceTetrahedron>(id); return;
		case ROID_PRISM:		create_piecewise_constant_set<ReferencePrism>(id); return;
		case ROID_PYRAMID:		create_piecewise_constant_set<ReferencePyramid>(id); return;
		case ROID_HEXAHEDRON:	create_piecewise_constant_set<ReferenceHexahedron>(id); return;
		case ROID_OCTAHEDRON:	create_piecewise_constant_set<ReferenceOctahedron>(id); return;
		default: return;
	}
}


template <typename TRefElem>
void LocalFiniteElementProvider::create_crouxeiz_raviart_set(const LFEID& id)
{
	static constexpr int dim = TRefElem::dim;
	if(id == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
		register_set(id, ConstSmartPtr<LocalShapeFunctionSet<dim> >(
					 new LocalShapeFunctionSetWrapper<CrouzeixRaviartLSFS<TRefElem> >));
	return;
}

void LocalFiniteElementProvider::
create_crouxeiz_raviart_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(roid){
		// no space for dim <= 1
		case ROID_TRIANGLE:		create_crouxeiz_raviart_set<ReferenceTriangle>(id); return;
		case ROID_QUADRILATERAL:create_crouxeiz_raviart_set<ReferenceQuadrilateral>(id); return;
		case ROID_TETRAHEDRON:	create_crouxeiz_raviart_set<ReferenceTetrahedron>(id); return;
		case ROID_PRISM:		create_crouxeiz_raviart_set<ReferencePrism>(id); return;
		case ROID_PYRAMID:		create_crouxeiz_raviart_set<ReferencePyramid>(id); return;
		case ROID_HEXAHEDRON:	create_crouxeiz_raviart_set<ReferenceHexahedron>(id); return;
		default: return;
	}
}

void LocalFiniteElementProvider::
create_set(ReferenceObjectID roid, const LFEID& id)
{
	try{
		switch(id.type())
		{
			case LFEID::LAGRANGE: 			create_lagrange_set(roid, id); break;
			case LFEID::MINI:   			create_mini_bubble_set(roid, id); break;
			case LFEID::NEDELEC:			create_nedelec_set(roid, id); break;
			case LFEID::PIECEWISE_CONSTANT:	create_piecewise_constant_set(roid, id); break;
			case LFEID::CROUZEIX_RAVIART:	create_crouxeiz_raviart_set(roid, id); break;
			default: return;
		}
	}
	UG_CATCH_THROW("LocalFiniteElementProvider: Creation of set "<<id<<
					" for "<<roid<<" failed.");
}

void LocalFiniteElementProvider::
create_set(const LFEID& id)
{
	switch(id.type())
	{
	//	this spaces can create also sub-elem spaces, since they are continuous
		case LFEID::LAGRANGE:
		case LFEID::MINI:

			for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
				const ReferenceObjectID roid = (ReferenceObjectID) r;
				const int dim = ReferenceElementDimension(roid);

				if(dim <= id.dim())
					create_set(roid, id);
			}
			break;

	//	this spaces can not create also sub-elem spaces, since they are discontinuous
		case LFEID::PIECEWISE_CONSTANT:
		case LFEID::CROUZEIX_RAVIART:
		case LFEID::NEDELEC:

			for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
				const ReferenceObjectID roid = (ReferenceObjectID) r;
				const int dim = ReferenceElementDimension(roid);

				if(dim == id.dim())
					create_set((ReferenceObjectID)roid, id);
			}
			break;

		default: return;
	}
}

template <int rdim, int dim>
void LocalFiniteElementProvider::
create_sub_dof_set(ReferenceObjectID roid, const LFEID& id)
{
	if(dim != id.dim())
		UG_THROW("Dimension must match here, internal error. ("<<dim<<","<<id.dim()<<")");

	// we like to have a DimLocalDoFSet for roid in dimension dim.
	// Say roid has dimension rdim. If rdim == dim this must be implemented.
	// If rdim < dim, this can be deduced from the implemented patterns for
	// elements in dim if and only if some implementation is for a element that
	// contains roid as a sub element

	for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
		const ReferenceObjectID elemRoid = (ReferenceObjectID) r;
		const int elemDim = ReferenceElementDimension(elemRoid);

		// we only take elements that are in the dimension of the space
		if(elemDim != dim) continue;

		// try to get the implementation, if not present skip
		ConstSmartPtr<DimLocalDoFSet<dim> > set = get_dof_ptr<dim>(elemRoid, id);
		if(set.invalid()) continue;

		// see if element contains the roid
		const ReferenceElement& rRefElem = ReferenceElementProvider::get(elemRoid);
		if(rRefElem.num(roid) == 0) continue;

		try{
			// create the sub-dof-pattern
			ConstSmartPtr<DimLocalDoFSet<rdim> > subSet =
				ConstSmartPtr<DimLocalDoFSet<rdim> >(new SubLocalDoFSet<rdim>(roid, *set) );

			// try to get an already registerd one
			ConstSmartPtr<DimLocalDoFSet<rdim> > givenSubSet = get_dof_ptr<rdim>(roid, id, false);

			// if already one set given, check equality
			if(givenSubSet.valid()){
				if(*givenSubSet != *subSet)
					UG_THROW("LocalFiniteElementProvider::create_sub_dof_set:\n"
							"Creating DimLocalDoFSet for "<<roid<<" and "<<id<<
							".\nAlready registered Set does not match with computed,"
							" but this indicates that the \nLocalDoFSets have not "
							"been implemented correctly. Check implementation.\n"
							"Sets are: \nGiven:\n"<<*givenSubSet<<"\nvs. New:\n"<<*subSet);
			}
			else {
				// if correct, register set
				register_set<rdim>(id, subSet);
			}
		}
		UG_CATCH_THROW("LocalFiniteElementProvider::create_sub_dof_set: Cannot "
						"create SubDoFSet for "<<roid<<" and space "<<id<<" when "
						"trying to deduce from LocalDoFSet for "<<elemRoid<<".");
	}
}

void LocalFiniteElementProvider::
create_dof_set(ReferenceObjectID roid, const LFEID& id)
{
	const int dim = id.dim();
	const int rdim = ReferenceElementDimension(roid);

	if(rdim == dim) {
		create_set(roid, id);
		return;
	}

	if(rdim > dim)
		UG_THROW("Wrong dimension for SubDoFs: "<<rdim<<", "<<dim);

	switch(dim){
		case 1:
			switch(rdim){
				case 0: create_sub_dof_set<0, 1>(roid, id); return;
				default: return;
			}
		case 2:
			switch(rdim){
				case 0: create_sub_dof_set<0, 2>(roid, id); return;
				case 1: create_sub_dof_set<1, 2>(roid, id); return;
				default: return;
			}
		case 3:
			switch(rdim){
				case 0: create_sub_dof_set<0, 3>(roid, id); return;
				case 1: create_sub_dof_set<1, 3>(roid, id); return;
				case 2: create_sub_dof_set<2, 3>(roid, id); return;
				default: return;
			}
		default: return;
	}
}


LocalFiniteElementProvider::
LocalFiniteElementProvider()
{};

LocalFiniteElementProvider::
~LocalFiniteElementProvider()
{};

const LocalDoFSet& LocalFiniteElementProvider::
get_dofs(ReferenceObjectID roid, const LFEID& id, bool bCreate)
{
//	init provider and search for identifier
	using Map = std::map<LFEID, LocalDoFSets>;
	Map::const_iterator iter = inst().m_mLocalDoFSets.find(id);

//	if not found
	if(iter == m_mLocalDoFSets.end() || (iter->second)[roid].invalid()){
		if(bCreate){
			create_dof_set(roid, id);
			return get_dofs(roid, id, false);
		}
		UG_THROW("LocalFiniteElementProvider: Cannot create LocalDoFSet for finite "
				"element type "<<id<<" and "<<roid);
	}

//	return dof set
	return *((iter->second)[roid]);
}

const CommonLocalDoFSet& LocalFiniteElementProvider::
get_dofs(const LFEID& id, bool bCreate)
{
//	init provider and search for identifier
	using Map = std::map<LFEID, CommonLocalDoFSet>;
	Map::const_iterator iter = inst().m_mCommonDoFSet.find(id);

//	if not found
	if(iter == m_mCommonDoFSet.end()){
		if(bCreate)	{
			create_set(id);
			return get_dofs(id, false);
		}
		UG_THROW("LocalFiniteElementProvider: Cannot create CommonLocalDoFSet for "<<id);
	}

//	return the common set
	return iter->second;
}

bool LocalFiniteElementProvider::continuous(const LFEID& id, bool bCreate)
{
	std::map<LFEID, bool>::iterator iter = m_mContSpace.find(id);
	if(iter == m_mContSpace.end())
	{
	//	try to create the set
		if(bCreate){
			create_set(id);
			return continuous(id, false);
		}

		UG_THROW("LocalFiniteElementProvider::continuous: no shape function "
				"set "<<id<<" registered.");
	}

	return (*iter).second;
}

void LocalFiniteElementProvider::register_set(const LFEID& id, ConstSmartPtr<LocalDoFSet> set)
{
//	reference object id
	const ReferenceObjectID roid = set->roid();

//	register
	m_mLocalDoFSets[id][roid] = set;

//	for creation of CommonLocalDoFSet: skip if not the dimension of the space
	if(set->dim() != id.dim()) return;

//	add this local dof set
	try{
		m_mCommonDoFSet[id].add(*set);
	}
	catch(UGError& err)
	{
	//	write error message
		std::stringstream ss;
		ss<<"LocalFiniteElementProvider::register_set(): "
				"Cannot build CommonLocalDoFSet for type: "<<id<<" when adding "
				" Reference element type "<<roid<<".\n"<<
				"CommonLocalDoFSet is:\n" << m_mCommonDoFSet[id]<<
				"LocalDoFSet is:\n" << *set;
		err.push_msg(ss.str(),__FILE__,__LINE__);

	//	remove set
		m_mCommonDoFSet.erase(id);

		throw(err);
	}
}

std::map<LFEID, LocalFiniteElementProvider::LocalDoFSets>
LocalFiniteElementProvider::m_mLocalDoFSets =
std::map<LFEID, LocalFiniteElementProvider::LocalDoFSets>();

std::map<LFEID, CommonLocalDoFSet>
LocalFiniteElementProvider::m_mCommonDoFSet =
std::map<LFEID, CommonLocalDoFSet>();

std::map<LFEID, bool>
LocalFiniteElementProvider::m_mContSpace =
std::map<LFEID, bool>();

} // namespace ug

