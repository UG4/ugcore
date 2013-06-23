/*
 * local_shape_function_set.cpp
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#include "local_shape_function_set.h"

// include spaces
#include "lagrange/lagrangep1.h"
#include "lagrange/lagrange.h"
#include "crouzeix-raviart/crouzeix_raviart.h"
#include "piecewise_constant/piecewise_constant.h"
#include "mini/mini.h"
#include "nedelec/nedelec.h"


namespace ug{

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
		typedef TImpl ImplType;

	public:
	///	Domain position type
		typedef typename ImplType::position_type position_type;

	///	Shape type
		typedef typename ImplType::shape_type shape_type;

	///	Gradient type
		typedef typename ImplType::grad_type grad_type;

	public:
	///	constructor
		LocalShapeFunctionSetWrapper(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		virtual LFEID type() const {return ImplType::type();}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		virtual bool continuous() const {return ImplType::continuous();}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		virtual size_t num_sh() const {return ImplType::num_sh();}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		virtual bool position(size_t i, position_type& pos) const
		{
			return ImplType::position(i, pos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual shape_type shape(size_t i, const position_type& x) const
		{
			return ImplType::shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual void shape(shape_type& s, size_t i, const position_type& x) const
		{
			typedef BaseLocalShapeFunctionSet<TImpl, TImpl::dim,
				  	  	  	  	  	  	  	   typename TImpl::shape_type,
				  	  	  	  	  	  	  	   typename TImpl::grad_type> baseType;

				  	  	  	  	  	  	  	   baseType::shape(s, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* vShape, const position_type& x) const
		{
			ImplType::shapes(vShape, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(std::vector<shape_type>& vShape, const position_type& x) const
		{
			ImplType::shapes(vShape, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(std::vector<std::vector<shape_type> >& vvShape,
							const std::vector<position_type>& vLocPos) const
		{
			ImplType::shapes(vvShape, vLocPos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual void grad(grad_type& g, size_t i, const position_type& x) const
		{
			ImplType::grad(g, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* vGrad, const position_type& x) const
		{
			ImplType::grads(vGrad, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(std::vector<grad_type>& vGrad, const position_type& x) const
		{
			ImplType::grads(vGrad, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(std::vector<std::vector<grad_type> >& vvGrad,
						   const std::vector<position_type>& vLocPos) const
		{
			ImplType::grads(vvGrad, vLocPos);
		}
};

////////////////////////////////////////////////////////////////////////////////
//	Provider for LocalShapeFunctionSets
////////////////////////////////////////////////////////////////////////////////

template <typename TRefElem>
void LocalShapeFunctionSetProvider::create_lagrange_set(const LFEID& id)
{
//	reference object id
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

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
		register_set(id, roid, set);
		return;
	}
//	if refdim < dim, the restriction of to the subelement is the lagrange space
//	of the refdim
	else if (dim < id.dim()){
	//  get (and maybe create) lagrange space for dim == refdim (always exist)
		ConstSmartPtr<LocalShapeFunctionSet<dim> > set =
			getptr<dim>(roid, LFEID(id.type(), dim, id.order()));

	//	register this space as space on subelement
		register_set(id, roid, set);
	}

	// else: for refdim > dim there exist no restrictions
}

void LocalShapeFunctionSetProvider::
create_lagrange_set(ReferenceObjectID roid, const LFEID& id)
{
//	only order >= 1 available
	if(id.order() < 1) return;

	try{
	//	switch type
		switch(roid)
		{
			case ROID_EDGE:			create_lagrange_set<ReferenceEdge>(id); return;
			case ROID_TRIANGLE:		create_lagrange_set<ReferenceTriangle>(id); return;
			case ROID_QUADRILATERAL:create_lagrange_set<ReferenceQuadrilateral>(id); return;
			case ROID_TETRAHEDRON:	create_lagrange_set<ReferenceTetrahedron>(id); return;
			case ROID_PRISM:		create_lagrange_set<ReferencePrism>(id); return;
			case ROID_HEXAHEDRON:	create_lagrange_set<ReferenceHexahedron>(id); return;

			case ROID_PYRAMID:
				// only space available for order 1
				if(id.order() != 1) return;
				{
					register_set(LFEID(LFEID::LAGRANGE, 3, 1), ROID_PYRAMID,
					             ConstSmartPtr<LocalShapeFunctionSet<3> >
								(new LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid> >));
				}

			default: return;
		}
	}
	UG_CATCH_THROW("LocalShapeFunctionSetProvider: Creation of Lagrange set "
					<<id<<" for "<<roid<<" failed.");
}


template <typename TRefElem>
void LocalShapeFunctionSetProvider::create_generic_sets(const LFEID& id)
{
//	reference object id
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;
	static const int dim = TRefElem::dim;

	switch(id.type()){
		case LFEID::PIECEWISE_CONSTANT:
			if(id == LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0))
				register_set(id, roid,
							 ConstSmartPtr<LocalShapeFunctionSet<dim> >(
							 new LocalShapeFunctionSetWrapper<PiecewiseConstantLSFS<TRefElem> >));
			return;
		case LFEID::CROUZEIX_RAVIART:
			if(id == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
				register_set(id, roid,
							 ConstSmartPtr<LocalShapeFunctionSet<dim> >(
							 new LocalShapeFunctionSetWrapper<CrouzeixRaviartLSFS<TRefElem> >));
			return;
		case LFEID::NEDELEC:
			if(id == LFEID(LFEID::NEDELEC, dim, 1))
				register_set(id, roid,
							 ConstSmartPtr<LocalShapeFunctionSet<dim, MathVector<dim>, MathMatrix<dim,dim> > >(
							 new LocalShapeFunctionSetWrapper<NedelecLSFS<TRefElem> >));
			return;

		case LFEID::MINI:
			// \todo: P1 bubble spaces restricted to a sub-geometric object (e.g. side)
			//		  are the usual P1 lagrange spaces. Add them here.
			if(id == LFEID(LFEID::MINI, dim, 1))
				register_set(id, roid,
							 ConstSmartPtr<LocalShapeFunctionSet<dim> >(
							 new LocalShapeFunctionSetWrapper<MiniBubbleLSFS<TRefElem> >));
			return;
		default: return;
	}
}

void LocalShapeFunctionSetProvider::
create_generic_sets(ReferenceObjectID roid, const LFEID& id)
{
	try{
	//	switch type
		switch(roid)
		{
			case ROID_EDGE:			create_generic_sets<ReferenceEdge>(id); return;
			case ROID_TRIANGLE:		create_generic_sets<ReferenceTriangle>(id); return;
			case ROID_QUADRILATERAL:create_generic_sets<ReferenceQuadrilateral>(id); return;
			case ROID_TETRAHEDRON:	create_generic_sets<ReferenceTetrahedron>(id); return;
			case ROID_PRISM:		create_generic_sets<ReferencePrism>(id); return;
			case ROID_PYRAMID:		create_generic_sets<ReferencePyramid>(id); return;
			case ROID_HEXAHEDRON:	create_generic_sets<ReferenceHexahedron>(id); return;
			default: return;
		}
	}
	UG_CATCH_THROW("LocalShapeFunctionSetProvider: Creation of set "<<id<<
					" for "<<roid<<" failed.");
}

void LocalShapeFunctionSetProvider::
create_set(ReferenceObjectID roid, const LFEID& id)
{
	switch(id.type())
	{
		case LFEID::LAGRANGE:
			create_lagrange_set(roid, id);
			break;

		case LFEID::PIECEWISE_CONSTANT:
		case LFEID::CROUZEIX_RAVIART:
		case LFEID::MINI:
		case LFEID::NEDELEC:
			create_generic_sets(roid, id);
			break;

		default: return;
	}
}

void LocalShapeFunctionSetProvider::
create_set(const LFEID& id)
{
	for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
		create_set((ReferenceObjectID)roid, id);
}


LocalShapeFunctionSetProvider::
LocalShapeFunctionSetProvider()
{};

LocalShapeFunctionSetProvider::
~LocalShapeFunctionSetProvider()
{};

bool LocalShapeFunctionSetProvider::continuous(const LFEID& type, bool bCreate)
{
	std::map<LFEID, bool>::iterator iter = m_mContSpace.find(type);
	if(iter == m_mContSpace.end())
	{
		if(bCreate)
		{
		//	try to create the set
			create_set(type);

		//	next try to return the set
			return continuous(type, false);
		}

		UG_THROW("LocalShapeFunctionSetProvider::continuous: no shape function "
				"set "<<type<<" registered.");
	}

	return (*iter).second;
}

std::map<LFEID, bool>
LocalShapeFunctionSetProvider::m_mContSpace = std::map<LFEID, bool>();

} // namespace ug

