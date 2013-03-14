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

std::map<LFEID, bool>&
LocalShapeFunctionSetProvider::get_continuous_map()
{
//	get type of map
	typedef std::map<LFEID, bool> Vec;

//	create static map
	static Vec sShapeFunctionSetMap;

//	return map
	return sShapeFunctionSetMap;
};

bool LocalShapeFunctionSetProvider::continuous(const LFEID& type, bool bCreate)
{
	std::map<LFEID, bool>& contMap = inst().get_continuous_map();
	std::map<LFEID, bool>::iterator iter = contMap.find(type);
	if(iter == contMap.end())
	{
		if(bCreate)
		{
		//	try to create the set
			dynamically_create_set(type);

		//	next try to return the set
			return continuous(type, false);
		}

		UG_THROW("LocalShapeFunctionSetProvider::continuous: no shape function "
				"set "<<type<<" registered.");
	}

	return (*iter).second;
}

template <typename TRefElem>
void LocalShapeFunctionSetProvider::init_standard_sets()
{
//	reference object id
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem> > sSetLagrangeP1;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 2> > sSetLagrangeP2;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 3> > sSetLagrangeP3;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 4> > sSetLagrangeP4;

//	insert into map: P1 Lagrange
	LFEID type1(LFEID::LAGRANGE, 1);
	register_set(type1, roid, sSetLagrangeP1);

//	insert into map: P2 Lagrange
	LFEID type2(LFEID::LAGRANGE, 2);
	register_set(type2, roid, sSetLagrangeP2);

//	insert into map: P3 Lagrange
	LFEID type3(LFEID::LAGRANGE, 3);
	register_set(type3, roid, sSetLagrangeP3);

//	insert into map: P4 Lagrange
	LFEID type4(LFEID::LAGRANGE, 4);
	register_set(type4, roid, sSetLagrangeP4);


//	insert into map: Crouzeix-Raviart
	static LocalShapeFunctionSetWrapper<CrouzeixRaviartLSFS<TRefElem> > sSetCrouzeixRaviart;
	LFEID typeCR(LFEID::CROUZEIX_RAVIART, 1);
	register_set(typeCR, roid, sSetCrouzeixRaviart);

//	insert into map: Piecewise constant
	static LocalShapeFunctionSetWrapper<PiecewiseConstantLSFS<TRefElem> > sSetPiecewiseConstant;
	LFEID typePC(LFEID::PIECEWISE_CONSTANT, 0);
	register_set(typePC, roid, sSetPiecewiseConstant);

//	insert into map: MINI element
	static LocalShapeFunctionSetWrapper<MiniBubbleLSFS<TRefElem> > sSetMiniBubble;
	LFEID typeMB(LFEID::MINI, 1);
	register_set(typeMB, roid, sSetMiniBubble);

//	insert into map: Nedelec element
	static LocalShapeFunctionSetWrapper<NedelecLSFS<TRefElem> > sSetNedelec;
	LFEID typeNedelec(LFEID::NEDELEC, 1);
	register_set(typeNedelec, roid, sSetNedelec);

}

template <typename TRefElem>
void LocalShapeFunctionSetProvider::init_flex_lagrange(size_t order)
{
//	reference object id
	static const ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	create static Sets
	LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >* sSetFlexLagrange
		= new LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >;
	sSetFlexLagrange->set_order(order);

//	insert into map: Lagrange
	LFEID type(LFEID::LAGRANGE, order);
	register_set(type, roid, *sSetFlexLagrange);

//	remember creation
	get_dynamic_allocated_vector<TRefElem::dim>().push_back(sSetFlexLagrange);
}

void LocalShapeFunctionSetProvider::
dynamically_create_set(ReferenceObjectID roid, LFEID id)
{
//	Lagrange space
	if(id.type() == LFEID::LAGRANGE)
	{
	//	only order >= 1 available
		if(id.order() < 1) return;

		try{
	//	switch type
		switch(roid)
		{
			case ROID_EDGE:			init_flex_lagrange<ReferenceEdge>(id.order()); return;
			case ROID_TRIANGLE:		init_flex_lagrange<ReferenceTriangle>(id.order()); return;
			case ROID_QUADRILATERAL:init_flex_lagrange<ReferenceQuadrilateral>(id.order()); return;
			case ROID_TETRAHEDRON:	init_flex_lagrange<ReferenceTetrahedron>(id.order()); return;
			case ROID_PRISM:		init_flex_lagrange<ReferencePrism>(id.order()); return;
			case ROID_HEXAHEDRON:	init_flex_lagrange<ReferenceHexahedron>(id.order()); return;
			default: UG_THROW("LocalShapeFunctionSetProvider: Roid="<<roid<<" unknown.");
		}

		}
		UG_CATCH_THROW("Dynamic Allocation of set failed.");
	}
}

void LocalShapeFunctionSetProvider::
dynamically_create_set(LFEID id)
{
//	Lagrange space
	if(id.type() == LFEID::LAGRANGE)
	{
	//	only order >= 1 available
		if(id.order() < 1) return;

		try{
	//	switch type
		init_flex_lagrange<ReferenceEdge>(id.order()); return;
		init_flex_lagrange<ReferenceTriangle>(id.order()); return;
		init_flex_lagrange<ReferenceQuadrilateral>(id.order()); return;
		init_flex_lagrange<ReferenceTetrahedron>(id.order()); return;
		init_flex_lagrange<ReferencePrism>(id.order()); return;
		init_flex_lagrange<ReferenceHexahedron>(id.order()); return;
		}
		UG_CATCH_THROW("Dynamic Allocation of set failed.");
	}
}

LocalShapeFunctionSetProvider::
LocalShapeFunctionSetProvider()
{
	static bool init = false;

	if(!init)
	{
	//	remember initialization
		init = true;

	//	clear created holder
		get_dynamic_allocated_vector<1>().clear();
		get_dynamic_allocated_vector<2>().clear();
		get_dynamic_allocated_vector<3>().clear();

	//	register all element types that allow higher orders
		try{
			init_standard_sets<ReferenceEdge>();
		}UG_CATCH_THROW("Cannot register standard Edge trial spaces.");
		try{
			init_standard_sets<ReferenceTriangle>();
		}UG_CATCH_THROW("Cannot register standard Triangle trial spaces.");
		try{
			init_standard_sets<ReferenceQuadrilateral>();
		}UG_CATCH_THROW("Cannot register standard Quadrilateral trial spaces.");
		try{
			init_standard_sets<ReferenceTetrahedron>();
		}UG_CATCH_THROW("Cannot register standard Tetrahedron trial spaces.");
		try{
			init_standard_sets<ReferencePrism>();
		}UG_CATCH_THROW("Cannot register standard Prism trial spaces.");
		try{
			init_standard_sets<ReferenceHexahedron>();
		}UG_CATCH_THROW("Cannot register standard Hexahedron trial spaces.");

	//	register 1st order pyramid
		LFEID type1(LFEID::LAGRANGE, 1);
		static LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid> > sSetLagrangeP1;
		try{
			register_set(type1, ROID_PYRAMID, sSetLagrangeP1);
		}UG_CATCH_THROW("Cannot register Pyramid P1 Lagrange trial spaces.");
		
	//	register pyramid for Crouzeix-Raviart element
		LFEID type2(LFEID::CROUZEIX_RAVIART, 1);
		static LocalShapeFunctionSetWrapper<CrouzeixRaviartLSFS<ReferencePyramid> > sSetCrouzeixRaviart;
		try{
			register_set(type2, ROID_PYRAMID, sSetCrouzeixRaviart);
		}UG_CATCH_THROW("Cannot register Pyramid Crouzeix-Raviart trial spaces.");	
		
	//	register pyramid for piecewise constant element
		LFEID type3(LFEID::PIECEWISE_CONSTANT, 0);
		static LocalShapeFunctionSetWrapper<PiecewiseConstantLSFS<ReferencePyramid> > sSetPiecewiseConstant;
		try{
			register_set(type3, ROID_PYRAMID, sSetPiecewiseConstant);
		}UG_CATCH_THROW("Cannot register Pyramid piecewise constant trial spaces.");

	}
};

LocalShapeFunctionSetProvider::
~LocalShapeFunctionSetProvider()
{
	std::vector<LocalShapeFunctionSet<1>*>& dynVec1 =
											get_dynamic_allocated_vector<1>();
	std::vector<LocalShapeFunctionSet<2>*>& dynVec2 =
											get_dynamic_allocated_vector<2>();
	std::vector<LocalShapeFunctionSet<3>*>& dynVec3 =
											get_dynamic_allocated_vector<3>();

	for(size_t i = 0; i < dynVec1.size(); ++i) delete dynVec1[i];
	for(size_t i = 0; i < dynVec2.size(); ++i) delete dynVec2[i];
	for(size_t i = 0; i < dynVec3.size(); ++i) delete dynVec3[i];
};


} // namespace ug

