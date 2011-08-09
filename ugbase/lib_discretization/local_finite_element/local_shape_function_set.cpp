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

namespace ug{

template <typename TRefElem>
bool LocalShapeFunctionSetProvider::init_standard_sets()
{
//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem> > sSetLagrangeP1;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 2> > sSetLagrangeP2;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 3> > sSetLagrangeP3;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 4> > sSetLagrangeP4;

//	insert into map: P1 Lagrange
	LFEID type1(LFEID::LAGRANGE, 1);
	if(!register_set(type1, sSetLagrangeP1))
		return false;

//	insert into map: P2 Lagrange
	LFEID type2(LFEID::LAGRANGE, 2);
	if(!register_set(type2, sSetLagrangeP2))
		return false;

//	insert into map: P3 Lagrange
	LFEID type3(LFEID::LAGRANGE, 3);
	if(!register_set(type3, sSetLagrangeP3))
		return false;

//	insert into map: P4 Lagrange
	LFEID type4(LFEID::LAGRANGE, 4);
	if(!register_set(type4, sSetLagrangeP4))
		return false;

//	return success
	return true;
}

template <typename TRefElem>
bool LocalShapeFunctionSetProvider::init_flex_lagrange(size_t order)
{
//	create static Sets
	LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >* sSetFlexLagrange
		= new LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >;
	sSetFlexLagrange->set_order(order);

//	insert into map: Lagrange
	LFEID type(LFEID::LAGRANGE, order);
	if(!register_set(type, *sSetFlexLagrange))
		return false;

//	remember creation
	get_dynamic_allocated_vector<TRefElem::dim>().push_back(sSetFlexLagrange);

//	return success
	return true;
}

void LocalShapeFunctionSetProvider::
dynamically_create_set(ReferenceObjectID roid, LFEID id)
{
//	Lagrange space
	if(id.type() == LFEID::LAGRANGE)
	{
	//	only order >= 1 available
		if(id.order() < 1) return;

	//	switch type
		switch(roid)
		{
			case ROID_EDGE:
				if(!init_flex_lagrange<ReferenceEdge>(id.order()))
					throw(UGFatalError("Dynamic Allocation of set failed."));
				return;
			case ROID_TRIANGLE:
				if(!init_flex_lagrange<ReferenceTriangle>(id.order()))
					throw(UGFatalError("Dynamic Allocation of set failed."));
				return;
			case ROID_QUADRILATERAL:
				if(!init_flex_lagrange<ReferenceQuadrilateral>(id.order()))
					throw(UGFatalError("Dynamic Allocation of set failed."));
				return;
			case ROID_TETRAHEDRON:
				if(!init_flex_lagrange<ReferenceTetrahedron>(id.order()))
					throw(UGFatalError("Dynamic Allocation of set failed."));
				return;
			case ROID_PRISM:
				if(!init_flex_lagrange<ReferencePrism>(id.order()))
					throw(UGFatalError("Dynamic Allocation of set failed."));
				return;
			case ROID_HEXAHEDRON:
				if(!init_flex_lagrange<ReferenceHexahedron>(id.order()))
					throw(UGFatalError("Dynamic Allocation of set failed."));
				return;
			default: return;
		}
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

	//	clear all maps
		clear_maps<ReferenceEdge>();
		clear_maps<ReferenceTriangle>();
		clear_maps<ReferenceQuadrilateral>();
		clear_maps<ReferenceTetrahedron>();
		clear_maps<ReferencePrism>();
		clear_maps<ReferencePyramid>();
		clear_maps<ReferenceHexahedron>();

	//	clear created holder
		get_dynamic_allocated_vector<1>().clear();
		get_dynamic_allocated_vector<2>().clear();
		get_dynamic_allocated_vector<3>().clear();

	//	register all element types that allow higher orders
		if(!init_standard_sets<ReferenceEdge>())
			throw(UGFatalError("Cannot register standard Edge trial spaces."));
		if(!init_standard_sets<ReferenceTriangle>())
			throw(UGFatalError("Cannot register standard Triangle trial spaces."));
		if(!init_standard_sets<ReferenceQuadrilateral>())
			throw(UGFatalError("Cannot register standard Quadrilateral trial spaces."));
		if(!init_standard_sets<ReferenceTetrahedron>())
			throw(UGFatalError("Cannot register standard Tetrahedron trial spaces."));
		if(!init_standard_sets<ReferencePrism>())
			throw(UGFatalError("Cannot register standard Prism trial spaces."));
		if(!init_standard_sets<ReferenceHexahedron>())
			throw(UGFatalError("Cannot register standard Hexahedron trial spaces."));

	//	register 1st order pyramid
		LFEID type1(LFEID::LAGRANGE, 1);
		static LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid> > sSetLagrangeP1;
		if(!register_set(type1, sSetLagrangeP1))
			throw(UGFatalError("Cannot register Pyramid P1 Lagrange trial spaces."));
	}
};

LocalShapeFunctionSetProvider::
~LocalShapeFunctionSetProvider()
{
	std::vector<DimLocalShapeFunctionSet<1>*>& dynVec1 =
											get_dynamic_allocated_vector<1>();
	std::vector<DimLocalShapeFunctionSet<2>*>& dynVec2 =
											get_dynamic_allocated_vector<2>();
	std::vector<DimLocalShapeFunctionSet<3>*>& dynVec3 =
											get_dynamic_allocated_vector<3>();

	for(size_t i = 0; i < dynVec1.size(); ++i) delete dynVec1[i];
	for(size_t i = 0; i < dynVec2.size(); ++i) delete dynVec2[i];
	for(size_t i = 0; i < dynVec3.size(); ++i) delete dynVec3[i];
};


} // namespace ug

