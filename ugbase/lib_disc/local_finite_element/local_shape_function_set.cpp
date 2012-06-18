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

namespace ug{

template <typename TRefElem>
void LocalShapeFunctionSetProvider::init_standard_sets()
{
//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem> > sSetLagrangeP1;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 2> > sSetLagrangeP2;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 3> > sSetLagrangeP3;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 4> > sSetLagrangeP4;

//	insert into map: P1 Lagrange
	LFEID type1(LFEID::LAGRANGE, 1);
	register_set(type1, sSetLagrangeP1);

//	insert into map: P2 Lagrange
	LFEID type2(LFEID::LAGRANGE, 2);
	register_set(type2, sSetLagrangeP2);

//	insert into map: P3 Lagrange
	LFEID type3(LFEID::LAGRANGE, 3);
	register_set(type3, sSetLagrangeP3);

//	insert into map: P4 Lagrange
	LFEID type4(LFEID::LAGRANGE, 4);
	register_set(type4, sSetLagrangeP4);


//	insert into map: Crouzeix-Raviart
	static LocalShapeFunctionSetWrapper<CrouzeixRaviartLSFS<TRefElem> > sSetCrouzeixRaviart;
	LFEID typeCR(LFEID::CROUZEIX_RAVIART, 1);
	register_set(typeCR, sSetCrouzeixRaviart);

}

template <typename TRefElem>
void LocalShapeFunctionSetProvider::init_flex_lagrange(size_t order)
{
//	create static Sets
	LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >* sSetFlexLagrange
		= new LocalShapeFunctionSetWrapper<FlexLagrangeLSFS<TRefElem> >;
	sSetFlexLagrange->set_order(order);

//	insert into map: Lagrange
	LFEID type(LFEID::LAGRANGE, order);
	register_set(type, *sSetFlexLagrange);

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
			case ROID_EDGE:
					init_flex_lagrange<ReferenceEdge>(id.order());
				return;
			case ROID_TRIANGLE:
					init_flex_lagrange<ReferenceTriangle>(id.order());
				return;
			case ROID_QUADRILATERAL:
				init_flex_lagrange<ReferenceQuadrilateral>(id.order());
				return;
			case ROID_TETRAHEDRON:
				init_flex_lagrange<ReferenceTetrahedron>(id.order());
				return;
			case ROID_PRISM:
				init_flex_lagrange<ReferencePrism>(id.order());
				return;
			case ROID_HEXAHEDRON:
				init_flex_lagrange<ReferenceHexahedron>(id.order());
				return;
			default: return;
		}

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
			register_set(type1, sSetLagrangeP1);
		}UG_CATCH_THROW("Cannot register Pyramid P1 Lagrange trial spaces.");
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

