/*
 * local_shape_function_set.cpp
 *
 *  Created on: 17.02.2010
 *      Author: andreasvogel
 */

#include "local_shape_function_set_provider.h"

// include spaces
#include "lagrange/lagrangep1.h"
#include "lagrange/lagrange.h"

namespace ug{

const LocalShapeFunctionSetBase&
LocalShapeFunctionSetProvider::get(LSFSID id, ReferenceObjectID roid)
{
//	init provider and search for identifier
	BaseMap::const_iterator iter = inst().m_baseMap.find(id);

//	if not found
	if(iter == m_baseMap.end())
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Base Trial Space Type "<<id<<" requested for"
				" Reference Element type " <<roid<<".\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	get vector
	const std::vector<const LocalShapeFunctionSetBase*>& vBase = iter->second;

//	check that space registered
	if(vBase[roid] == NULL)
	{
		UG_LOG("ERROR in 'LocalShapeFunctionSetProvider::get': "
				"Unknown Base Trial Space for Type "<<id<<" and Reference"
				" element type "<<roid<<" requested.\n");
		throw(UGFatalError("Trial Space type unknown"));
	}

//	return shape function set
	return *(vBase[roid]);
}

template <typename TRefElem>
bool
LocalShapeFunctionSetProvider::
init_standard_local_shape_function_sets()
{
//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem, 1> > sSetLagrangeP1;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 2> > sSetLagrangeP2;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 3> > sSetLagrangeP3;
	static LocalShapeFunctionSetWrapper<LagrangeLSFS<TRefElem, 4> > sSetLagrangeP4;

//	insert into map: P1 Lagrange
	LSFSID type1(LSFSID::LAGRANGE, 1);
	if(!register_local_shape_function_set(type1, sSetLagrangeP1))
		return false;

//	insert into map: P2 Lagrange
	LSFSID type2(LSFSID::LAGRANGE, 2);
	if(!register_local_shape_function_set(type2, sSetLagrangeP2))
		return false;

//	insert into map: P3 Lagrange
	LSFSID type3(LSFSID::LAGRANGE, 3);
	if(!register_local_shape_function_set(type3, sSetLagrangeP3))
		return false;

//	insert into map: P4 Lagrange
	LSFSID type4(LSFSID::LAGRANGE, 4);
	if(!register_local_shape_function_set(type4, sSetLagrangeP4))
		return false;

//	return success
	return true;
}

LocalShapeFunctionSetProvider::
LocalShapeFunctionSetProvider()
{
	static bool init = false;

	if(!init)
	{
	//	remember initialization
		init = true;

	//	register all element types that allow higher orders
		if(!init_standard_local_shape_function_sets<ReferenceEdge>())
			throw(UGFatalError("Cannot register standard Edge trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceTriangle>())
			throw(UGFatalError("Cannot register standard Triangle trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceQuadrilateral>())
			throw(UGFatalError("Cannot register standard Quadrilateral trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceTetrahedron>())
			throw(UGFatalError("Cannot register standard Tetrahedron trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferencePrism>())
			throw(UGFatalError("Cannot register standard Prism trial spaces."));
		if(!init_standard_local_shape_function_sets<ReferenceHexahedron>())
			throw(UGFatalError("Cannot register standard Hexahedron trial spaces."));

	//	register 1st order pyramid
		LSFSID type1(LSFSID::LAGRANGE, 1);
		static LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid, 1> > sSetLagrangeP1;
		if(!register_local_shape_function_set(type1, sSetLagrangeP1))
			throw(UGFatalError("Cannot register Pyramid P1 Lagrange trial spaces."));

	}
};

std::map<LSFSID, std::vector<const LocalShapeFunctionSetBase*> >
LocalShapeFunctionSetProvider::m_baseMap =
		std::map<LSFSID, std::vector<const LocalShapeFunctionSetBase*> >();


} // namespace ug

