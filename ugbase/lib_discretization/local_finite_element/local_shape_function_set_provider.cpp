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

template <typename TRefElem>
bool LocalShapeFunctionSetProvider::init_standard_sets()
{
//	create static Sets
	static LocalShapeFunctionSetWrapper<LagrangeP1<TRefElem, 1> > sSetLagrangeP1;
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

LocalShapeFunctionSetProvider::
LocalShapeFunctionSetProvider()
{
	static bool init = false;

	if(!init)
	{
	//	remember initialization
		init = true;

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
		static LocalShapeFunctionSetWrapper<LagrangeP1<ReferencePyramid, 1> > sSetLagrangeP1;
		if(!register_set(type1, sSetLagrangeP1))
			throw(UGFatalError("Cannot register Pyramid P1 Lagrange trial spaces."));

	}
};

} // namespace ug

