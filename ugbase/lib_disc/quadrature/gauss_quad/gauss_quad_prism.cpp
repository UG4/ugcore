//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference prism.


#include "../quadrature.h"
#include "gauss_quad_prism.h"
#include "common/util/provider.h"

namespace ug{

GaussQuadrature<ReferencePrism, 0>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.000000000000000000;
	m_vPoint[0][1] = 0.000000000000000000;
	m_vPoint[0][2] = 0.000000000000000000;

	m_vPoint[1][0] = 1.000000000000000000;
	m_vPoint[1][1] = 0.000000000000000000;
	m_vPoint[1][2] = 0.000000000000000000;

	m_vPoint[2][0] = 0.000000000000000000;
	m_vPoint[2][1] = 1.000000000000000000;
	m_vPoint[2][2] = 0.000000000000000000;

	m_vPoint[3][0] = 0.000000000000000000;
	m_vPoint[3][1] = 0.000000000000000000;
	m_vPoint[3][2] = 1.000000000000000000;

	m_vPoint[4][0] = 1.000000000000000000;
	m_vPoint[4][1] = 0.000000000000000000;
	m_vPoint[4][2] = 1.000000000000000000;

	m_vPoint[5][0] = 0.000000000000000000;
	m_vPoint[5][1] = 1.000000000000000000;
	m_vPoint[5][2] = 1.000000000000000000;

	m_vWeight[0] = 0.5 * 0.166666666666666660;
	m_vWeight[1] = 0.5 * 0.166666666666666660;
	m_vWeight[2] = 0.5 * 0.166666666666666660;
	m_vWeight[3] = 0.5 * 0.166666666666666660;
	m_vWeight[4] = 0.5 * 0.166666666666666660;
	m_vWeight[5] = 0.5 * 0.166666666666666660;
}

GaussQuadrature<ReferencePrism, 2>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.666666666666666660;
	m_vPoint[0][1] = 0.166666666666666660;
	m_vPoint[0][2] = 0.211324865405187000;

	m_vPoint[1][0] = 0.166666666666666660;
	m_vPoint[1][1] = 0.666666666666666660;
	m_vPoint[1][2] = 0.211324865405187000;

	m_vPoint[2][0] = 0.166666666666666660;
	m_vPoint[2][1] = 0.166666666666666660;
	m_vPoint[2][2] = 0.211324865405187000;

	m_vPoint[3][0] = 0.666666666666666660;
	m_vPoint[3][1] = 0.166666666666666660;
	m_vPoint[3][2] = 0.788675134594813000;

	m_vPoint[4][0] = 0.166666666666666660;
	m_vPoint[4][1] = 0.666666666666666660;
	m_vPoint[4][2] = 0.788675134594813000;

	m_vPoint[5][0] = 0.166666666666666660;
	m_vPoint[5][1] = 0.166666666666666660;
	m_vPoint[5][2] = 0.788675134594813000;

	m_vWeight[0] = 0.5 * 0.166666666666666660;
	m_vWeight[1] = 0.5 * 0.166666666666666660;
	m_vWeight[2] = 0.5 * 0.166666666666666660;
	m_vWeight[3] = 0.5 * 0.166666666666666660;
	m_vWeight[4] = 0.5 * 0.166666666666666660;
	m_vWeight[5] = 0.5 * 0.166666666666666660;
}




template <>
FlexGaussQuadrature<ReferencePrism>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 0:{
		const static GaussQuadrature<ReferencePrism, 0>& q0 
			= Provider<GaussQuadrature<ReferencePrism, 0> >::get();

		m_order = q0.order();
		m_numPoints = q0.size();
		m_pvPoint = q0.points();
		m_pvWeight = q0.weights();
		}break;

	case 2:{
		const static GaussQuadrature<ReferencePrism, 2>& q2 
			= Provider<GaussQuadrature<ReferencePrism, 2> >::get();

		m_order = q2.order();
		m_numPoints = q2.size();
		m_pvPoint = q2.points();
		m_pvWeight = q2.weights();
		}break;

	default: UG_ASSERT(0, "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}



// register rules
template <>
bool RegisterGaussQuadRule<ReferencePrism>(QuadratureRuleProvider<ReferencePrism::dim>& factory)
{
	static FlexGaussQuadrature<ReferencePrism> gaussQuadratureReferencePrism_0(0);
	static FlexGaussQuadrature<ReferencePrism> gaussQuadratureReferencePrism_2(2);

	bool success = true;
	factory.register_rule<ReferencePrism>(gaussQuadratureReferencePrism_0);
	factory.register_rule<ReferencePrism>(gaussQuadratureReferencePrism_2);

	return success;
};

}; // namespace ug

