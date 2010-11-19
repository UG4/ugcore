//This file is parsed from UG 3.9.


#include "quadrature.h"

namespace ug{

template <>
GaussQuadrature<ReferencePrism>::GaussQuadrature(int order)
{
	switch(order)
	{
	case 0:
		m_order = 0;
		m_num_points = 6;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.000000000000000000;
		m_points[0][1] = 0.000000000000000000;
		m_points[0][2] = 0.000000000000000000;
		m_points[1][0] = 1.000000000000000000;
		m_points[1][1] = 0.000000000000000000;
		m_points[1][2] = 0.000000000000000000;
		m_points[2][0] = 0.000000000000000000;
		m_points[2][1] = 1.000000000000000000;
		m_points[2][2] = 0.000000000000000000;
		m_points[3][0] = 0.000000000000000000;
		m_points[3][1] = 0.000000000000000000;
		m_points[3][2] = 1.000000000000000000;
		m_points[4][0] = 1.000000000000000000;
		m_points[4][1] = 0.000000000000000000;
		m_points[4][2] = 1.000000000000000000;
		m_points[5][0] = 0.000000000000000000;
		m_points[5][1] = 1.000000000000000000;
		m_points[5][2] = 1.000000000000000000;

		m_weights[0] = 0.5 *  0.166666666666666660;
		m_weights[1] = 0.5 *  0.166666666666666660;
		m_weights[2] = 0.5 *  0.166666666666666660;
		m_weights[3] = 0.5 *  0.166666666666666660;
		m_weights[4] = 0.5 *  0.166666666666666660;
		m_weights[5] = 0.5 *  0.166666666666666660;
		break;

	case 2:
		m_order = 2;
		m_num_points = 6;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.666666666666666660;
		m_points[0][1] = 0.166666666666666660;
		m_points[0][2] = 0.211324865405187000;
		m_points[1][0] = 0.166666666666666660;
		m_points[1][1] = 0.666666666666666660;
		m_points[1][2] = 0.211324865405187000;
		m_points[2][0] = 0.166666666666666660;
		m_points[2][1] = 0.166666666666666660;
		m_points[2][2] = 0.211324865405187000;
		m_points[3][0] = 0.666666666666666660;
		m_points[3][1] = 0.166666666666666660;
		m_points[3][2] = 0.788675134594813000;
		m_points[4][0] = 0.166666666666666660;
		m_points[4][1] = 0.666666666666666660;
		m_points[4][2] = 0.788675134594813000;
		m_points[5][0] = 0.166666666666666660;
		m_points[5][1] = 0.166666666666666660;
		m_points[5][2] = 0.788675134594813000;

		m_weights[0] = 0.5 *  0.166666666666666660;
		m_weights[1] = 0.5 *  0.166666666666666660;
		m_weights[2] = 0.5 *  0.166666666666666660;
		m_weights[3] = 0.5 *  0.166666666666666660;
		m_weights[4] = 0.5 *  0.166666666666666660;
		m_weights[5] = 0.5 *  0.166666666666666660;
		break;

	default: assert(0 && "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}

template <>
bool RegisterQuadratureRule(QuadratureRuleProvider<ReferencePrism>& factory)
{
	static GaussQuadrature<ReferencePrism> gaussQuadratureReferencePrism_0(0);
	static GaussQuadrature<ReferencePrism> gaussQuadratureReferencePrism_2(2);

	bool success = true;
	success &= factory.register_rule(gaussQuadratureReferencePrism_0);
	success &= factory.register_rule(gaussQuadratureReferencePrism_2);
	return success;
}

}; // namespace ug

/*
 // register quadratures at factory
namespace {
using namespace ug;

template <>
std::vector<const QuadratureRule<ReferencePrism>* > QuadratureRuleProvider<ReferencePrism>::m_rules =
	std::vector<const QuadratureRule<ReferencePrism>* >();
};
*/
