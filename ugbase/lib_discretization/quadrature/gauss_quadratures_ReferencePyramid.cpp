//This file is parsed from UG 3.9.


#include "quadrature.h"

namespace ug{

template <>
GaussQuadrature<ReferencePyramid>::GaussQuadrature(int order)
{
	switch(order)
	{
	case 2:
		m_order = 2;
		m_num_points = 8;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.585410200000000000;
		m_points[0][1] = 0.728196600000000000;
		m_points[0][2] = 0.138196600000000000;
		m_points[1][0] = 0.138196600000000000;
		m_points[1][1] = 0.728196600000000000;
		m_points[1][2] = 0.138196600000000000;
		m_points[2][0] = 0.138196600000000000;
		m_points[2][1] = 0.276309200000000000;
		m_points[2][2] = 0.585410200000000000;
		m_points[3][0] = 0.138196600000000000;
		m_points[3][1] = 0.276309200000000000;
		m_points[3][2] = 0.138196600000000000;
		m_points[4][0] = 0.728196600000000000;
		m_points[4][1] = 0.138196600000000000;
		m_points[4][2] = 0.138196600000000000;
		m_points[5][0] = 0.728196600000000000;
		m_points[5][1] = 0.585410200000000000;
		m_points[5][2] = 0.138196600000000000;
		m_points[6][0] = 0.276309200000000000;
		m_points[6][1] = 0.138196600000000000;
		m_points[6][2] = 0.585410200000000000;
		m_points[7][0] = 0.276309200000000000;
		m_points[7][1] = 0.138196600000000000;
		m_points[7][2] = 0.138196600000000000;

		m_weights[0] = 1./3. *  0.125000000000000000;
		m_weights[1] = 1./3. *  0.125000000000000000;
		m_weights[2] = 1./3. *  0.125000000000000000;
		m_weights[3] = 1./3. *  0.125000000000000000;
		m_weights[4] = 1./3. *  0.125000000000000000;
		m_weights[5] = 1./3. *  0.125000000000000000;
		m_weights[6] = 1./3. *  0.125000000000000000;
		m_weights[7] = 1./3. *  0.125000000000000000;
		break;

	default: assert(0 && "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}
}; // namespace ug

 // register quadratures at factory
namespace {
using namespace ug;

template <>
std::vector<const QuadratureRule<ReferencePyramid>* > QuadratureRuleFactory<ReferencePyramid>::m_rules =
	std::vector<const QuadratureRule<ReferencePyramid>* >();

GaussQuadrature<ReferencePyramid> gaussQuadratureReferencePyramid_2(2);

static const bool registered_2 = QuadratureRuleFactory<ReferencePyramid>::instance().register_rule(gaussQuadratureReferencePyramid_2);

};
