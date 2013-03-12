//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference pyramid.


#include "../quadrature.h"
#include "gauss_quad_pyramid.h"
#include "common/util/provider.h"

namespace ug{

GaussQuadrature<ReferencePyramid, 2>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.585410200000000000;
	m_vPoint[0][1] = 0.728196600000000000;
	m_vPoint[0][2] = 0.138196600000000000;

	m_vPoint[1][0] = 0.138196600000000000;
	m_vPoint[1][1] = 0.728196600000000000;
	m_vPoint[1][2] = 0.138196600000000000;

	m_vPoint[2][0] = 0.138196600000000000;
	m_vPoint[2][1] = 0.276309200000000000;
	m_vPoint[2][2] = 0.585410200000000000;

	m_vPoint[3][0] = 0.138196600000000000;
	m_vPoint[3][1] = 0.276309200000000000;
	m_vPoint[3][2] = 0.138196600000000000;

	m_vPoint[4][0] = 0.728196600000000000;
	m_vPoint[4][1] = 0.138196600000000000;
	m_vPoint[4][2] = 0.138196600000000000;

	m_vPoint[5][0] = 0.728196600000000000;
	m_vPoint[5][1] = 0.585410200000000000;
	m_vPoint[5][2] = 0.138196600000000000;

	m_vPoint[6][0] = 0.276309200000000000;
	m_vPoint[6][1] = 0.138196600000000000;
	m_vPoint[6][2] = 0.585410200000000000;

	m_vPoint[7][0] = 0.276309200000000000;
	m_vPoint[7][1] = 0.138196600000000000;
	m_vPoint[7][2] = 0.138196600000000000;

	m_vWeight[0] = 1./3. * 0.125000000000000000;
	m_vWeight[1] = 1./3. * 0.125000000000000000;
	m_vWeight[2] = 1./3. * 0.125000000000000000;
	m_vWeight[3] = 1./3. * 0.125000000000000000;
	m_vWeight[4] = 1./3. * 0.125000000000000000;
	m_vWeight[5] = 1./3. * 0.125000000000000000;
	m_vWeight[6] = 1./3. * 0.125000000000000000;
	m_vWeight[7] = 1./3. * 0.125000000000000000;
}




template <>
FlexGaussQuadrature<ReferencePyramid>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 2:{
		const static GaussQuadrature<ReferencePyramid, 2>& q2 
			= Provider<GaussQuadrature<ReferencePyramid, 2> >::get();

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
bool RegisterGaussQuadRule<ReferencePyramid>(QuadratureRuleProvider<ReferencePyramid::dim>& factory)
{
	static FlexGaussQuadrature<ReferencePyramid> gaussQuadratureReferencePyramid_2(2);

	bool success = true;
	factory.register_rule<ReferencePyramid>(gaussQuadratureReferencePyramid_2);

	return success;
};

}; // namespace ug

