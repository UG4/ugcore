//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference edge.


#include "gauss_quad.h"

namespace ug{

class GaussQuadratureVertex : public QuadratureRule<0>
{
	public:
	/// Constructor
		GaussQuadratureVertex(){
			m_order = 20; // something large
			m_numPoints = 1;
			m_vPoint[0][0] = 0.000000000000000000;
			m_vWeight[0] =  1.000000000000000000;
			m_pvPoint = &m_vPoint[0];
			m_pvWeight = &m_vWeight[0];
		}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[1];

	/// weights
		number m_vWeight[1];
};

}; // namespace ug

