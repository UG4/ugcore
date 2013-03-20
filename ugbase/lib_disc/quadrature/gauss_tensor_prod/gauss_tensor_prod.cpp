/*
 * gauss_tensor_prod.cpp
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau
 */


#include "../quadrature.h"
#include "common/util/provider.h"
#include "gauss_tensor_prod.h"
#include "../gauss_legendre/gauss_legendre.h"
#include "../gauss_jacobi/gauss_jacobi10.h"
#include "../gauss_jacobi/gauss_jacobi20.h"

namespace ug {

GaussQuadratureTriangle::GaussQuadratureTriangle(int order)
{
	GaussLegendre quadRule = GaussLegendre(order);
	GaussJacobi10 quadRule10 = GaussJacobi10(order);

	m_order = std::min(quadRule.order(), quadRule10.order());
	m_numPoints = quadRule10.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	int cnt = 0;
	for(int i = 0; i < quadRule.size(); i++){
		for(int j = 0; j < quadRule10.size(); j++, cnt++){
			pvPoint[cnt][0] = quadRule10.point(j)[0];
			pvPoint[cnt][1] = (1.0 - quadRule10.point(j)[0] ) * quadRule.point(i)[0];
			pvWeight[cnt] = quadRule.weight(i) * quadRule10.weight(j);
		}
	}

};

GaussQuadratureTriangle::~GaussQuadratureTriangle()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadratureQuadrilateral::GaussQuadratureQuadrilateral(int order)
{
	GaussLegendre quadRule = GaussLegendre(order);

	m_order = std::min(quadRule.order(), quadRule.order());
	m_numPoints = quadRule.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	int cnt  = 0;
	for(int i = 0; i < quadRule.size(); i ++) {
		for(int j = 0; j < quadRule.size(); j++, cnt++) {
			pvPoint[cnt][0] = quadRule.point(i)[0];
			pvPoint[cnt][1] = quadRule.point(j)[0];
			pvWeight[cnt] = quadRule.weight(i) * quadRule.weight(j);
		}
	}

};

GaussQuadratureQuadrilateral::~GaussQuadratureQuadrilateral()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadratureHexahedron::GaussQuadratureHexahedron(int order)
{
	GaussLegendre quadRule = GaussLegendre(order);

	m_order = std::min(quadRule.order(), std::min(quadRule.order(), quadRule.order()));
	m_numPoints = quadRule.size() * quadRule.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	int cnt  = 0;
	for(int i = 0; i < quadRule.size(); i ++) {
		for(int j = 0; j < quadRule.size(); j++) {
			for(int k = 0; k < quadRule.size(); k++, cnt++) {
				pvPoint[cnt][0] = quadRule.point(i)[0];
				pvPoint[cnt][1] = quadRule.point(j)[0];
				pvPoint[cnt][2] = quadRule.point(k)[0];
				pvWeight[cnt] = quadRule.weight(i) * quadRule.weight(j) * quadRule.weight(k);
			}
		}
	}

};

GaussQuadratureHexahedron::~GaussQuadratureHexahedron()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadratureTetrahedron::GaussQuadratureTetrahedron(int order)
{
	GaussLegendre quadRule = GaussLegendre(order);
	GaussJacobi10 quadRule10 = GaussJacobi10(order);
	GaussJacobi20 quadRule20 = GaussJacobi20(order);

	m_order = std::min(quadRule.order(), std::min(quadRule10.order(), quadRule20.order()));
	m_numPoints = quadRule.size() * quadRule10.size() * quadRule20.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	int cnt = 0;
	for(int i = 0; i < quadRule20.size(); i++) {
		for(int j = 0; j < quadRule10.size(); j++) {
			for(int k = 0; k < quadRule.size(); k++, cnt++) {
				pvPoint[cnt][0] = quadRule20.point(i)[0];
				pvPoint[cnt][1] = (1.0 - quadRule20.point(i)[0] ) * quadRule10.point(j)[0];
				pvPoint[cnt][2] = (1.0 - quadRule20.point(i)[0]) * (1.0 - quadRule10.point(j)[0]) * quadRule.point(k)[0];
				pvWeight[cnt] = quadRule20.weight(i) * quadRule10.weight(j) * quadRule.weight(k);
			}
		}
	}

};

GaussQuadratureTetrahedron::~GaussQuadratureTetrahedron()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

GaussQuadraturePrism::GaussQuadraturePrism(int order)
{
	GaussLegendre quadRule = GaussLegendre(order);
	GaussJacobi10 quadRule10 = GaussJacobi10(order);

	m_order = std::min(quadRule.order(), quadRule10.order());
	m_numPoints = quadRule10.size() * quadRule.size() * quadRule.size();
	MathVector<dim>* pvPoint; number* pvWeight;
	m_pvPoint = pvPoint = new position_type[m_numPoints];
	m_pvWeight = pvWeight = new weight_type[m_numPoints];

	int cnt = 0;
	for(int i = 0; i < quadRule10.size(); i++) {
		for(int j = 0; j < quadRule.size(); j++) {
			for(int k = 0; k < quadRule.size(); k++, cnt++) {
				pvPoint[cnt][0] = quadRule10.point(i)[0];
				pvPoint[cnt][1] = (1.0 - quadRule10.point(i)[0]) * quadRule.point(j)[0];
				pvPoint[cnt][2] = quadRule.point(k)[0];
				pvWeight[cnt] = quadRule10.weight(i) * quadRule.weight(j) * quadRule.weight(k);
			}
		}
	}

};

GaussQuadraturePrism::~GaussQuadraturePrism()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
}

}
