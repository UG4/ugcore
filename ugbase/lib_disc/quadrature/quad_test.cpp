/*
 * quad_test.cpp
 *
 *  Created on: 13.03.2013
 *      Author: andreasvogel
 */

#include "quad_test.h"
#include "gauss/gauss_quad.h"
#include "quadrature.h"
#include "quadrature_provider.h"

#include <iostream>
using namespace std;

namespace ug{

void PrintQuadRule(const QuadratureRule<2>& q)
{
	cout << "QuadRule\n";
	cout << "--------\n";
	cout << "Order: " << q.order() << "\n";
	cout << "Size:  " << q.size() << "\n";
	for(size_t i = 0; i < q.size(); ++i)
		cout << "Weight "<<i<<": "<<q.weight(i) <<"\n";
	for(size_t i = 0; i < q.size(); ++i)
		cout << "Point  "<<i<<": "<<q.point(i) <<"\n";

}

void TestQuadRule()
{
	typedef GaussQuadrature<ReferenceTriangle, 1> Q;
	cout << "Order: " << Q::order() << "\n";
	cout << "Size:  " << Q::size() << "\n";
	for(size_t i = 0; i < Q::size(); ++i)
		cout << "Weight "<<i<<": "<<Q::weight(i) <<"\n";
	for(size_t i = 0; i < Q::size(); ++i)
		cout << "Point  "<<i<<": "<<Q::point(i) <<"\n";

	GaussQuadrature<ReferenceTriangle, 2> q2;
	cout << "Order: " << q2.order() << "\n";
	cout << "Size:  " << q2.size() << "\n";
	for(size_t i = 0; i < q2.size(); ++i)
		cout << "Weight "<<i<<": "<<q2.weight(i) <<"\n";
	for(size_t i = 0; i < q2.size(); ++i)
		cout << "Point  "<<i<<": "<<q2.point(i) <<"\n";

	PrintQuadRule(FlexGaussQuadrature<ReferenceTriangle>(1));
	PrintQuadRule(FlexGaussQuadrature<ReferenceTriangle>(2));
	PrintQuadRule(FlexGaussQuadrature<ReferenceTriangle>(3));

	PrintQuadRule(QuadratureRuleProvider<2>::get_rule(ROID_TRIANGLE, 3));
	PrintQuadRule(QuadratureRuleProvider<2>::get_rule<ReferenceTriangle>(3));

}

} // end namespace ug
