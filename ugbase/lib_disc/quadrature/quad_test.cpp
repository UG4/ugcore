/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
	using Q = GaussQuadrature<ReferenceTriangle, 1>;
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

	PrintQuadRule(QuadratureRuleProvider<2>::get(ROID_TRIANGLE, 3));
	PrintQuadRule(QuadratureRuleProvider<2>::get<ReferenceTriangle>(3));

}

} // end namespace ug
