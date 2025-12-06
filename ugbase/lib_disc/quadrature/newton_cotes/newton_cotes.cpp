/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
 * Author: Lisa Grau
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

#include "newton_cotes.h"
#include "../quadrature.h"
//ø #include "common/util/provider.h"

namespace ug
{

//constructor
NewtonCotes::NewtonCotes(size_t order)
{

	if (order == 0) order = 1;
	m_order = order;
	m_numPoints = order+1;
	position_type* pvPoint= new position_type[m_numPoints];
	weight_type* pvWeight= new weight_type[m_numPoints];
	m_pvPoint = pvPoint;
	m_pvWeight = pvWeight;

	switch(m_numPoints)
	{
		case 0:
		case 1:
		case 2:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 1.;

			pvWeight[0] = 0.5;
			pvWeight[1] = 0.5;
		break;

		case 3:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.5;
			pvPoint[2][0] = 1.;

			pvWeight[0] = 0.1666666666666666666666666666666666666667;
			pvWeight[1] = 0.6666666666666666666666666666666666666667;
			pvWeight[2] = 0.1666666666666666666666666666666666666667;
		break;

		case 4:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.3333333333333333333333333333333333333333;
			pvPoint[2][0] = 0.6666666666666666666666666666666666666667;
			pvPoint[3][0] = 1.;

			pvWeight[0] = 0.125;
			pvWeight[1] = 0.375;
			pvWeight[2] = 0.375;
			pvWeight[3] = 0.125;
		break;

		case 5:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.25;
			pvPoint[2][0] = 0.5;
			pvPoint[3][0] = 0.75;
			pvPoint[4][0] = 1.;

			pvWeight[0] = 0.07777777777777777777777777777777777777778;
			pvWeight[1] = 0.3555555555555555555555555555555555555556;
			pvWeight[2] = 0.1333333333333333333333333333333333333333;
			pvWeight[3] = 0.3555555555555555555555555555555555555556;
			pvWeight[4] = 0.07777777777777777777777777777777777777778;
		break;

		case 6:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.2;
			pvPoint[2][0] = 0.4;
			pvPoint[3][0] = 0.6;
			pvPoint[4][0] = 0.8;
			pvPoint[5][0] = 1.;

			pvWeight[0] = 0.06597222222222222222222222222222222222222;
			pvWeight[1] = 0.2604166666666666666666666666666666666667;
			pvWeight[2] = 0.1736111111111111111111111111111111111111;
			pvWeight[3] = 0.1736111111111111111111111111111111111111;
			pvWeight[4] = 0.2604166666666666666666666666666666666667;
			pvWeight[5] = 0.06597222222222222222222222222222222222222;
		break;

		case 7:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.1666666666666666666666666666666666666667;
			pvPoint[2][0] = 0.3333333333333333333333333333333333333333;
			pvPoint[3][0] = 0.5;
			pvPoint[4][0] = 0.6666666666666666666666666666666666666667;
			pvPoint[5][0] = 0.8333333333333333333333333333333333333333;
			pvPoint[6][0] = 1.;

			pvWeight[0] = 0.04880952380952380952380952380952380952381;
			pvWeight[1] = 0.2571428571428571428571428571428571428571;
			pvWeight[2] = 0.03214285714285714285714285714285714285714;
			pvWeight[3] = 0.3238095238095238095238095238095238095238;
			pvWeight[4] = 0.03214285714285714285714285714285714285714;
			pvWeight[5] = 0.2571428571428571428571428571428571428571;
			pvWeight[6] = 0.04880952380952380952380952380952380952381;
		break;

		case 8:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.1428571428571428571428571428571428571429;
			pvPoint[2][0] = 0.2857142857142857142857142857142857142857;
			pvPoint[3][0] = 0.4285714285714285714285714285714285714286;
			pvPoint[4][0] = 0.5714285714285714285714285714285714285714;
			pvPoint[5][0] = 0.7142857142857142857142857142857142857143;
			pvPoint[6][0] = 0.8571428571428571428571428571428571428571;
			pvPoint[7][0] = 1.;

			pvWeight[0] = 0.04346064814814814814814814814814814814815;
			pvWeight[1] = 0.2070023148148148148148148148148148148148;
			pvWeight[2] = 0.0765625;
			pvWeight[3] = 0.172974537037037037037037037037037037037;
			pvWeight[4] = 0.172974537037037037037037037037037037037;
			pvWeight[5] = 0.0765625;
			pvWeight[6] = 0.2070023148148148148148148148148148148148;
			pvWeight[7] = 0.04346064814814814814814814814814814814815;
		break;

		case 9:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.125;
			pvPoint[2][0] = 0.25;
			pvPoint[3][0] = 0.375;
			pvPoint[4][0] = 0.5;
			pvPoint[5][0] = 0.625;
			pvPoint[6][0] = 0.75;
			pvPoint[7][0] = 0.875;
			pvPoint[8][0] = 1.;

			pvWeight[0] = 0.03488536155202821869488536155202821869489;
			pvWeight[1] = 0.2076895943562610229276895943562610229277;
			pvWeight[2] = -0.03273368606701940035273368606701940035273;
			pvWeight[3] = 0.3702292768959435626102292768959435626102;
			pvWeight[4] = -0.1601410934744268077601410934744268077601;
			pvWeight[5] = 0.3702292768959435626102292768959435626102;
			pvWeight[6] = -0.03273368606701940035273368606701940035273;
			pvWeight[7] = 0.2076895943562610229276895943562610229277;
			pvWeight[8] = 0.03488536155202821869488536155202821869489;
		break;

		case 10:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.1111111111111111111111111111111111111111;
			pvPoint[2][0] = 0.2222222222222222222222222222222222222222;
			pvPoint[3][0] = 0.3333333333333333333333333333333333333333;
			pvPoint[4][0] = 0.4444444444444444444444444444444444444444;
			pvPoint[5][0] = 0.5555555555555555555555555555555555555556;
			pvPoint[6][0] = 0.6666666666666666666666666666666666666667;
			pvPoint[7][0] = 0.7777777777777777777777777777777777777778;
			pvPoint[8][0] = 0.8888888888888888888888888888888888888889;
			pvPoint[9][0] = 1.;

			pvWeight[0] = 0.03188616071428571428571428571428571428571;
			pvWeight[1] = 0.1756808035714285714285714285714285714286;
			pvWeight[2] = 0.01205357142857142857142857142857142857143;
			pvWeight[3] = 0.2158928571428571428571428571428571428571;
			pvWeight[4] = 0.06448660714285714285714285714285714285714;
			pvWeight[5] = 0.06448660714285714285714285714285714285714;
			pvWeight[6] = 0.2158928571428571428571428571428571428571;
			pvWeight[7] = 0.01205357142857142857142857142857142857143;
			pvWeight[8] = 0.1756808035714285714285714285714285714286;
			pvWeight[9] = 0.03188616071428571428571428571428571428571;
		break;

		case 11:
			pvPoint[0][0] = 0;
			pvPoint[1][0] = 0.1;
			pvPoint[2][0] = 0.2;
			pvPoint[3][0] = 0.3;
			pvPoint[4][0] = 0.4;
			pvPoint[5][0] = 0.5;
			pvPoint[6][0] = 0.6;
			pvPoint[7][0] = 0.7;
			pvPoint[8][0] = 0.8;
			pvPoint[9][0] = 0.9;
			pvPoint[10][0] = 1.;

			pvWeight[0] = 0.02683414836192613970391748169525947303725;
			pvWeight[1] = 0.1775359414248303137192026080914969803859;
			pvWeight[2] = -0.08104357062690396023729357062690396023729;
			pvWeight[3] = 0.4549462882796216129549462882796216129549;
			pvWeight[4] = -0.4351551226551226551226551226551226551227;
			pvWeight[5] = 0.7137646304312970979637646304312970979638;
			pvWeight[6] = -0.4351551226551226551226551226551226551227;
			pvWeight[7] = 0.4549462882796216129549462882796216129549;
			pvWeight[8] = -0.08104357062690396023729357062690396023729;
			pvWeight[9] = 0.1775359414248303137192026080914969803859;
			pvWeight[10] = 0.02683414836192613970391748169525947303725;
		break;

		default:
			UG_THROW("NewtonCotes: Rule for order " << order  << " not supported.");
	}
}

NewtonCotes::~NewtonCotes()
{
	delete[] m_pvPoint;
	delete[] m_pvWeight;
};
}
