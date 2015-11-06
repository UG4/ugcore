/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ROTATING_CONE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ROTATING_CONE__

#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"

namespace ug{

class RotatingCone2d
	: public StdGlobPosData<RotatingCone2d, number, 2>
{
	public:
		RotatingCone2d(double _eps, double _cx, double _cy, double _ax, double _ay,
		             double _nu, double _delta)
			: eps(_eps), cx(_cx), cy(_cy), ax(_ax), ay(_ay), nu(_nu), delta(_delta)
		{};

		inline void evaluate(number& val, const MathVector<2>& pos, number time, int si) const
		{
			const number t = time;
			const number x = pos[0];
			const number y = pos[1];

			const number xRot = cos(nu*t) * (x-cx) - sin(nu*t) * (y-cy);
			const number yRot = sin(nu*t) * (x-cx) + cos(nu*t) * (y-cy);

			const number expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t);
			const number scale = delta/(delta+4*eps*t);

			val = scale * exp(expo);
		}

	protected:
		double eps, cx, cy, ax, ay, nu, delta;
};

}
#endif // __H__UG__LIB_DISC__SPATIAL_DISC__ROTATING_CONE__
