/*
 * rotating_cone.h
 *
 *      Author: andreasvogel
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
