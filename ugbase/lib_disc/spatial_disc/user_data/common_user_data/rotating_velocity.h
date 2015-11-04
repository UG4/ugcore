
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ROTATING_VELOCITY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ROTATING_VELOCITY__

#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"

namespace ug{

class RotatingVelocity2d
	: public StdGlobPosData<RotatingVelocity2d, MathVector<2>, 2>
{
	public:
		RotatingVelocity2d(double _cx, double _cy, double _nu)
			: cx(_cx), cy(_cy), nu(_nu)
		{};

		inline void evaluate(MathVector<2>& val, const MathVector<2>& pos, number time, int si) const
		{
			const number x = pos[0];
			const number y = pos[1];

			val[0] = nu*(y - cx);
			val[1] = nu*(cy - x);
		}

	protected:
		double cx, cy, nu;
};

}
#endif
