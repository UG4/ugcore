/*
 * rhs.h
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__RHS__
#define __H__LIBDISCRETIZATION__RHS__

#include "lib_grid/lib_grid.h"
#include "lib_grid/geometric_objects.h"
#include "numericalsolution.h"
#include "common/math/ugmath.h"
#include <string>

namespace ug {

class RHS {

	public:
		RHS(std::string name);

		void set_name(std::string name);
		std::string name();

		virtual void compute_defect_at_ip(vector2 LocIP, vector3 GlobIP, double& DefectValue) = 0;

		virtual ~RHS()
		{}

	protected:
		std::string m_name;
};

class Scalarf : public RHS {

	protected:
		typedef bool (*RHSFunction)(MathVector<3>, number&);

	public:
		Scalarf(std::string str, RHSFunction rhs):RHS(str)
		{
			m_RHSFunction = rhs;
		};

		void compute_defect_at_ip(MathVector<2> LocIP, MathVector<3> GlobIP, double& DefectValue)
		{
			number f;

			m_RHSFunction(GlobIP, f);

			DefectValue = f;
		}

		~Scalarf()
		{}

	protected:
		RHSFunction m_RHSFunction;

};

} /* end namespace libDiscretization */



#endif /* __H__LIBDISCRETIZATION__RHS__ */
