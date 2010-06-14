/*
 * rhs.h
 *
 *  Created on: 04.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__RHS__
#define __H__LIBDISCRETIZATION__RHS__

#include "lib_grid/lib_grid.h"
#include "numericalsolution.h"
#include "common/math/ugmath.h"
#include <string>

namespace ug {

template <int d>
class RHS {

	public:
		RHS(std::string name);

		void set_name(std::string name);
		std::string name();

		virtual void compute_defect_at_ip(MathVector<d> GlobIP, double& DefectValue) = 0;

		virtual ~RHS()
		{}

	protected:
		std::string m_name;
};

template <int d>
class Scalarf : public RHS<d> {

	protected:
		typedef bool (*RHSFunction)(MathVector<d>, number&);

	public:
		Scalarf(std::string str, RHSFunction rhs) : RHS<d>(str)
		{
			m_RHSFunction = rhs;
		};

		void compute_defect_at_ip(MathVector<d> GlobIP, double& DefectValue)
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

#include "rhs_impl.h"

#endif /* __H__LIBDISCRETIZATION__RHS__ */
