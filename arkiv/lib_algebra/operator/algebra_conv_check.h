/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

/*
 * algebraic version of composite_conv_check.h by M. Breit
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_CONVERGENCE_CHECK__
#define __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_CONVERGENCE_CHECK__

#include <ostream>
#include <sstream>
#include <string>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>

#include "common/common.h"
#include "common/stopwatch.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/common/function_group.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Composite convergence check
////////////////////////////////////////////////////////////////////////////////

/** AlgebraicConvCheck
 *
 * This is an implementation of the convergence check interface,
 * that makes it possible to define required defect reductions on
 * the individual functions constituting the overall solution.
 */
template <typename TVector>
class AlgebraicConvCheck : public IConvergenceCheck<TVector>
{
	public:
	/// constructors
	/// \{
		AlgebraicConvCheck(size_t ncmp);
		AlgebraicConvCheck(size_t ncmp, int maxSteps, number minDefect, number relReduction);
		AlgebraicConvCheck(size_t ncmp, int maxSteps, number minDefect, number relReduction, bool verbose);
	/// \}

		///	clones this instance
		virtual SmartPtr<IConvergenceCheck<TVector> > clone();
		virtual std::string config_string() const {return std::string("AlgebraicConvCheck");}

		/// defect control
		void start_defect(number initialDefect);
		void start(const TVector& d);
		void update_defect(number newDefect);
		void update(const TVector& d);
		bool iteration_ended();
		bool post();

		/// information about current status
		int step() const {return m_currentStep;}
		number defect() const {return defect_all();};
		number reduction() const {return defect_all()/initial_defect_all();};
		number rate() const {return defect_all()/last_defect_all();}
		number avg_rate() const {return std::pow(defect_all()/initial_defect_all(), 1.0/m_currentStep);}

		/// output
		int get_offset() const {return m_offset;};
		void set_offset(int offset){m_offset = offset;};
		void set_symbol(char symbol){m_symbol = symbol;};
		void set_name(std::string name) {m_name = name;};
		void set_info(std::string info) {m_info = info;};


	/// sets maximum number of iteration steps
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}

	///	sets check for single component
		inline void set_component_check(const size_t cmp,
								 const number abs,
								 const number red)
		{
			m_vCmpInfo[cmp].minDefect = abs;
			m_vCmpInfo[cmp].relReduction = red;
		}

	///	sets check for all components
		void set_component_checks(const number abs, const number red)
		{
			for(size_t cmp = 0; cmp < m_vCmpInfo.size(); ++cmp)
				set_component_check(cmp, abs, red);
		}

	///	sets if verbose
		void set_verbose(bool level) {m_verbose = level;};

	///	enables time measurement
		void set_time_measurement(bool yesOrNo) {m_bTimeMeas = yesOrNo;};

	///	prints a line using prefixes
		void print_line(std::string line);


	/// statistics
		void get_statistics(double *first, double *last, int &niter) const
		{
			for (size_t cmp = 0; cmp < m_vCmpInfo.size(); cmp++)
			{
				const CmpInfo& cmpInfo = m_vCmpInfo[cmp];
				first[cmp] = cmpInfo.initDefect;
				last[cmp] = cmpInfo.currDefect;
			}
			niter = step();
		};
	protected:
		void print_offset();
		bool is_valid_number(number value);
		const std::string& fctName(size_t i) {return m_vCmpInfo[i].name;};

	/// calculates the 2-norm of the entries of the vector vec specified by index
		number norm(const TVector& vec, size_t cmp);

	protected:

		struct CmpInfo{
			CmpInfo(number minDef, number relRed) : minDefect(minDef), relReduction(relRed) {}
			std::string name; 	///< Name of components

			number initDefect;	///< Initial Defect of component
			number currDefect;	///< Current Defect of component
			number lastDefect;	///< Last Defect if component

			number minDefect;	///< Minimal required Defect of component
			number relReduction;///< Relative reduction required for component
			number weight; ///< weight for this component
		};

	///	info on components
		std::vector<CmpInfo> m_vCmpInfo;





	///	returns defect for all components
		number defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vCmpInfo.size(); ++fct)
				defect += (m_vCmpInfo[fct].currDefect*m_vCmpInfo[fct].currDefect);
			return sqrt(defect);
		}

	///	returns last defect for all components
		number last_defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vCmpInfo.size(); ++fct)
				defect += (m_vCmpInfo[fct].lastDefect*m_vCmpInfo[fct].lastDefect);
			return sqrt(defect);
		}

	///	returns initial defect for all components
		number initial_defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vCmpInfo.size(); ++fct)
				defect += m_vCmpInfo[fct].initDefect*m_vCmpInfo[fct].initDefect;
			return sqrt(defect);
		}

	protected:


		int m_maxSteps;///< maximum number of steps to be performed
		number m_minDefect;	///< Minimal required Defect of component
		number m_relReduction;///< Relative reduction required for component

		bool m_verbose;	///< verbose level
		int m_currentStep;///< current step


		int m_offset;///< number of spaces inserted before output
		char m_symbol;///< symbol for output appearance
		std::string m_name;///< name of iteration
		std::string m_info;///< info for iteration (e.g. preconditioner type)

		bool m_bTimeMeas;///< enables time measurement
		Stopwatch m_stopwatch;///< a stopwatch
};

} // end namespace ug

#include "algebra_conv_check_impl.h"

#endif