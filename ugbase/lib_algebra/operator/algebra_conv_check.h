/*
 * algebra_conv_check.h
 *
 * A. Naegel
 * (algebraic version of composite_conv_check.h by M. Breit)
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
#include <math.h>

#include "common/common.h"
#include "common/stopwatch.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_disc/dof_manager/surface_dof_distribution.h"
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
template <class TVector>
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

#endif /* __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK__ */
