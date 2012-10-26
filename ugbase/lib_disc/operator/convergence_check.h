/*
 * convergence_check.h
 *
 *      Author: M. Breit
 */

#ifndef __H__LIB_DISC__OPERATOR__CONVERGENCE_CHECK__
#define __H__LIB_DISC__OPERATOR__CONVERGENCE_CHECK__

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

/** CompositeConvCheck
 *
 * This is an implementation of the convergence check interface,
 * that makes it possible to define required defect reductions on
 * the individual functions constituting the overall solution.
 */
template <class TVector, class TDomain>
class CompositeConvCheck : public IConvergenceCheck<TVector>
{
	public:
	/// constructors
		CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > approx);

	/// sets maximum number of iteration steps
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}

	///	sets the approximation space
		void set_functions(const char* functionNames);
		void set_minimum_defect(const char* minDefect, number minDefectForRest = 1e-10);
		void set_reduction(const char* reduction, number reductionForRest = 1e-8);

	/// defect control
		void start_defect(number initialDefect);
		void start(const TVector& d);
		void update_defect(number newDefect);
		void update(const TVector& d);
		bool iteration_ended();
		bool post();

	/// information about current status
		number defect() const {return m_currentOverallDefect;};
		int step() const {return m_currentStep;}
		number reduction() const {return m_currentOverallDefect/m_initialOverallDefect;};

	/// information about current status for single component
		number defect(size_t fctIndex) const {return m_currentDefect[fctIndex];};
		number reduction(size_t fctIndex) const {return m_currentDefect[fctIndex]/m_initialDefect[fctIndex];};
		number previousDefect(size_t fctIndex) const {return m_lastDefect[fctIndex];};

	/// output
		int get_offset() const {return m_offset;};
		void set_offset(int offset){m_offset = offset;};
		void set_symbol(char symbol){m_symbol = symbol;};
		void set_name(std::string name) {m_name = name;};
		void set_info(std::string info) {m_info = info;};

	///	sets if verbose
		void set_verbose(bool level) {m_verbose = level;};

	///	enables time measurement
		void timeMeasurement(bool yesOrNo) {m_timeMeas = yesOrNo;};

	protected:
		void print_offset();
		bool is_valid_number(number value);
		std::string fctName(size_t fctIndex) {return m_fctName[fctIndex];};

	///	extracts multi-indices for a fct-comp on a element type
		template <typename TBaseElem>
		void extract_multi_indices();

	/// calculates the 2-norm of the entries of the vector vec specified by index
		number norm(const TVector& vec, std::vector<MultiIndex<2> > index);

	protected:
		// start defect
		std::vector<number> m_initialDefect;
		number m_initialOverallDefect;

		// current defect
		std::vector<number> m_currentDefect;
		number m_currentOverallDefect;

		// defect of the previous step
		std::vector<number> m_lastDefect;

		// current step
		int m_currentStep;

		// maximum number of steps to be performed
		int m_maxSteps;

		// absolute reduction to be reached for convergence
		std::vector<number> m_minDefect;

		// relative reduction to be reached for convergence
		std::vector<number> m_relReduction;

	protected:
		// verbose level
		bool m_verbose;

		// number of spaces inserted before output
		int m_offset;

		// symbol for output appearance
		char m_symbol;

		// name of iteration
		std::string m_name;

		// info for iteration (e.g. preconditioner type)
		std::string m_info;

	private:
		bool m_timeMeas;
		Stopwatch m_stopwatch;
		std::vector<std::vector<MultiIndex<2> > > m_vvMultiIndex;
		std::vector<std::string> m_fctName;
		FunctionGroup m_fctGrp;
		ConstSmartPtr<SurfaceDoFDistribution> m_dd;
};

} // end namespace ug

#include "convergence_check_impl.h"

#endif /* __H__LIB_DISC__OPERATOR__CONVERGENCE_CHECK__ */
