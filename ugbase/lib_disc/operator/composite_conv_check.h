/*
 * composite_conv_check.h
 *
 *      Author: M. Breit
 */

#ifndef __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK__
#define __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK__

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
	/// \{
		CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > approx);
		CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		                   int maxSteps, number minDefect, number relReduction);
	/// \}

	/// set level of grid, where defect vectors come from
		void set_level(int level);

	/// sets maximum number of iteration steps
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}

	///	sets default values for non-explicitly specified cmps
		void set_rest_check(number minDefect, number relReduction){
			m_bCheckRest = true;
			m_restMinDefect = minDefect; m_restRelReduction = relReduction;
			update_rest_check();
		}

	///	disables rest check
		void disable_rest_check(){m_bCheckRest = false; update_rest_check();}

	///	sets check for single component
	/// \{
		void set_component_check(const std::string& vFctName,
		                         const std::vector<number>& vMinDefect,
		                         const std::vector<number>& vRelReduction);

		void set_component_check(const std::vector<std::string>& vFctName,
		                         const std::vector<number>& vMinDefect,
		                         const std::vector<number>& vRelReduction);

		void set_component_check(const std::vector<std::string>& vFctName,
								 const number minDefect,
								 const number relReduction);

		void set_component_check(const std::string& fctName,
								 const number minDefect,
								 const number relReduction);
	/// \}

	///	sets check for all components in approximation space
		void set_all_component_check(const number minDefect,
		                             const number relReduction);

	///	sets check for group of components
	/// \{
		void set_group_check(const std::vector<std::string>& vFctName,
							 const number minDefect,
							 const number relReduction);

		void set_group_check(const std::string& fctNames,
							 const number minDefect,
							 const number relReduction);
	/// \}

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

	///	sets if verbose
		void set_verbose(bool level) {m_verbose = level;};

	///	enables time measurement
		void timeMeasurement(bool yesOrNo) {m_bTimeMeas = yesOrNo;};

	///	clones this instance
		virtual SmartPtr<IConvergenceCheck<TVector> > clone();

	///	prints a line using prefixes
		void print_line(std::string line);

	protected:
		void print_offset();
		bool is_valid_number(number value);
		const std::string& fctName(size_t i) {return m_CmpInfo[i].name;};

	///	extracts multi-indices for a fct-comp on a element type
		template <typename TBaseElem>
		void extract_multi_indices(ConstSmartPtr<DoFDistribution> dd);

	///	extracts multi-indices from dof distribution
		void extract_multi_indices(ConstSmartPtr<DoFDistribution> dd);

	/// calculates the 2-norm of the entries of the vector vec specified by index
		number norm(const TVector& vec, const std::vector<MultiIndex<2> >& index);

	protected:
	///	ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		struct NativCmpInfo{
			std::string name; 	///< Name of components

			number initDefect;	///< Initial Defect of component
			number currDefect;	///< Current Defect of component
			number lastDefect;	///< Last Defect if component

			std::vector<MultiIndex<2> > vMultiIndex; ///< associated indices
		};

	///	info on natural components
		std::vector<NativCmpInfo> m_vNativCmpInfo;
		size_t m_numAllDoFs;

	///	returns defect for all components
		number defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vNativCmpInfo.size(); ++fct)
				defect += pow(m_vNativCmpInfo[fct].currDefect, 2);
			return sqrt(defect);
		}

	///	returns last defect for all components
		number last_defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vNativCmpInfo.size(); ++fct)
				defect += pow(m_vNativCmpInfo[fct].lastDefect, 2);
			return sqrt(defect);
		}

	///	returns initial defect for all components
		number initial_defect_all() const {
			number defect = 0.0;
			for(size_t fct = 0; fct < m_vNativCmpInfo.size(); ++fct)
				defect += pow(m_vNativCmpInfo[fct].initDefect, 2);
			return sqrt(defect);
		}

	protected:
		struct CmpInfo
		{
			CmpInfo() : isRest(false) {}

			std::vector<int> vFct;	///< Index of components
			std::string name; 		///< Name of components

			number initDefect;	///< Initial Defect of component
			number currDefect;	///< Current Defect of component
			number lastDefect;	///< Last Defect if component

			number minDefect;	///< Minimal required Defect of component
			number relReduction;///< Relative reduction required for component

			bool isRest; 	///< Shows, that this is group of remaining cmps
		};

	///	infos for each component
		std::vector<CmpInfo> m_CmpInfo;

	///	default Values
		bool m_bCheckRest;
		number m_restMinDefect;
		number m_restRelReduction;
		void update_rest_check();

	/// current step
		int m_currentStep;

	/// maximum number of steps to be performed
		int m_maxSteps;

	protected:
		/// verbose level
		bool m_verbose;

		/// number of spaces inserted before output
		int m_offset;

		/// symbol for output appearance
		char m_symbol;

		/// name of iteration
		std::string m_name;

		/// info for iteration (e.g. preconditioner type)
		std::string m_info;

	protected:
		/// enables time measurement
		bool m_bTimeMeas;

		/// a stopwatch
		Stopwatch m_stopwatch;
};

} // end namespace ug

#include "composite_conv_check_impl.h"

#endif /* __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK__ */
