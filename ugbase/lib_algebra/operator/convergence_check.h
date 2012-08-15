/*
 * convergence_check.h
 *
 *  Created on: 18.07.2010
 *      Author: andreasvogel
 *
 *  Comment: Based on an idea by Arne Naegel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK__
#define __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK__

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
#include "lib_algebra/operator/interface/function_base.h"
#include "lib_disc/dof_manager/surface_dof_distribution.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/common/function_group.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Convergence check
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/** IConvergenceCheck
 *
 * This is the base class for a convergence checking object. An instance is
 * passed to an iterative solver to control the convergence.
 *
 */
class IConvergenceCheck
{
	public:
		///////////////////
		// defect control
		///////////////////

		/// sets the given start defect
		virtual void start_defect(number defect) = 0;

		/// computes the start defect and set it
		virtual void start(IFunctionBase& d) = 0;

		/// sets the update for the current defect
		virtual void update_defect(number defect) = 0;

		/// computes the defect and sets it a the next defect value
		virtual void update(IFunctionBase& d) = 0;

		/** iteration_ended
		 *
		 *	Checks if the iteration must be ended.
		 *	This can be due to convergence or divergence.
		 *
		 * \return 	true 		if iteration ended
		 * 			false 		if iteration can and must be continued until convergence
		 */
		virtual bool iteration_ended() = 0;

		/** post
		 *
		 * post-processes the iteration. Some informative outputs of the status of
		 * the iteration after finishing the iteration can be placed here
		 *
		 * \return 	true 		if iteration was successful
		 * 			false 		if iteration did not lead to a satisfying result
		 */
		virtual bool post() = 0;

		/////////////////////////////////////
		// informations about current status
		/////////////////////////////////////

		/// returns the current defect
		virtual number defect() const = 0;

		/// returns the current number of steps
		virtual int step() const = 0;

		// returns the current relative reduction
		virtual number reduction() const = 0;

		////////////////
		// output style
		////////////////

		/// sets the number of spaces printed before output information
		virtual void set_offset(int offset) = 0;

		/// get the current offset
		virtual int get_offset() const = 0;

		/// sets the symbol used for output
		virtual void set_symbol(char symbol) = 0;

		/// sets the name of the iteration
		virtual void set_name(std::string name) = 0;

		/// sets info string
		virtual void set_info(std::string name) = 0;

		/// virtual destructor
		virtual ~IConvergenceCheck() {};
};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Standard convergence check
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

/** StdConvCheck
 *
 * This is a standard implementation of the convergence check.
 * The function_type must provide the following function:
 *  - number two_norm()  --- giving the two norm of the function
 */
class StandardConvCheck : public IConvergenceCheck
{
	public:
		StandardConvCheck();

		StandardConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose);

		void set_verbose(bool level) {m_verbose = level;}
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}
		void set_minimum_defect(number minDefect) {m_minDefect = minDefect;}
		void set_reduction(number relReduction) {m_relReduction = relReduction;}

		void start_defect(number initialDefect);

		void start(IFunctionBase& d);

		void update_defect(number newDefect);

		void update(IFunctionBase& d);

		bool iteration_ended();

		bool post();

		number reduction() const {return m_currentDefect/m_initialDefect;};
		number defect() const {return m_currentDefect;};
		number previous_defect() const { return m_lastDefect; }
		int step() const {return m_currentStep;}

		int get_offset() const {return m_offset;}
		void set_offset(int offset){m_offset = offset;}
		void set_symbol(char symbol){m_symbol = symbol;}
		void set_name(std::string name) {m_name = name;}
		void set_info(std::string info) {m_info = info;}
		const std::vector<number> get_defects() const { return _defects;}

	protected:
		void print_offset();

		bool is_valid_number(number value);

	protected:
		// start defect
		number m_initialDefect;

		// current defect
		number m_currentDefect;

		// defect of the previous step
		number m_lastDefect;

		// current step
		int m_currentStep;

	protected:
		// maximum number of steps to be performed
		int m_maxSteps;

		// absolute reduction to be reached for convergence
		number m_minDefect;

		// relative reduction to be reached for convergence
		number m_relReduction;

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
		std::vector<number> _defects;
};





///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// Individual function convergence check
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

/** IndivFctConvCheck
 *
 * This is an implementation of the convergence check interface,
 * that makes it possible to define required defect reductions on
 * the individual functions constituting the overall solution.
 */
template <class TDomain>
class IndivFctConvCheck : public IConvergenceCheck
{
	public:
	/// constructors
		IndivFctConvCheck();
		//IndivFctConvCheck(ConstSmartPtr<SurfaceDoFDistribution> dofDistr);

	/// definitions for check
		// inherited from IConvergenceCheck
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}
		void set_minimum_defect(number minDefect);
		void set_reduction(number relReduction);

		// own
		void set_approxSpace(const ApproximationSpace<TDomain>& approx);
		void set_functions(const char* functionNames);
		void set_minimum_defect(const char* minDefect, number minDefectForRest = 1e-10);
		void set_reduction(const char* reduction, number reductionForRest = 1e-8);

	/// defect control
		// inherited from IConvergenceCheck
		void start_defect(number initialDefect);
		void start(IFunctionBase& d);
		void update_defect(number newDefect);
		void update(IFunctionBase& d);
		bool iteration_ended();
		bool post();

		// own

	/// information about current status
		// inherited from IConvergenceCheck
		number defect() const {return m_currentOverallDefect;};
		int step() const {return m_currentStep;}
		number reduction() const {return m_currentOverallDefect/m_initialOverallDefect;};

		// own
		number defect(size_t fctIndex) const {return m_currentDefect[fctIndex];};
		number reduction(size_t fctIndex) const {return m_currentDefect[fctIndex]/m_initialDefect[fctIndex];};
		number previousDefect(size_t fctIndex) const {return m_lastDefect[fctIndex];};

	/// output
		// inherited from IConvergenceCheck
		int get_offset() const {return m_offset;};
		void set_offset(int offset){m_offset = offset;};
		void set_symbol(char symbol){m_symbol = symbol;};
		void set_name(std::string name) {m_name = name;};
		void set_info(std::string info) {m_info = info;};

		// own
		void set_verbose(bool level) {m_verbose = level;};
		//const std::vector<number> get_defects() const { return _defects;}

	protected:
		void print_offset();
		bool is_valid_number(number value);
		std::string fctName(size_t fctIndex) {return m_fctName[fctIndex];};

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
		bool m_locked;
		std::vector<std::vector<size_t> > m_indices;
		std::vector<std::string> m_fctName;
		FunctionGroup m_fctGrp;
		ConstSmartPtr<SurfaceDoFDistribution> m_dd;
};


} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK__ */
