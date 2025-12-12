/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

/*
 *  Comment: Based on an idea by Arne Naegel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK__
#define __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK__

#include <ostream>
#include <sstream>
#include <string>
//#include <limits>
//#include <algorithm>
#include <iomanip>
#include <vector>
//#include <list>
#include <cmath>

//#include "common/common.h"
#include "common/stopwatch.h"
#include "lib_disc/function_spaces/approximation_space.h"
//#include "lib_disc/common/function_group.h"

namespace ug {

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
 * \tparam TVector 	vector type
 */
template <typename TVector>
class IConvergenceCheck
{
	public:
		/// sets the given start defect
		virtual void start_defect(number defect) = 0;

		/// computes the start defect and set it
		virtual void start(const TVector& d) = 0;

		/// sets the update for the current defect
		virtual void update_defect(number defect) = 0;

		/// computes the defect and sets it a the next defect value
		virtual void update(const TVector& d) = 0;

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

		// returns the current convergence rate
		virtual number rate() const = 0;

		// returns the averaged convergence rate
		virtual number avg_rate() const = 0;

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

		/// prints a line
		virtual void print_line(std::string line) = 0;

		/// clone the object
		virtual SmartPtr<IConvergenceCheck<TVector> > clone() = 0;

		/// virtual destructor
		virtual ~IConvergenceCheck() = default;

		///	returns information about configuration parameters
		/**
		 * this should return necessary information about parameters and possibly
		 * calling config_string of subcomponents.
		 *
		 * \returns std::string	necessary information about configuration parameters
		 */

		virtual std::string config_string() const = 0;
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
 *  - number norm()  --- giving the two norm of the function
 */
template <typename TVector>
class StdConvCheck : public IConvergenceCheck<TVector>
{
	public:
		StdConvCheck();

		StdConvCheck(int maxSteps, number minDefect, number relReduction);
		StdConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose);
		StdConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose,bool suppressUnsuccessful);

		void set_verbose(bool level) {m_verbose = level;}
		void set_maximum_steps(int maxSteps) {m_maxSteps = maxSteps;}
		void set_minimum_defect(number minDefect) {m_minDefect = minDefect;}
		void set_reduction(number relReduction) {m_relReduction = relReduction;}
		void set_supress_unsuccessful(bool bsupress){ m_supress_unsuccessful = bsupress; }

		void start_defect(number initialDefect) override;

		void start(const TVector& d) override;

		void update_defect(number newDefect) override;

		void update(const TVector& d) override;

		bool iteration_ended() override;

		bool post() override;

		std::string config_string() const override
		{
			std::stringstream ss;
			ss << "StdConvCheck( max steps = " << m_maxSteps << ", min defect = " << m_minDefect << ", relative reduction = " << m_relReduction << ")";
			return ss.str();
		}

		number reduction() const override {return m_currentDefect/m_initialDefect;};
		number defect() const override {return m_currentDefect;};
		number previous_defect() const { return m_lastDefect; }
		int step() const override {return m_currentStep;}
		number rate() const override {return m_currentDefect/m_lastDefect;};
		number avg_rate() const override {return std::pow((number)m_ratesProduct,(number)1.0/step());}

		int get_offset() const override {return m_offset;}
		void set_offset(int offset) override {m_offset = offset;}
		void set_symbol(char symbol) override {m_symbol = symbol;}
		void set_name(std::string name) override {m_name = name;}
		void set_info(std::string info) override {m_info = info;}
		const std::vector<number> get_defects() const { return _defects;}
		number get_defect(size_t i) const { return _defects[i];}
		void print_line(std::string line) override;

		SmartPtr<IConvergenceCheck<TVector> > clone() override {
			SmartPtr<StdConvCheck > newInst(new StdConvCheck);
			// use std assignment (implicit member-wise is fine here)
			*newInst = *this;
			return newInst;
		}

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

		// sum of the convergence rates over all steps
		number m_ratesProduct;

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
		
		// return true in post method if max nr of iterations is reached
		bool m_supress_unsuccessful;

	private:
		std::vector<number> _defects;
};

} // end namespace ug

#include "convergence_check_impl.h"

#endif