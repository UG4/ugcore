/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK_IMPL__
#define __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK_IMPL__

#include "convergence_check.h"
#include "common/util/string_util.h"

namespace ug{


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Standard convergence check							//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

template <typename TVector>
StdConvCheck<TVector>::
StdConvCheck()
 :	 m_initialDefect(0.0), m_currentDefect(0.0), m_lastDefect(0.0), m_currentStep(0),
  	 m_ratesProduct(1), m_maxSteps(200), m_minDefect(10e-8), m_relReduction(10e-10),
	 m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	 m_supress_unsuccessful(false)
	 {};

template <typename TVector>
StdConvCheck<TVector>::
StdConvCheck(int maxSteps, number minDefect, number relReduction)
 :	 m_initialDefect(0.0), m_currentDefect(0.0), m_lastDefect(0.0), m_currentStep(0),
  	 m_ratesProduct(1), m_maxSteps(maxSteps), m_minDefect(minDefect), m_relReduction(relReduction),
	 m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	 m_supress_unsuccessful(false)
	 {};

template <typename TVector>
StdConvCheck<TVector>::
StdConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose)
 :	 m_initialDefect(0.0), m_currentDefect(0.0), m_lastDefect(0.0), m_currentStep(0),
  	 m_ratesProduct(1), m_maxSteps(maxSteps), m_minDefect(minDefect), m_relReduction(relReduction),
	 m_verbose(verbose), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	 m_supress_unsuccessful(false)
	 {};

template <typename TVector>
StdConvCheck<TVector>::
StdConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose,bool supressUnsuccessful)
 :	 m_initialDefect(0.0), m_currentDefect(0.0), m_lastDefect(0.0), m_currentStep(0),
  	 m_ratesProduct(1), m_maxSteps(maxSteps), m_minDefect(minDefect), m_relReduction(relReduction),
	 m_verbose(verbose), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	 m_supress_unsuccessful(supressUnsuccessful)
	 {};

template <typename TVector>
void StdConvCheck<TVector>::start_defect(number initialDefect)
{
	//UG_LOG_ALL_PROCS("start defect is " << initialDefect << "\n")
	_defects.clear();
	m_initialDefect = initialDefect;
	m_currentDefect = m_initialDefect;
	m_currentStep = 0;
	m_ratesProduct = 1;

	if(m_verbose)
	{
		UG_LOG("\n");

	//  number of symbols to print before name and info
		int num_sym = 8;
		int num_line_length = 50;

		int max_length = std::max(m_name.length(), m_info.length());
		int space_left = std::max(num_line_length - max_length - num_sym, 0);

	//	print name line
		print_offset();
		UG_LOG(repeat(m_symbol, num_sym));
		int pre_space = (int)(max_length -(int)m_name.length()) / 2;
		UG_LOG(repeat(' ', pre_space));
		UG_LOG("  "<< m_name << "  ");
		UG_LOG(repeat(' ', max_length - pre_space -m_name.length()));
		UG_LOG(repeat(m_symbol, space_left));
		UG_LOG("\n");
	//	print info line
		print_offset();
		if(m_info.length() > 0)
		{
			UG_LOG(repeat(m_symbol, num_sym));
			UG_LOG("  "<< m_info << "  ");
			UG_LOG(repeat(' ', max_length-m_info.length()));
			UG_LOG(repeat(m_symbol, space_left))
			UG_LOG("\n");
		} else {
			UG_LOG("\n");
		}

	//	start iteration output
		print_offset(); UG_LOG("  Iter      Defect         Rate \n");
		print_offset(); UG_LOG(std::setw(4) << step() << ":    "
								<< std::scientific << defect() <<  "      -------\n");
	}
}

template <typename TVector>
void StdConvCheck<TVector>::start(const TVector& d)
{
	start_defect(d.norm());
}

template <typename TVector>
void StdConvCheck<TVector>::update_defect(number newDefect)
{
	m_lastDefect = m_currentDefect;
	m_currentDefect = newDefect;
	m_currentStep++;
	m_ratesProduct *= newDefect/m_lastDefect;

	if(m_verbose)
	{
		print_offset(); UG_LOG(std::setw(4) << step() << ":    " << std::scientific << defect() <<
							"    " << defect()/m_lastDefect << "\n");
		_defects.push_back(defect());
	}
}

template <typename TVector>
void StdConvCheck<TVector>::update(const TVector& d)
{
	update_defect(d.norm());
}

template <typename TVector>
bool StdConvCheck<TVector>::iteration_ended()
{
	if(!is_valid_number(m_currentDefect)) return true;
	if(step() >= m_maxSteps) return true;
	if(defect() < m_minDefect) return true;
	if(reduction() < m_relReduction) return true;
	return false;
}

template <typename TVector>
bool StdConvCheck<TVector>::post()
{
	bool success = false;

	if(defect() < m_minDefect)
	{
		if(m_verbose)
		{
			print_offset(); UG_LOG("Absolute defect norm " << m_minDefect << " reached after " << step() << " steps.\n");
		}
		success = true;
	};

	if(reduction() < m_relReduction)
	{
		if(m_verbose)
		{
			print_offset(); UG_LOG("Relative reduction " << m_relReduction << " reached after " << step() << " steps.\n");
		}
		success = true;
	};

	if (m_verbose && is_valid_number(m_currentDefect))
	{
		print_offset(); UG_LOG("Average reduction over " << step() << " steps: " << pow(reduction(), 1.0/step()) << "\n");
	}

	if(!success)
	{
		if (!is_valid_number(m_currentDefect))
			if(m_verbose)
			{
				print_offset(); UG_LOG("Current defect " << m_currentDefect << " is not a valid number.\n");
			}

		if(step() >= m_maxSteps){
			if(m_verbose)
			{
				print_offset(); UG_LOG("Maximum numbers of "<< m_maxSteps << " iterations reached without convergence.\n");
			}
			if (m_supress_unsuccessful) return true;
		}
	}

	if(m_verbose)
	{
		print_offset();
		UG_LOG(repeat(m_symbol, 5));
		if(success) {UG_LOG("  Iteration converged  ");}
		else {UG_LOG("  Iteration not successful  ");}
		UG_LOG(repeat(m_symbol, 5));
		UG_LOG("\n\n");
	}
	return success;
}

template <typename TVector>
void StdConvCheck<TVector>::print_offset()
{
	// step 1: whitespace
	UG_LOG(repeat(' ', m_offset));

	// step 2: print style character
	UG_LOG(m_symbol << " ");
}

template <typename TVector>
void StdConvCheck<TVector>::print_line(std::string line)
{
	print_offset();
	UG_LOG(line << "\n");
}


template <typename TVector>
bool StdConvCheck<TVector>::is_valid_number(number value)
{
	// (value >= std::numeric_limits<number>::min() ) == true if value > -infty
	// (value <= std::numeric_limits<number>::max() ) == true if value < infty
	// (value == value                         ) == true if value != NaN

	if (value == 0.0) return true;
	else return (value >= std::numeric_limits<number>::min()
				&& value <= std::numeric_limits<number>::max()
				&& value == value && value >= 0.0);
}

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__CONVERGENCE_CHECK_IMPL__ */
