/*
 * convergence_check.cpp
 *
 *  Created on: 13.03.2012
 *      Author: andreasvogel
 */

#include "convergence_check.h"
#include "common/util/string_util.h"

namespace ug{


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Standard convergence check							//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

StandardConvCheck::
StandardConvCheck()
 :	 m_initialDefect(0.0), m_currentDefect(0.0), m_lastDefect(0.0), m_currentStep(0),
	 m_maxSteps(200), m_minDefect(10e-8), m_relReduction(10e-10),
	 m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info("")
	 {};

StandardConvCheck::
StandardConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose)
 :	 m_initialDefect(0.0), m_currentDefect(0.0), m_lastDefect(0.0), m_currentStep(0),
	 m_maxSteps(maxSteps), m_minDefect(minDefect), m_relReduction(relReduction),
	 m_verbose(verbose), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info("")
	 {};

void StandardConvCheck::start_defect(number initialDefect)
{
	_defects.clear();
	m_initialDefect = initialDefect;
	m_currentDefect = m_initialDefect;
	m_currentStep = 0;

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

void StandardConvCheck::start(IFunctionBase& d)
{
	start_defect(d.two_norm());
}

void StandardConvCheck::update_defect(number newDefect)
{
	m_lastDefect = m_currentDefect;
	m_currentDefect = newDefect;
	m_currentStep++;

	if(m_verbose)
	{
		print_offset(); UG_LOG(std::setw(4) << step() << ":    " << std::scientific << defect() <<
							"    " << defect()/m_lastDefect << "\n");
		_defects.push_back(defect());
	}
}

void StandardConvCheck::update(IFunctionBase& d)
{
	update_defect(d.two_norm());
}

bool StandardConvCheck::iteration_ended()
{
	if(!is_valid_number(m_currentDefect)) return true;
	if(step() >= m_maxSteps) return true;
	if(defect() < m_minDefect) return true;
	if(reduction() < m_relReduction) return true;
	return false;
}

bool StandardConvCheck::post()
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

	if(!success)
	{
		if (!is_valid_number(m_currentDefect))
			if(m_verbose)
			{
				print_offset(); UG_LOG("Current defect << " << m_currentDefect << " is not a valid number.\n");
			}

		if(step() >= m_maxSteps)
			if(m_verbose)
			{
				print_offset(); UG_LOG("Maximum numbers of "<< m_maxSteps << " iterations reached without convergence.\n");
			};
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

void StandardConvCheck::print_offset()
{
	// step 1: whitespace
	UG_LOG(repeat(' ', m_offset));

	// step 2: print style character
	UG_LOG(m_symbol << " ");
}

bool StandardConvCheck::is_valid_number(number value)
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
