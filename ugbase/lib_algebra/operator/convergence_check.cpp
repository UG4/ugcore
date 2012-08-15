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




//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Individual function convergence check				//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
template <class TDomain>
IndivFctConvCheck<TDomain>::
IndivFctConvCheck()
 :	 m_initialDefect(0), m_initialOverallDefect(0.0),
  	 m_currentDefect(0), m_currentOverallDefect(0.0),
 	 m_lastDefect(0), m_currentStep(0), m_maxSteps(100),
 	 m_minDefect(0),
 	 m_relReduction(0),
	 m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	 m_locked(false),
	 m_dd(NULL)
{};

template <class TDomain>
void IndivFctConvCheck<TDomain>::set_approxSpace(const ApproximationSpace<TDomain>& approx)
{
	// do not call this method a second time after locking
	if (m_locked) UG_THROW("You can not set the approximation space twice.");

	m_dd = approx.surface_dof_distribution();

	//\TODO: check, if this is correct (a fortiori wrt 'hanging' and block vectors), please!
	// compute indices for faster access later
	m_indices.clear();
	typedef typename SurfaceDoFDistribution::dim_traits<TDomain::dim>::geometric_base_object Elem;
	typename SurfaceDoFDistribution::template traits<Elem>::const_iterator iter, iterBegin, iterEnd;

	for (size_t fi = 0; fi < m_dd->function_pattern().num_fct(); fi++)
	{
		//	get element iterator for current subset
		iterBegin = m_dd->template begin<Elem>();
		iterEnd = m_dd->template end<Elem>();

		// buffer for indices (will contain duplicates)
		std::vector<size_t> buffer;

		// loop over all elements
		for (iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get local indices
			LocalIndices ind;
			m_dd->indices(*iter, ind, true);	// false?

			for (size_t dof = 0; dof < ind.num_dof(fi); dof++)
			buffer.push_back(ind.index(fi,dof));
		}

		// remove duplicates from buffer and store
		std::sort(buffer.begin(), buffer.end());
		std::unique(buffer.begin(),buffer.end());
		m_indices.push_back(buffer);

		// store function name
		m_fctName.push_back(m_dd->function_pattern().name(fi));
	}
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::set_functions(const char* functionNames)
{
	// do not call this method a second time after locking
	if (m_locked) UG_THROW("You can not set the functions twice.");

	// get the functions specified by function names
	try {m_fctGrp = m_dd->fct_grp_by_name(functionNames);}
		UG_CATCH_THROW("At least one of the functions in '" << functionNames
						<< "' is not contained in the approximation space (or something else was wrong).");


	// remove unnecessary function names and indices / add 'rest' if necessary
	std::vector<std::vector<size_t> > finalIndices(0);
	std::vector<std::string> finalNames(0);
	std::vector<bool> used(m_dd->num_fct(), false);
	for (size_t i = 0; i < m_fctGrp.num_fct(); i++)
	{
		finalIndices.push_back(m_indices[m_fctGrp[i]]);
		finalNames.push_back(m_fctName[m_fctGrp[i]]);
		used[m_fctGrp[i]] = true;
	}

	std::vector<size_t> rest(0);
	std::stringstream ss; ss << "rest (";
	std::vector<std::string>::iterator itName = m_fctName.begin();
	std::vector<std::vector<size_t> >::iterator itInd = m_indices.begin();
	for (std::vector<bool>::iterator itUsed = used.begin(); itUsed != used.end(); itUsed++)
	{
		if (!*itUsed)
		{
			if (rest.size()) ss << ", ";
			ss << *itName;
			rest.insert(rest.end(), itInd->begin(), itInd->end());
		}
		itInd++;
		itName++;
	}
	ss << ")";

	if (rest.size())
	{
		finalIndices.push_back(rest);
		finalNames.push_back(ss.str());
	}

	m_indices = finalIndices;
	m_fctName = finalNames;

	m_locked = true;
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::set_minimum_defect(number minDefect)
{
	UG_THROW(	"This method is not intended to be used in this implementation\n"
				"of IConvergenceCheck. If you want to set a minimal absolute defect\n"
				"value for all functions, use:\n"
				"set_minimum_defect(\"\", minVal).");
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::set_reduction(number relReduction)
{
	UG_THROW(	"This method is not intended to be used in this implementation\n"
				"of IConvergenceCheck. If you want to set a defect reduction\n"
				"value for all functions, use:\n"
				"set_reduction(\"\", red).");
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::set_minimum_defect(const char* minDefect, number minDefectForRest)
{
	//	tokenize strings
	std::vector<std::string> tokens;
	TokenizeString(minDefect, tokens, ',');

	// check if number of values is correct
	if (tokens.size() != m_fctGrp.num_fct())
	{
		UG_THROW(	"The number of supplied values (" << tokens.size() << ") does not match the number\n"
					"of given function names (" << m_fctGrp.num_fct() << "); perhaps you have forgot to call\n"
					"IndivFctConvCheck::set_functions prior to this method.");
	}

	// save values as number
	m_minDefect.clear();
	for (size_t i = 0; i < tokens.size(); i++)
	{
		std::istringstream stm;
		m_minDefect.push_back(0.0);
		stm.str(tokens[i]);
		stm >> m_minDefect.back();
	}

	// set minDefectForRest, if needed
	if (m_fctGrp.num_fct() < m_dd->num_fct())
		m_minDefect.push_back(minDefectForRest);
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::set_reduction(const char* reduction, number reductionForRest)
{
	//	tokenize strings
	std::vector<std::string> tokens;
	TokenizeString(reduction, tokens, ',');

	// check if number of values is correct
	if (tokens.size() != m_fctGrp.num_fct())
	{
		UG_THROW(	"The number of supplied values (" << tokens.size() << ") does not match the number\n"
					"of given function names (" << m_fctGrp.num_fct() << "); perhaps you have forgot to call\n"
					"IndivFctConvCheck::set_functions prior to this method.");
	}

	// save values as number
	m_relReduction.clear();
	for (size_t i = 0; i < tokens.size(); i++)
	{
		std::istringstream stm;
		m_relReduction.push_back(0.0);
		stm.str(tokens[i]);
		stm >> m_relReduction.back();
	}

	// set minDefectForRest, if needed
	if (m_fctGrp.num_fct() < m_dd->num_fct())
		m_relReduction.push_back(reductionForRest);
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::start_defect(number initialDefect)
{
	UG_THROW(	"This method cannot be used to set defect values,\n"
				"since obviously this class is meant for an individual\n"
				"defect calculation of more than one function\n"
				"(use start(IFunctionBase& d) instead).");
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::start(IFunctionBase& d)
{
	m_initialDefect.clear();

	// calculate the defect's 2-norm for each function
	for (size_t i = 0; i < m_indices.size(); i++)
	{
		m_initialDefect.push_back(d.two_norm(m_indices[i]));
		m_initialOverallDefect += m_initialDefect.back()*m_initialDefect.back();
	}
	m_initialOverallDefect = sqrt(m_initialOverallDefect);
	m_currentDefect = m_initialDefect;
	m_currentOverallDefect = m_initialOverallDefect;
	m_currentStep = 0;

	if (m_verbose)
	{
		UG_LOG("\n");

		//  number of symbols to print before name and info
		int num_sym = 18;
		int num_line_length = 80;

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
		if (m_info.length() > 0)
		{
			UG_LOG(repeat(m_symbol, num_sym));
			UG_LOG("  "<< m_info << "  ");
			UG_LOG(repeat(' ', max_length-m_info.length()));
			UG_LOG(repeat(m_symbol, space_left))
			UG_LOG("\n");
		}
		else
		{
			UG_LOG("\n");
		}

	//	start iteration output
		print_offset(); UG_LOG("  Iter      Defect           Rate         Reduction        Required      Function\n");

		print_offset(); UG_LOG(std::setw(4) << step() << ":    "
								<< std::scientific << defect(0) <<  "      --------        -------       "
								<< m_relReduction[0] << "    " << std::setw(0) << fctName(0) << "\n");
		for (size_t i = 1; i < m_fctName.size(); i++)
		{
			print_offset(); UG_LOG("         "
									<< std::scientific << defect(i) <<  "      --------        -------       "
									<< m_relReduction[i] << "    " << fctName(i) << "\n");
		}
	}
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::update_defect(number newDefect)
{
	UG_THROW(	"This method cannot be used to update defect values,\n"
				"since obviously this class is meant for an individual\n"
				"defect calculation of more than one function\n"
				"(use update(IFunctionBase& d) instead).");
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::update(IFunctionBase& d)
{
	m_currentOverallDefect = 0.0;
	m_lastDefect = m_currentDefect;

	// calculate the defect's 2-norm for each function
	for (size_t i = 0; i < m_indices.size(); i++)
	{
		m_currentDefect[i] = d.two_norm(m_indices[i]);
		m_currentOverallDefect += m_currentDefect[i] * m_currentDefect[i];
	}
	m_currentOverallDefect = sqrt(m_currentOverallDefect);

	m_currentStep++;

	if (m_verbose)
	{
		print_offset(); UG_LOG(std::setw(4) << step() << ":    "
										<< std::scientific << defect(0) << "    " << defect(0)/previousDefect(0)
										<< "    " << reduction(0) << "    "
										<< m_relReduction[0] << "    " << std::setw(0) << fctName(0) << "\n");
		for (size_t i = 1; i < m_fctName.size(); i++)
		{
			print_offset(); UG_LOG("         " << std::scientific << defect(i) << "    " << defect(i)/previousDefect(i)
								<< "    " << reduction(i) << "    "
								<< m_relReduction[i] << "    " << fctName(i) << "\n");
		}
	}
}


template <class TDomain>
bool IndivFctConvCheck<TDomain>::iteration_ended()
{
	bool ended = true;

	for (size_t i = 0; i < m_currentDefect.size(); i++)
	{
		if (!is_valid_number(m_currentDefect[i])) return true;
		if (step() >= m_maxSteps) return true;
		ended = ended && (defect(i) < m_minDefect[i] || reduction(i) < m_relReduction[i]);
	}

	return ended;
}


template <class TDomain>
bool IndivFctConvCheck<TDomain>::post()
{
	bool success = true;

	bool allValid = true;
	std::vector<bool> valid (m_fctName.size(), true);
	std::vector<bool> minDef (m_fctName.size(), true);
	std::vector<bool> red (m_fctName.size(), true);

	for (size_t i = 0; i < m_fctName.size(); i++)
	{
		if (!is_valid_number(m_currentDefect[i]))
		{
			valid[i] = false;
			allValid = false;
		}
		minDef[i] = defect(i) < m_minDefect[i];
		red[i]    = reduction(i) < m_relReduction[i];
	}

	success = allValid && step() < m_maxSteps;

	if (m_verbose)
	{
		if (!success)
		{
			if (!allValid)
				for (size_t i = 0; i < m_fctName.size(); i++)
					if (!valid[i])
					{
						print_offset(); UG_LOG("Current defect for '" << m_fctName[i] << "' is not a valid number.\n");
					}

			if (step() >= m_maxSteps)
			{
				print_offset(); UG_LOG("Maximum numbers of "<< m_maxSteps << " iterations reached without convergence.\n");
			}
		}
		else
		{
			for (size_t i = 0; i < m_fctName.size(); i++)
			{
				if (minDef[i])
				{
					print_offset(); UG_LOG("Absolute defect norm of " << m_minDefect[i] << " for '"
								<< m_fctName[i] << "' reached after " << step() << " steps.\n");
				}

				if (red[i])
				{
					print_offset(); UG_LOG("Relative reduction of " << m_relReduction[i] << " for '"
								<< m_fctName[i] << "' reached after " << step() << " steps.\n");
				}
			}
		}

		print_offset();
		UG_LOG(repeat(m_symbol, 5));
		if (success) {UG_LOG("  Iteration converged  ");}
		else {UG_LOG("  Iteration not successful  ");}
		UG_LOG(repeat(m_symbol, 5));
		UG_LOG("\n\n");
	}

	return success;
}


template <class TDomain>
void IndivFctConvCheck<TDomain>::print_offset()
{
	// step 1: whitespace
	UG_LOG(repeat(' ', m_offset));

	// step 2: print style character
	UG_LOG(m_symbol << " ");
}


template <class TDomain>
bool IndivFctConvCheck<TDomain>::is_valid_number(number value)
{
	if (value == 0.0) return true;
	else return (value >= std::numeric_limits<number>::min()
				&& value <= std::numeric_limits<number>::max()
				&& value == value && value >= 0.0);
}


// explicit class declarations
template class IndivFctConvCheck<Domain1d>;
template class IndivFctConvCheck<Domain2d>;
template class IndivFctConvCheck<Domain3d>;



} // end namespace ug
