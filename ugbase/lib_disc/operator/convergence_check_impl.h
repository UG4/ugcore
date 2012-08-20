/*
 * convergence_check.cpp
 *
 *      Author: M. Breit
 */

#ifndef __H__LIB_DISC__OPERATOR__CONVERGENCE_CHECK_IMPL__
#define __H__LIB_DISC__OPERATOR__CONVERGENCE_CHECK_IMPL__

#include "convergence_check.h"
#include "common/util/string_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Individual function convergence check
////////////////////////////////////////////////////////////////////////////////

template <class TVector, class TDomain>
CompositeConvCheck<TVector, TDomain>::
CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > approx)
 :	 m_initialDefect(0), m_initialOverallDefect(0.0),
  	 m_currentDefect(0), m_currentOverallDefect(0.0),
 	 m_lastDefect(0), m_currentStep(0), m_maxSteps(100),
 	 m_minDefect(0),
 	 m_relReduction(0),
	 m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	 m_timeMeas(true),
	 m_dd(NULL)
{
	m_dd = approx->surface_dof_distribution();

	// compute indices for faster access later
	m_vvMultiIndex.clear();
	m_vvMultiIndex.resize(m_dd->num_fct());
	if(m_dd->max_dofs(VERTEX)) extract_multi_indices<VertexBase>();
	if(m_dd->max_dofs(FACE)) extract_multi_indices<EdgeBase>();
	if(m_dd->max_dofs(EDGE)) extract_multi_indices<Face>();
	if(m_dd->max_dofs(VOLUME)) extract_multi_indices<Volume>();

	// store function name
	for (size_t fi = 0; fi < m_dd->num_fct(); fi++)
		m_fctName.push_back(m_dd->function_pattern().name(fi));
}

template <class TVector, class TDomain>
template <typename TBaseElem>
void CompositeConvCheck<TVector, TDomain>::
extract_multi_indices()
{
	typename SurfaceDoFDistribution::template traits<TBaseElem>::const_iterator iter, iterBegin, iterEnd;

	//	get element iterator for current subset
	iterBegin = m_dd->template begin<TBaseElem>();
	iterEnd = m_dd->template end<TBaseElem>();

	// loop over all elements
	std::vector<MultiIndex<2> > vInd;
	for (iter = iterBegin; iter != iterEnd; ++iter)
	{
		for (size_t fi = 0; fi < m_dd->num_fct(); fi++)
		{
			// inner multi indices of the grid object for a function component
			m_dd->inner_multi_indices(*iter, fi, vInd);

			// remember multi indices
			for (size_t dof = 0; dof < vInd.size(); dof++)
				m_vvMultiIndex[fi].push_back(vInd[dof]);
		}
	}

	// note: no duplicate indices possible
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::set_functions(const char* functionNames)
{
	// get the functions specified by function names
	try {m_fctGrp = m_dd->fct_grp_by_name(functionNames);}
		UG_CATCH_THROW("At least one of the functions in '" << functionNames
						<< "' is not contained in the approximation space (or something else was wrong).");


	// remove unnecessary function names and indices / add 'rest' if necessary
	std::vector<std::vector<MultiIndex<2> > > finalIndices(0);
	std::vector<std::string> finalNames(0);
	std::vector<bool> used(m_dd->num_fct(), false);
	for (size_t i = 0; i < m_fctGrp.num_fct(); i++)
	{
		finalIndices.push_back(m_vvMultiIndex[m_fctGrp[i]]);
		finalNames.push_back(m_fctName[m_fctGrp[i]]);
		used[m_fctGrp[i]] = true;
	}

	std::vector<MultiIndex<2> > rest(0);
	std::stringstream ss; ss << "rest (";
	std::vector<std::string>::iterator itName = m_fctName.begin();
	std::vector<std::vector<MultiIndex<2> > >::iterator itInd = m_vvMultiIndex.begin();
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

	m_vvMultiIndex = finalIndices;
	m_fctName = finalNames;
}


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::set_minimum_defect(const char* minDefect, number minDefectForRest)
{
	//	tokenize strings
	std::vector<std::string> tokens;
	TokenizeString(minDefect, tokens, ',');

	// check if number of values is correct
	if (tokens.size() != m_fctGrp.num_fct())
	{
		UG_THROW(	"The number of supplied values (" << tokens.size() << ") does not match the number\n"
					"of given function names (" << m_fctGrp.num_fct() << "); perhaps you have forgot to call\n"
					"CompositeConvCheck::set_functions prior to this method.");
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


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::set_reduction(const char* reduction, number reductionForRest)
{
	//	tokenize strings
	std::vector<std::string> tokens;
	TokenizeString(reduction, tokens, ',');

	// check if number of values is correct
	if (tokens.size() != m_fctGrp.num_fct())
	{
		UG_THROW(	"The number of supplied values (" << tokens.size() << ") does not match the number\n"
					"of given function names (" << m_fctGrp.num_fct() << "); perhaps you have forgot to call\n"
					"CompositeConvCheck::set_functions prior to this method.");
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


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::start_defect(number initialDefect)
{
	UG_THROW(	"This method cannot be used to set defect values,\n"
				"since obviously this class is meant for an individual\n"
				"defect calculation of more than one function\n"
				"(use start(TVector& d) instead).");
}


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::start(const TVector& vec)
{
	// start time measurement
	if (m_timeMeas)	m_stopwatch.start();

	m_initialDefect.clear();

	// calculate the defect's 2-norm for each function
	for (size_t i = 0; i < m_vvMultiIndex.size(); i++)
	{
		number compDefect = 0.0;
		for(size_t dof = 0; dof < m_vvMultiIndex[i].size(); ++dof)
		{
			const number val = DoFRef(vec, m_vvMultiIndex[i][dof]);
			compDefect += val*val;
		}
		compDefect = sqrt(compDefect);
		m_initialDefect.push_back(compDefect);
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


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::update_defect(number newDefect)
{
	UG_THROW(	"This method cannot be used to update defect values,\n"
				"since obviously this class is meant for an individual\n"
				"defect calculation of more than one function\n"
				"(use update(TVector& d) instead).");
}


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::update(const TVector& vec)
{
	m_currentOverallDefect = 0.0;
	m_lastDefect = m_currentDefect;

	// calculate the defect's 2-norm for each function
	for (size_t i = 0; i < m_vvMultiIndex.size(); i++)
	{
		number compDefect = 0.0;
		for(size_t dof = 0; dof < m_vvMultiIndex[i].size(); ++dof)
		{
			const number val = DoFRef(vec, m_vvMultiIndex[i][dof]);
			compDefect += val*val;
		}
		compDefect = sqrt(compDefect);
		m_currentDefect[i] = compDefect;
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


template <class TVector, class TDomain>
bool CompositeConvCheck<TVector, TDomain>::iteration_ended()
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


template <class TVector, class TDomain>
bool CompositeConvCheck<TVector, TDomain>::post()
{
	if (m_timeMeas) m_stopwatch.stop();

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
		std::stringstream tmsg;
		if (m_timeMeas)
		{
			number time = m_stopwatch.ms()/1000.0;
			tmsg << " (t: " << setprecision(3) << time << "s;  t/it: " << time / step() << "s)";
		}
		if (success) {UG_LOG("  Iteration converged" << tmsg.str() << "  ");}
		else {UG_LOG("  Iteration not successful" << tmsg.str() << "  ");}
		UG_LOG(repeat(m_symbol, 5));
		UG_LOG("\n\n");
	}

	return success;
}


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::print_offset()
{
	// step 1: whitespace
	UG_LOG(repeat(' ', m_offset));

	// step 2: print style character
	UG_LOG(m_symbol << " ");
}


template <class TVector, class TDomain>
bool CompositeConvCheck<TVector, TDomain>::is_valid_number(number value)
{
	if (value == 0.0) return true;
	else return (value >= std::numeric_limits<number>::min()
				&& value <= std::numeric_limits<number>::max()
				&& value == value && value >= 0.0);
}

} // end namespace ug


#endif /* __H__LIB_DISC__OPERATOR__CONVERGENCE_CHECK_IMPL__ */
