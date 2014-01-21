/*
 * composite_conv_check.cpp
 *
 *      Author: M. Breit
 */

#ifndef __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK_IMPL__
#define __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK_IMPL__

#include "composite_conv_check.h"
#include "common/util/string_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Composite convergence check
////////////////////////////////////////////////////////////////////////////////

template <class TVector, class TDomain>
CompositeConvCheck<TVector, TDomain>::
CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
 :	m_spApprox(spApproxSpace),
  	m_bCheckRest(true), m_restMinDefect(1-12), m_restRelReduction(1-10),
	m_currentStep(0), m_maxSteps(100),
	m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	m_bTimeMeas(true)
{
	set_level(GridLevel::TOP);
}

template <class TVector, class TDomain>
CompositeConvCheck<TVector, TDomain>::
CompositeConvCheck(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                   int maxSteps, number minDefect, number relReduction)
 :	m_spApprox(spApproxSpace),
  	m_bCheckRest(true), m_restMinDefect(minDefect), m_restRelReduction(relReduction),
	m_currentStep(0), m_maxSteps(maxSteps),
	m_verbose(true), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	m_bTimeMeas(true)
{
	set_level(GridLevel::TOP);
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::set_level(int level)
{
	ConstSmartPtr<DoFDistribution> dd = m_spApprox->dof_distribution(GridLevel(level, GridLevel::SURFACE));
	extract_dof_indices(dd);

	update_rest_check();
}


template <class TVector, class TDomain>
template <typename TBaseElem>
void CompositeConvCheck<TVector, TDomain>::
extract_dof_indices(ConstSmartPtr<DoFDistribution> dd)
{
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterBegin, iterEnd;

	//	get element iterator for current subset
	iterBegin = dd->template begin<TBaseElem>();
	iterEnd = dd->template end<TBaseElem>();

	// loop over all elements
	std::vector<DoFIndex> vInd;
	for (iter = iterBegin; iter != iterEnd; ++iter)
	{
		for (size_t fi = 0; fi < dd->num_fct(); fi++)
		{
			// inner multi indices of the grid object for a function component
			dd->inner_dof_indices(*iter, fi, vInd);

			// remember multi indices
			for (size_t dof = 0; dof < vInd.size(); dof++)
				m_vNativCmpInfo[fi].vMultiIndex.push_back(vInd[dof]);
		}
	}

	// note: no duplicate indices possible
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
extract_dof_indices(ConstSmartPtr<DoFDistribution> dd)
{
	// compute indices for faster access later
	m_vNativCmpInfo.clear();
	m_vNativCmpInfo.resize(dd->num_fct());
	if(dd->max_dofs(VERTEX)) extract_dof_indices<VertexBase>(dd);
	if(dd->max_dofs(FACE)) extract_dof_indices<EdgeBase>(dd);
	if(dd->max_dofs(EDGE)) extract_dof_indices<Face>(dd);
	if(dd->max_dofs(VOLUME)) extract_dof_indices<Volume>(dd);

	for(size_t i = 0; i < m_vNativCmpInfo.size(); ++i)
		m_vNativCmpInfo[i].name = dd->name(i);

	m_numAllDoFs = 0;
	for(size_t i = 0; i < m_vNativCmpInfo.size(); ++i)
		m_numAllDoFs += m_vNativCmpInfo[i].vMultiIndex.size();
}


template <class TVector, class TDomain>
number CompositeConvCheck<TVector, TDomain>::
norm(const TVector& vec, const std::vector<DoFIndex>& vMultiIndex)
{
#ifdef UG_PARALLEL

	// 	make vector d additive unique
	if (!const_cast<TVector*>(&vec)->change_storage_type(PST_UNIQUE))
		UG_THROW("CompositeConvCheck::norm(): Cannot change ParallelStorageType to unique.");
#endif

	double norm = 0.0;
	for (size_t dof = 0; dof < vMultiIndex.size(); ++dof)
	{
		const number val = DoFRef(vec, vMultiIndex[dof]);
		norm += (double) (val*val);
	}

#ifdef UG_PARALLEL
	// sum squared local norms
	norm = vec.layouts()->proc_comm().allreduce(norm, PCL_RO_SUM);
#endif

	// return global norm
	return sqrt((number) norm);
}


template <class TVector, class TDomain>
SmartPtr<IConvergenceCheck<TVector> > CompositeConvCheck<TVector, TDomain>::clone()
{
	SmartPtr<CompositeConvCheck<TVector, TDomain> > newInst(
						new CompositeConvCheck<TVector, TDomain>(m_spApprox));

	// use std assignment (implicit member-wise is fine here)
	*newInst = *this;
	return newInst;
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_component_check(const std::vector<std::string>& vFctName,
                    const std::vector<number>& vMinDefect,
                    const std::vector<number>& vRelReduction)
{
	if(vFctName.size() != vMinDefect.size() || vFctName.size() != vRelReduction.size())
		UG_THROW("CompositeConvCheck: Please specify one value for each cmp.");

	for(size_t cmp = 0; cmp < vFctName.size(); ++cmp)
		set_component_check(vFctName[cmp], vMinDefect[cmp], vRelReduction[cmp]);

}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_component_check(const std::string& vFctName,
                    const std::vector<number>& vMinDefect,
                    const std::vector<number>& vRelReduction)
{
	set_component_check(TokenizeTrimString(vFctName), vMinDefect, vRelReduction);
}


template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_component_check(const std::vector<std::string>& vFctName,
					 const number minDefect,
					 const number relReduction)
{
	for(size_t cmp = 0; cmp < vFctName.size(); ++cmp)
		set_component_check(vFctName[cmp], minDefect, relReduction);
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_all_component_check(const number minDefect,
                        const number relReduction)
{
	for(size_t fct = 0; fct < m_spApprox->num_fct(); ++fct){
		std::string name = m_spApprox->name(fct);
		set_component_check(name, minDefect, relReduction);
	}
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_component_check(const std::string& fctNames,
                    const number minDefect,
                    const number relReduction)
{
	std::vector<std::string> vName = TokenizeTrimString(fctNames);

	for(size_t i = 0; i < vName.size(); ++i)
	{
		const int fct = m_spApprox->fct_id_by_name(vName[i].c_str());

		// search for single component
		size_t cmp = 0;
		for(; cmp < m_CmpInfo.size(); ++cmp){

			// not checking rest
			if(m_CmpInfo[cmp].isRest) continue;

			// only single valued checked
			if(m_CmpInfo[cmp].vFct.size() != 1) continue;

			// check if component found
			if(m_CmpInfo[cmp].vFct[0] == fct)
				break;
		}

		// if not found add new one
		if(cmp == m_CmpInfo.size())
			m_CmpInfo.push_back(CmpInfo());

		// set values
		m_CmpInfo[cmp].isRest = false;
		m_CmpInfo[cmp].vFct.clear();
		m_CmpInfo[cmp].vFct.push_back(fct);
		m_CmpInfo[cmp].name = vName[i];
		m_CmpInfo[cmp].minDefect = minDefect;
		m_CmpInfo[cmp].relReduction = relReduction;
	}

	// update rest values
	update_rest_check();
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_group_check(const std::vector<std::string>& vFctName,
                const number minDefect,
                const number relReduction)
{
	std::vector<int> vFct(vFctName.size());
	for(size_t i = 0; i < vFctName.size(); ++i){
		std::string name = TrimString(vFctName[i]);
		vFct[i] = m_spApprox->fct_id_by_name(name.c_str());
	}

	// search for group
	size_t cmp = 0;
	for(; cmp < m_CmpInfo.size(); ++cmp){

		// not checking rest
		if(m_CmpInfo[cmp].isRest) continue;

		// only single valued checked
		if(m_CmpInfo[cmp].vFct.size() != vFct.size()) continue;

		// check if component found
		bool bFound = true;
		for(size_t i = 0; i < vFct.size(); ++i){
			if(m_CmpInfo[cmp].vFct[i] != vFct[i])
				bFound = false;
		}

		if(bFound)	break;
	}

	// if not found add new one
	if(cmp == m_CmpInfo.size())
		m_CmpInfo.push_back(CmpInfo());

	std::string name("");
	for(size_t fct = 0; fct < vFct.size(); ++fct){
		if(!name.empty()) name.append(", ");
		name.append(m_vNativCmpInfo[vFct[fct]].name);
	}

	// set values
	m_CmpInfo[cmp].isRest = false;
	m_CmpInfo[cmp].vFct = vFct;
	m_CmpInfo[cmp].name = name;
	m_CmpInfo[cmp].minDefect = minDefect;
	m_CmpInfo[cmp].relReduction = relReduction;

	// update rest values
	update_rest_check();
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
set_group_check(const std::string& fctNames,
                const number minDefect,
                const number relReduction)
{
	set_group_check(TokenizeTrimString(fctNames), minDefect, relReduction);
}

template <class TVector, class TDomain>
void CompositeConvCheck<TVector, TDomain>::
update_rest_check()
{
	// remove old rest
	for(size_t cmp = 0; cmp < m_CmpInfo.size(); ++cmp)
		if(m_CmpInfo[cmp].isRest)
			m_CmpInfo.erase(m_CmpInfo.begin()+cmp);

	// check if rest required
	if(!m_bCheckRest) return;

	// detect used functions
	std::vector<bool> vUsed(m_vNativCmpInfo.size(), false);
	for(size_t cmp = 0; cmp < m_CmpInfo.size(); ++cmp)
	{
		for(size_t i = 0; i < m_CmpInfo[cmp].vFct.size(); ++i)
			vUsed[ m_CmpInfo[cmp].vFct[i] ] = true;
	}

	std::vector<int> vFct;
	std::string name("");

	for(size_t fct = 0; fct < vUsed.size(); ++fct){
		if(vUsed[fct]) continue;

		vFct.push_back(fct);
		if(!name.empty()) name.append(", ");
		name.append(m_vNativCmpInfo[fct].name);
	}

	// if no rest --> not added
	if(vFct.empty()) return;

	// add new rest
	m_CmpInfo.push_back(CmpInfo());
	m_CmpInfo.back().isRest = true;
	m_CmpInfo.back().vFct = vFct;
	m_CmpInfo.back().name = name;
	m_CmpInfo.back().minDefect = m_restMinDefect;
	m_CmpInfo.back().relReduction = m_restRelReduction;
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
	// assert correct number of dofs
	if (vec.size() != m_numAllDoFs)
	{
		UG_THROW("Number of dofs in CompositeConvCheck does not match"
				"number of dofs given in vector from algorithm (" << m_numAllDoFs
				<< ", but " << vec.size() << " given). \nMake sure that you set "
						"the right grid level via set_level().");
	}

	// start time measurement
	if (m_bTimeMeas)	m_stopwatch.start();

	// update native defects
	for (size_t fct = 0; fct < m_vNativCmpInfo.size(); fct++){
		m_vNativCmpInfo[fct].initDefect = norm(vec, m_vNativCmpInfo[fct].vMultiIndex);
		m_vNativCmpInfo[fct].currDefect = m_vNativCmpInfo[fct].initDefect;
	}

	// update grouped defects
	for (size_t cmp = 0; cmp < m_CmpInfo.size(); cmp++){
		CmpInfo& cmpInfo = m_CmpInfo[cmp];

		cmpInfo.currDefect = 0.0;
		for(size_t i = 0; i < cmpInfo.vFct.size(); ++i)
			cmpInfo.currDefect += pow(m_vNativCmpInfo[cmpInfo.vFct[i]].currDefect, 2);
		cmpInfo.currDefect = sqrt(cmpInfo.currDefect);

		cmpInfo.initDefect = cmpInfo.currDefect;
	}

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
			UG_LOG(repeat(m_symbol, space_left));
			UG_LOG("\n");
		}
		else
		{
			UG_LOG("\n");
		}

	//	start iteration output
		print_offset();
		UG_LOG("  Iter      Defect        Required       Rate        "
				"Reduction     Required     Component(s)\n");

		for (size_t cmp = 0; cmp < m_CmpInfo.size(); cmp++)
		{
			CmpInfo& cmpInfo = m_CmpInfo[cmp];

			print_offset();
			if(cmp != 0) {UG_LOG("         " );}
			else {UG_LOG(std::setw(4) << step() << ":    ");}

			UG_LOG(std::scientific << cmpInfo.currDefect <<  "    ");
			UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.minDefect <<   "    " << std::setprecision(6) );
			UG_LOG(std::scientific << "  -----  " << "    ");
			UG_LOG(std::scientific << "  --------  " << "    ");
			UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.relReduction << "    " << std::setprecision(6));
			UG_LOG(std::scientific << cmpInfo.name << "\n");
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
	// assert correct number of dofs
	if (vec.size() != m_numAllDoFs)
	{
		UG_THROW("Number of dofs in CompositeConvCheck does not match"
				 "number of dofs given in vector from algorithm (" << m_numAllDoFs <<
				 ", but " << vec.size() << " given). \nMake sure that you set the "
				"right grid level via set_level().");
	}

	// update native defects
	for (size_t fct = 0; fct < m_vNativCmpInfo.size(); fct++){
		m_vNativCmpInfo[fct].lastDefect = m_vNativCmpInfo[fct].currDefect;
		m_vNativCmpInfo[fct].currDefect = norm(vec, m_vNativCmpInfo[fct].vMultiIndex);
	}

	// update grouped defects
	for (size_t cmp = 0; cmp < m_CmpInfo.size(); cmp++){
		CmpInfo& cmpInfo = m_CmpInfo[cmp];

		cmpInfo.lastDefect = cmpInfo.currDefect;

		cmpInfo.currDefect = 0.0;
		for(size_t i = 0; i < cmpInfo.vFct.size(); ++i)
			cmpInfo.currDefect += pow(m_vNativCmpInfo[cmpInfo.vFct[i]].currDefect, 2);
		cmpInfo.currDefect = sqrt(cmpInfo.currDefect);
	}

	m_currentStep++;

	if (m_verbose)
	{
		for (size_t cmp = 0; cmp < m_CmpInfo.size(); cmp++)
		{
			CmpInfo& cmpInfo = m_CmpInfo[cmp];

			print_offset();
			if(cmp != 0) {UG_LOG("         " );}
			else {UG_LOG(std::setw(4) << step() << ":    ");}

			UG_LOG(std::scientific << cmpInfo.currDefect <<  "    ");
			UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.minDefect <<   "    " << std::setprecision(6));
			if(cmpInfo.lastDefect != 0.0){
				UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.currDefect / cmpInfo.lastDefect << "    "<< std::setprecision(6));
			} else {
				UG_LOG(std::scientific << "  -----  " << "    ");
			}
			if(cmpInfo.initDefect != 0.0){
				UG_LOG(std::scientific << cmpInfo.currDefect / cmpInfo.initDefect << "    ");
			} else {
				UG_LOG(std::scientific << "  --------  " << "    ");
			}
			UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.relReduction << "    " << std::setprecision(6));
			UG_LOG(std::scientific << cmpInfo.name << "\n");
		}
	}
}


template <class TVector, class TDomain>
bool CompositeConvCheck<TVector, TDomain>::iteration_ended()
{
	if (step() >= m_maxSteps) return true;

	bool ended = true;
	for (size_t cmp = 0; cmp < m_CmpInfo.size(); cmp++)
	{
		CmpInfo& cmpInfo = m_CmpInfo[cmp];

		if (!is_valid_number(cmpInfo.currDefect)) return true;

		bool cmpFinished = false;
		if(cmpInfo.currDefect < cmpInfo.minDefect) cmpFinished = true;
		if(cmpInfo.initDefect != 0.0)
			if((cmpInfo.currDefect/cmpInfo.initDefect) < cmpInfo.relReduction)
				cmpFinished = true;

		ended = ended && cmpFinished;
	}

	return ended;
}


template <class TVector, class TDomain>
bool CompositeConvCheck<TVector, TDomain>::post()
{
	if (m_bTimeMeas) m_stopwatch.stop();

	bool success = true;

	if (step() >= m_maxSteps){
		print_offset();
		UG_LOG("Maximum numbers of "<< m_maxSteps <<
		       	  " iterations reached without convergence.\n");
		success = false;
	}

	for (size_t cmp = 0; cmp < m_CmpInfo.size(); cmp++)
	{
		CmpInfo& cmpInfo = m_CmpInfo[cmp];

		if (!is_valid_number(cmpInfo.currDefect)){
			success = false;
			if(m_verbose){
				print_offset();
				UG_LOG("Current defect for '" << cmpInfo.name <<
					   "' is not a valid number.\n");
			}
		}

		bool cmpFinished = false;
		if (cmpInfo.currDefect < cmpInfo.minDefect)
		{
			print_offset();
			UG_LOG("Absolute defect    of " << cmpInfo.minDefect << " for '"
			       << cmpInfo.name << "' reached after " << step() << " steps.\n");
			cmpFinished = true;
		}

		if(cmpInfo.initDefect != 0.0){
			if (cmpInfo.currDefect/cmpInfo.initDefect < cmpInfo.relReduction)
			{
				print_offset();
				UG_LOG("Relative reduction of " << cmpInfo.relReduction << " for '"
					   << cmpInfo.name << "' reached after " << step() << " steps.\n");
				cmpFinished = true;
			}
		}

		if(!cmpFinished) success = false;
	}

	if(m_verbose){
		print_offset();
		UG_LOG(repeat(m_symbol, 5));
		std::stringstream tmsg;
		if (m_bTimeMeas)
		{
			number time = m_stopwatch.ms()/1000.0;
			tmsg << " (t: " << std::setprecision(3) << time << "s;  t/it: "
					<< time / step() << "s)";
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
void CompositeConvCheck<TVector, TDomain>::print_line(std::string line)
{
	print_offset();
	UG_LOG(line << "\n");
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


#endif /* __H__LIB_DISC__OPERATOR__COMPOSITE_CONVERGENCE_CHECK_IMPL__ */

