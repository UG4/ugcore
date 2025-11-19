/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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
 *  algebraic version of composite_conv_check.h by M. Breit
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_CONV_CHECK_IMPL__
#define __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_CONV_CHECK_IMPL__

#include "algebra_conv_check.h"
#include "common/util/string_util.h"
#include "lib_algebra/adapter/scalar_subvector_adapter.hh"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Algebraic (composite) convergence check
////////////////////////////////////////////////////////////////////////////////

template <class TVector>
AlgebraicConvCheck<TVector>::AlgebraicConvCheck(size_t ncmp)
:	
	m_vCmpInfo(ncmp, CmpInfo(1e-12, 1e-10)),
	m_maxSteps(100), m_minDefect(1e-12), m_relReduction(1e-10), m_verbose(true),
 	m_currentStep(0), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	m_bTimeMeas(true)
{}


template <class TVector>
AlgebraicConvCheck<TVector>::
AlgebraicConvCheck(size_t ncmp, int maxSteps, number minDefect, number relReduction)
:	 m_maxSteps(maxSteps), m_minDefect(minDefect), m_relReduction(relReduction), m_verbose(true),
 	 m_currentStep(0), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	m_bTimeMeas(true), m_vCmpInfo(ncmp, CmpInfo(minDefect, relReduction))
{
	//const CmpInfo(minDefect, relReduction);

}

template <class TVector>
AlgebraicConvCheck<TVector>::
AlgebraicConvCheck(size_t ncmp, int maxSteps, number minDefect, number relReduction, bool verbose)
:	 m_maxSteps(maxSteps), m_minDefect(minDefect), m_relReduction(relReduction), m_verbose(verbose),
 	 m_currentStep(0), m_offset(0), m_symbol('%'), m_name("Iteration"), m_info(""),
	m_bTimeMeas(true)
{
	//const CmpInfo(minDefect, relReduction);
	//m_vCmpInfo(ncmp, a);
}





template <class TVector>
SmartPtr<IConvergenceCheck<TVector> > AlgebraicConvCheck<TVector>::clone()
{
	SmartPtr<AlgebraicConvCheck<TVector> > newInst(new AlgebraicConvCheck<TVector>(this->m_vCmpInfo.size()));

	// use std assignment (implicit member-wise is fine here)
	*newInst = *this;
	return newInst;
}



template <class TVector>
void AlgebraicConvCheck<TVector>::start_defect(number initialDefect)
{
	UG_THROW(	"This method cannot be used to set defect values,\n"
				"since obviously this class is meant for an individual\n"
				"defect calculation of more than one function\n"
				"(use start(TVector& d) instead).");
}


template <class TVector>
void AlgebraicConvCheck<TVector>::start(const TVector& vec)
{
	// start time measurement
	if (m_bTimeMeas)	m_stopwatch.start();

	// update defects
	for (size_t fct = 0; fct < m_vCmpInfo.size(); fct++){
		m_vCmpInfo[fct].currDefect = m_vCmpInfo[fct].initDefect = norm(vec, fct);
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
			UG_LOG(repeat(m_symbol, space_left))
			UG_LOG("\n");
		}
		else
		{
			UG_LOG("\n");
		}

	//	start iteration output
		print_offset();
		UG_LOG("  Iter      Defect        Rate        \n");

		for (size_t cmp = 0; cmp < m_vCmpInfo.size(); cmp++)
		{
			CmpInfo& cmpInfo = m_vCmpInfo[cmp];

			print_offset();
			if(cmp != 0) {UG_LOG("         " )}
			else {UG_LOG(std::setw(4) << step() << ":    ");}


			UG_LOG(std::scientific << cmpInfo.currDefect <<  "    ")
			UG_LOG(std::scientific << "  -----  " << "    c=" << cmp << std::endl)

			//UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.minDefect <<   "    " << std::setprecision(6) )

			//UG_LOG(std::scientific << "  --------  " << "    ")
			//UG_LOG(std::scientific << std::setprecision(3) << cmpInfo.relReduction << "    " << std::setprecision(6))
			//UG_LOG(std::scientific << cmpInfo.name << "\n");
		}
	}
}


template <class TVector>
void AlgebraicConvCheck<TVector>::update_defect(number newDefect)
{
	UG_THROW(	"This method cannot be used to update defect values,\n"
				"since obviously this class is meant for an individual\n"
				"defect calculation of more than one function\n"
				"(use update(TVector& d) instead).");
}


template <class TVector>
void AlgebraicConvCheck<TVector>::update(const TVector& vec)
{
	// update defects
	for (size_t fct = 0; fct < m_vCmpInfo.size(); fct++){
		m_vCmpInfo[fct].lastDefect = m_vCmpInfo[fct].currDefect;
		m_vCmpInfo[fct].currDefect = norm(vec, fct);
	}

	//vec.print(); // debug only!

	m_currentStep++;

	if (m_verbose)
	{
		for (size_t cmp = 0; cmp < m_vCmpInfo.size(); cmp++)
		{
			CmpInfo& cmpInfo = m_vCmpInfo[cmp];

			print_offset();
			if(cmp != 0) {UG_LOG("         " )}
			else {UG_LOG(std::setw(4) << step() << ":    ");}

			UG_LOG(std::scientific << cmpInfo.currDefect <<  "    ")
			if(cmpInfo.lastDefect != 0.0){
				UG_LOG(std::scientific << std::setprecision(6) << cmpInfo.currDefect / cmpInfo.lastDefect << "    c="<< cmp<< std::setprecision(6)  << std::endl)
			} else {
				UG_LOG(std::scientific << "  -----  " << "    c="<< cmp  << std::endl)
			}
		}
	}
}


template <class TVector>
bool AlgebraicConvCheck<TVector>::iteration_ended()
{
	if (step() >= m_maxSteps) return true;

	bool ended = true;
	for (size_t cmp = 0; cmp < m_vCmpInfo.size(); cmp++)
	{
		CmpInfo& cmpInfo = m_vCmpInfo[cmp];

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


template <class TVector>
bool AlgebraicConvCheck<TVector>::post()
{
	if (m_bTimeMeas) m_stopwatch.stop();

	bool success = true;

	if (step() >= m_maxSteps){
		print_offset();
		UG_LOG("Maximum numbers of "<< m_maxSteps <<
		       	  " iterations reached without convergence.\n");
		success = false;
	}

	for (size_t cmp = 0; cmp < m_vCmpInfo.size(); cmp++)
	{
		CmpInfo& cmpInfo = m_vCmpInfo[cmp];

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


template <class TVector>
void AlgebraicConvCheck<TVector>::print_offset()
{
	// step 1: whitespace
	UG_LOG(repeat(' ', m_offset));

	// step 2: print style character
	UG_LOG(m_symbol << " ");
}

template <class TVector>
void AlgebraicConvCheck<TVector>::print_line(std::string line)
{
	print_offset();
	UG_LOG(line << "\n");
}


template <class TVector>
bool AlgebraicConvCheck<TVector>::is_valid_number(number value)
{
	if (value == 0.0) return true;
	else return (value >= std::numeric_limits<number>::min()
				&& value <= std::numeric_limits<number>::max()
				&& value == value && value >= 0.0);
}





template <class TVector>
number AlgebraicConvCheck<TVector>::
norm(const TVector& vec, size_t cmp)
{
#ifdef UG_PARALLEL

	// 	make vector d additive unique
	if (!const_cast<TVector*>(&vec)->change_storage_type(PST_UNIQUE))
		UG_THROW("AlgebraicConvCheck::norm(): Cannot change ParallelStorageType to unique.");
#endif


	ConstScalarSubVectorAdapter<TVector, typename ug::Vector<double> > dummy (vec,cmp);
	double norm = VecNormSquared(dummy);
/*	for (size_t dof = 0; dof < vMultiIndex.size(); ++dof)
	{
		const number val = BlockRef(vec, vMultiIndex[dof]);
		norm += (double) (val*val);
	}
*/
#ifdef UG_PARALLEL
	// sum squared local norms
	norm = vec.layouts()->proc_comm().allreduce(norm, PCL_RO_SUM);
#endif

	// return global norm
	return sqrt((number) norm);
}

} // end namespace ug


#endif