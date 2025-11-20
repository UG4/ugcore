/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_mg_stats_impl
#define __H__UG_mg_stats_impl

#include "common/util/stringify.h"
#include "lib_disc/function_spaces/grid_function_util.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
MGStats<TDomain, TAlgebra>::
MGStats() :
	m_statsRow(-1),
	m_lastLvlWritten(-1),
	m_maxLvl(0),
	m_filenamePrefix("mgstats"),
	m_exitOnError(false),
	m_writeErrVecs(false),
	m_writeErrDiffs(false)
{
	for(int i = 0; i < NUM_STAGES; ++i)
		m_stageIsActive[i] = true;
	m_stageIsActive[INVALID] = false;
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
set_exit_on_error(bool exitOnError)
{
	m_exitOnError = exitOnError;
}


template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
set_write_err_vecs(bool writeErrVecs)
{
	m_writeErrVecs = writeErrVecs;
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
set_write_err_diffs(bool writeErrDiffs)
{
	m_writeErrDiffs = writeErrDiffs;
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
set_filename_prefix(const char* filename)
{
	m_filenamePrefix = filename;
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
set_active_stages(const std::vector<int>& activeStages)
{
	for(int i = 0; i < NUM_STAGES; ++i)
		m_stageIsActive[i] = false;

	for(size_t i = 0; i < activeStages.size(); ++i){
		if(activeStages[i] >= 0 && (int)activeStages[i] < NUM_STAGES)
			m_stageIsActive[activeStages[i]] = true;
	}

	m_stageIsActive[INVALID] = false;
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
save_stats_to_file()
{
	std::string filename = mkstr(m_filenamePrefix << ".log");
	save_stats_to_file(filename.c_str());
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
save_stats_to_file(const char* filename)
{
	#ifdef UG_PARALLEL
	if(pcl::ProcRank() == 0){
	#endif

		std::ofstream out(filename);
		UG_COND_THROW(!out, "MGStats: Couldn't open '" << filename << "' for writing.");
		out << m_stats << std::endl;
		out.close();

	#ifdef UG_PARALLEL
	}
	#endif
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
print()
{
	UG_LOG(m_stats << std::endl);
}

template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
clear()
{
	m_stats.clear();
}


template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
set_defect(grid_func_t& gf, int lvl, Stage stage)
{
	if(lvl < 0)
		return;

	if(!m_stageIsActive[stage])
		return;

	level_required(lvl);
	FuncEntry& fe = m_funcs[lvl];
	
	if(stage == INVALID){
		fe.stage = INVALID;
		return;
	}

	if(fe.func.invalid() || fe.func->num_indices() != gf.num_indices()){
		fe.func = gf.clone_without_values();
		fe.tmpFunc = gf.clone_without_values();
		fe.stage = INVALID;
		fe.norm = 0;
	}

//	this is just an assumption that a new cycle started
	if(stage <= fe.stage && lvl >= m_lastLvlWritten){
		++m_statsRow;
		m_lastLvlWritten = -1;
		write_header(lvl);
	//	write a line to cout so that the output and the stats can be connected
		UG_LOG("    >>> mg-stats row: " << m_statsRow << std::endl);
	}

	++m_statsRow;

	const int r = m_statsRow;

	m_lastLvlWritten = lvl;

//	write the current row index
	m_stats(r, 0) = mkstr(r);
	m_stats(r, 1) = stage_norm_name(stage);
	
	
//	compute statistics
//	copy gf to fe.tmpFunc to avoid any side effects on the multigrid iteration
	fe.tmpFunc->assign(gf);
	grid_func_t& newFunc = *fe.tmpFunc;
	const number norm = newFunc.norm();
	m_stats(r, 2 + m_maxLvl - lvl) = mkstr(norm);

	if((fe.stage < stage) && (norm > fe.norm)){
	//	during the last operation the defect grew.
	//	output the old and new grid functions and their difference
		std::string filename;

		if(m_writeErrVecs){
			filename = mkstr(m_filenamePrefix << "-row-" << r << "-lvl-"
	                            << lvl << "--" << fe.stage << "-"
	                            << stage_name(fe.stage) << ".vec");

			SaveVectorForConnectionViewer((*fe.func), filename.c_str());


			filename = mkstr(m_filenamePrefix << "-row-" << r << "-lvl-"
			                    << lvl << "--" << stage << "-"
			                    << stage_name(stage) << ".vec");

			SaveVectorForConnectionViewer(gf, filename.c_str());
		}

		if(m_writeErrDiffs){
			filename = mkstr(m_filenamePrefix << "-row-" << r << "-lvl-"
			                    << lvl << "-diff--" << stage << "-"
			                    << stage_name(fe.stage) << "-"
			                    << stage_name(stage) << ".vec");
			SaveVectorDiffForConnectionViewer(gf, *fe.func, filename.c_str());
		}

		std::string errMsg =
			mkstr("MGStats: Defect deteriorated on level " << lvl
			      << " between stages '" << stage_name(fe.stage) << "' and '"
			      << stage_name(stage) << "'. (MGStats row: " << r << ")");

		if(m_exitOnError){
			save_stats_to_file();
			print();
			UG_THROW(errMsg);
		}
		else{
			UG_LOG(errMsg << std::endl);
		}
	}

//	swap grid function of fe
	sp_grid_func_t sptmp = fe.func;
	fe.func = fe.tmpFunc;
	fe.tmpFunc = sptmp;
	fe.norm = norm;
	fe.stage = stage;
}


template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
level_required(int lvl)
{
	if((int)m_funcs.size() <= lvl)
		m_funcs.resize(lvl + 1);
}


template <typename TDomain, typename TAlgebra>
void MGStats<TDomain, TAlgebra>::
write_header(int maxLvl)
{
	m_maxLvl = maxLvl;
	m_stats(m_statsRow, 0) = "row";
	m_stats(m_statsRow, 1) = "stage";
	for(int i = 0; i < maxLvl; ++i){
		m_stats(m_statsRow, i+2) = mkstr("lvl " << maxLvl - i);
	}
}


}//	end of namespace

#endif