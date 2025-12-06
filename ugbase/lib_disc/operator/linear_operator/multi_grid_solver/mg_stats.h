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

#ifndef __H__UG_mg_stats
#define __H__UG_mg_stats

#include <string>
#include <vector>
#include "common/util/table.h"

namespace ug{

///	Records statistics on how individual parts of a multigrid method worked
template <typename TDomain, typename TAlgebra>
class MGStats {
	public:
		using grid_func_t = GridFunction<TDomain, TAlgebra>;
		using sp_grid_func_t = SmartPtr<grid_func_t>;


	///	Defines at which stage data is recorded in a given multigrid cycle
		enum Stage {
			BEFORE_PRE_SMOOTH,
			AFTER_PRE_SMOOTH,
			BEFORE_POST_SMOOTH,
			AFTER_POST_SMOOTH,
			INVALID					// always last!
		};

		static constexpr int NUM_STAGES = INVALID + 1;

		MGStats();

	///	If enabled, a deterioration of the norm of the defect leads to an error
	/**	disabled by default.*/
		void set_exit_on_error(bool exitOnError);

	///	If enabled, involved defects are written to file if the defect deteriorates
	/**	disabled by default.*/
		void set_write_err_vecs(bool writeErrVecs);

	///	If enabled, a diff bettween defects involved is written to file if the defect deteriorates
	/**	disabled by default.*/
		void set_write_err_diffs(bool writeErrDiffs);

	///	sets the active stages. All other stages will be ignored.
	/**	\param activeStages		vector containing constants from MGStats::Stage.
	 *							Each listed stage will be considered active, all
	 *							non-listed stages will be ignored.
	 */
		void set_active_stages(const std::vector<int>& activeStages);

	///	sets the prefix with which files are written
	/**	This concerns the log-file and grid-function files.
	 * Default is 'mgstats'.*/
		void set_filename_prefix(const char* filename);

	///	saves current stats to filenamePrefix.log
		void save_stats_to_file() const;

	///	saves current stats to the specified file
		void save_stats_to_file(const char* filename) const;

	///	prints the current stats
		void print() const;

	///	clears the current stats
		void clear();

	///	set the defect on a certain level for a given stage
	/**	If the defect for the previous stage was set for the same level,
	 * norms are compared and a diff can be computed.
	 *
	 * If 'exitOnError' is enabled, the method will print the current stats and
	 * throw an instance of UGError of the norm of the defect deteriorated.
	 *
	 * \note	In parallel environments this method is a synchronization point.
	 *			Make sure to call it from all processes which hold 'gf'.*/
		void set_defect(grid_func_t& gf, int lvl, Stage stage);

	///	returns the name of a given stage as a string
		const char* stage_name(Stage stage) {
			switch (stage) {
				case BEFORE_PRE_SMOOTH:  return "BEFORE_PRE_SMOOTH";
				case AFTER_PRE_SMOOTH:   return "AFTER_PRE_SMOOTH";
				case BEFORE_POST_SMOOTH: return "BEFORE_POST_SMOOTH";
				case AFTER_POST_SMOOTH:  return "AFTER_POST_SMOOTH";
				default:                 return "INVALID";
			}
		}

	///	returns the name of the norm of a given stage as a string
		const char* stage_norm_name(Stage stage) {
			switch (stage) {
				case BEFORE_PRE_SMOOTH:  return "|bef pre smth|";
				case AFTER_PRE_SMOOTH:   return "|aft pre smth|";
				case BEFORE_POST_SMOOTH: return "|bef post smth|";
				case AFTER_POST_SMOOTH:  return "|aft post smth|";
				default:                 return "|INVALID|";
			}
		}

	private:
		void level_required(int lvl);
		void write_header(int maxLvl);

		struct FuncEntry{
			FuncEntry() : stage(INVALID), norm(0) {}
			sp_grid_func_t	func;
			sp_grid_func_t	tmpFunc;
			Stage			stage;
			number			norm;
		};

		std::vector<FuncEntry>	m_funcs;
		bool		m_stageIsActive[NUM_STAGES];
		StringTable	m_stats;
		int			m_statsRow;
		int			m_lastLvlWritten;
		int			m_maxLvl;
		std::string	m_filenamePrefix;
		bool		m_exitOnError;
		bool		m_writeErrVecs;
		bool		m_writeErrDiffs;
};


}//	end of namespace

////////////////////////////////
//	include implementation
#include "mg_stats_impl.hpp"


#endif