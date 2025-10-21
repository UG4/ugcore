/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_COUNT__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_COUNT__

#include "lib_grid/tools/surface_view.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "dof_distribution_info.h"

namespace ug{

class DoFCount : public DoFDistributionInfoProvider
{
	public:
		static const unsigned char ES_MAX = (ES_V_MASTER | ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE);
		static const unsigned char SS_MAX = SurfaceView::ALL;

		static const int ALL_FCT = -1;
		static const int ALL_SUBSET = -1;
		static const unsigned char ALL_ES = ES_MAX + 1;
		static const unsigned char ALL_SS = SS_MAX + 1;
		static const unsigned char UNIQUE_ES = ES_MAX + 2;
		static const unsigned char UNIQUE_SS = SS_MAX + 2;

	protected:
		struct Cnt{

			Cnt();
			void add(uint64 num, SurfaceView::SurfaceState ss, unsigned char is);
			void collect_values(std::vector<uint64>& vNum) const;
			void set_values(const std::vector<uint64>& vNum, size_t& cnt);

			struct PCnt {

				PCnt();

				// dofs exactly matching this state
				uint64 num(unsigned char is) const;

				// dofs that contain this state
				uint64 num_contains(unsigned char is) const;

				void collect_values(std::vector<uint64>& vNum) const;
				void set_values(const std::vector<uint64>& vNum, size_t& cnt);

				std::vector<uint64> vNumIS;
			};

			uint64 num(SurfaceView::SurfaceState ss, unsigned char is) const;
			uint64 num_contains(SurfaceView::SurfaceState ss, unsigned char is) const;

			std::vector<PCnt> vNumSS;
		};

	public:
		DoFCount() {};
		DoFCount(const GridLevel& gl, ConstSmartPtr<DoFDistributionInfo> spDDInfo);

		/// sums values over all procs (reduced to 'proc', allreduce for -1)
		void sum_values_over_procs(int proc = -1);
		void collect_values(std::vector<uint64>& vNum) const;
		void set_values(const std::vector<uint64>& vNum);

		void add(int fct, int si, SurfaceView::SurfaceState ss, unsigned char is, uint64 numDoF);

		const GridLevel& grid_level() const {return m_gridLevel;}

		uint64 num(int fct, int si, SurfaceView::SurfaceState ss, unsigned char is) const;
		uint64 num_contains(int fct, int si, SurfaceView::SurfaceState ss, unsigned char is) const;

	protected:
		GridLevel m_gridLevel;
		std::vector<std::vector<Cnt> > vvCmpSubset;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_COUNT__ */
