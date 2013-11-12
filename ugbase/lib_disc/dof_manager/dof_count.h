/*
 * dof_count.h
 *
 *  Created on: 06.11.2013
 *      Author: andreasvogel
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
		static const byte ES_MAX = (ES_V_MASTER | ES_H_MASTER | ES_H_SLAVE | ES_V_SLAVE);
		static const byte SS_MAX = SurfaceView::ALL;

		static const int ALL_FCT = -1;
		static const int ALL_SUBSET = -1;
		static const byte ALL_ES = ES_MAX + 1;
		static const byte ALL_SS = SS_MAX + 1;
		static const byte UNIQUE_ES = ES_MAX + 2;
		static const byte UNIQUE_SS = SS_MAX + 2;

	protected:
		struct Cnt{

			Cnt();
			void add(uint64 num, SurfaceView::SurfaceState ss, byte is);
			void collect_values(std::vector<uint64>& vNum) const;
			void set_values(const std::vector<uint64>& vNum, size_t& cnt);

			struct PCnt {

				PCnt();

				// dofs exactly matching this state
				uint64 num(byte is) const;

				// dofs that contain this state
				uint64 num_contains(byte is) const;

				void collect_values(std::vector<uint64>& vNum) const;
				void set_values(const std::vector<uint64>& vNum, size_t& cnt);

				std::vector<uint64> vNumIS;
			};

			uint64 num(SurfaceView::SurfaceState ss, byte is) const;
			uint64 num_contains(SurfaceView::SurfaceState ss, byte is) const;

			std::vector<PCnt> vNumSS;
		};

	public:
		DoFCount() {};
		DoFCount(const GridLevel& gl, ConstSmartPtr<DoFDistributionInfo> spDDInfo);

		/// sums values over all procs (reduced to 'proc', allreduce for -1)
		void sum_values_over_procs(int proc = -1);
		void collect_values(std::vector<uint64>& vNum) const;
		void set_values(const std::vector<uint64>& vNum);

		void add(int fct, int si, SurfaceView::SurfaceState ss, byte is, uint64 numDoF);

		const GridLevel& grid_level() const {return m_gridLevel;}

		uint64 num(int fct, int si, SurfaceView::SurfaceState ss, byte is) const;
		uint64 num_contains(int fct, int si, SurfaceView::SurfaceState ss, byte is) const;

	protected:
		GridLevel m_gridLevel;
		std::vector<std::vector<Cnt> > vvCmpSubset;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_COUNT__ */
