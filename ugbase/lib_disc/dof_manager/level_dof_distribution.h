/*
 * level_dof_distribution.h
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__LEVEL_DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__LEVEL_DOF_DISTRIBUTION__

#include "dof_distribution.h"
#include "mg_dof_distribution.h"
#include "lib_disc/domain_traits.h"

namespace ug{

class LevelMGDoFDistribution : public MGDoFDistribution
{
	friend class LevelDoFDistribution;

	public:
	///	constructor
		LevelMGDoFDistribution(SmartPtr<MultiGrid> spMG,
		                       SmartPtr<MGSubsetHandler> spMGSH,
		                       ConstSmartPtr<DoFDistributionInfo> spDDInfo,
		                       bool bGrouped);

	///	register to manage a DoFDistribution for a given level
	/**	Those dof-distributions will be informed, whenever redistribute_dofs is
	 * executed.
	 * This method is e.g. used by LevelDoFDistribution.*/
		void manage_dof_distribution(DoFDistribution* mdd, int lvl);

	public:
	///	initializes the indices
		void reinit();

	protected:
	///	initializes the indices
		template <typename TBaseElem>
		void reinit();

	protected:
		///	returns the number of indices on whole level
		size_t num_indices(const int lev) const {
			return m_vLev[lev].numIndex;
		}

		///	returns the number of indices on a level and a subset
		size_t num_indices(const int lev, const int si) const {
			return m_vLev[lev].vNumIndexOnSubset[si];
		}

		///	permutes the indices on a grid level
		void permute_indices(const std::vector<size_t>& vNewInd, int lev);

		/// permutes the indices on a grid level for an base element type
		template <typename TBaseElem>
		void permute_indices(const std::vector<size_t>& vNewInd, int lev);

	///	adjusts storage for requested level
		void level_required(int level);

	///	informations for each level
		std::vector<LevInfo> m_vLev;

		LevInfo& lev_info(int lev) {return m_vLev[lev];}
		const LevInfo& lev_info(int lev) const {return m_vLev[lev];}

		std::vector<DoFDistribution*>	m_managingDoFDists;

#ifdef UG_PARALLEL
	///	returns the algebra layouts
		ConstSmartPtr<AlgebraLayouts> layouts(const int lev) const {return lev_info(lev).layouts();}

	///	returns the algebra layouts
		SmartPtr<AlgebraLayouts> layouts(const int lev) {return lev_info(lev).layouts();}

		void create_layouts_and_communicator(int l);

		void create_index_layout(IndexLayout& layout, InterfaceNodeTypes keyType, int l);

		template <typename TBaseElem>
		void add_indices_from_layouts(IndexLayout& indexLayout, InterfaceNodeTypes keyType, int l);
#endif
};



class LevelDoFDistribution :  public DoFDistribution
{
	public:
	///	constructor
		LevelDoFDistribution(SmartPtr<LevelMGDoFDistribution> spLevMGDD,
		                     SmartPtr<SurfaceView> spSurfView,
		                     int level);

		virtual ~LevelDoFDistribution();

	/// return the number of dofs distributed
		size_t num_indices() const {return m_spMGDD->num_indices(grid_level().level());}

	/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return m_spMGDD->num_indices(grid_level().level(), si);}

	///	returns adjacency graph if available
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	///	renames the indices
		void permute_indices(const std::vector<size_t>& vIndNew);

#ifdef UG_PARALLEL
	public:
	///	returns the algebra layouts
		ConstSmartPtr<AlgebraLayouts> layouts() const {return m_spMGDD->layouts(grid_level().level());}

	// \TODO: Non-const access should be private or be removed
	public:
	///	returns the algebra layouts
		SmartPtr<AlgebraLayouts> layouts() {return m_spMGDD->layouts(grid_level().level());}
#endif

	protected:
		template <typename TBaseElem>
		void get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	protected:
	///	MultiGrid Level DoF Distribution
		SmartPtr<LevelMGDoFDistribution> m_spMGDD;

	///	MultiGrid Subset Handler
		SmartPtr<SurfaceView> m_spSurfView;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__LEVEL_DOF_DISTRIBUTION__ */
