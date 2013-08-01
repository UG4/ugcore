/*
 * surface_dof_distribution.h
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__SURFACE_DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__SURFACE_DOF_DISTRIBUTION__

#include "dof_distribution.h"
#include "mg_dof_distribution.h"
#include "lib_disc/domain_traits.h"

namespace ug{

class SurfaceDoFDistribution : public MGDoFDistribution,
							   public DoFDistribution
{
	public:
	///	constructor
		SurfaceDoFDistribution(SmartPtr<MultiGrid> spMG,
		                       SmartPtr<MGSubsetHandler> spMGSH,
		                       ConstSmartPtr<DoFDistributionInfo> spDDInfo,
		                       SmartPtr<SurfaceView> spSurfView,
		                       int level, bool bGrouped);

	protected:
	///	initializes the indices
		template <typename TBaseElem>
		void reinit();

	public:
	///	initializes the indices and permute values
		void reinit();

	/// return the number of dofs distributed
		size_t num_indices() const {return lev_info().numIndex;}

	/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return lev_info().vNumIndexOnSubset[si];}

	///	returns adjacency graph if available
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	///	renames the indices
		void permute_indices(const std::vector<size_t>& vIndNew);
		
	///	returns parent != NULL, if it is a copy and shadowed by elem
		template <typename TBaseElem>
		TBaseElem* parent_if_shadowed_copy(TBaseElem* elem) const;

	protected:
		/// permutes the indices for an base element type
		template <typename TBaseElem>
		void permute_indices(const std::vector<size_t>& vNewInd);

		template <typename TBaseElem>
		void get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

		LevInfo& lev_info() {return m_levInfo;}
		const LevInfo& lev_info() const {return m_levInfo;}

	protected:
	///	MultiGrid Subset Handler
		SmartPtr<SurfaceView> m_spSurfView;

	///	DoF Info
		LevInfo m_levInfo;

	///	level
		int m_level;

#ifdef UG_PARALLEL
	protected:
		void create_layouts_and_communicator();

		void create_index_layout(IndexLayout& layout, int keyType);

		template <typename TBaseElem>
		void add_indices_from_layouts(IndexLayout& indexLayout, int keyType);

	public:
	///	returns the algebra layouts
		ConstSmartPtr<AlgebraLayouts> layouts() const {return lev_info().layouts();}

	// \TODO: Non-const access should be private or be removed
	public:
	///	returns the algebra layouts
		SmartPtr<AlgebraLayouts> layouts() {return lev_info().layouts();}
#endif
};


template <typename TBaseElem>
TBaseElem* SurfaceDoFDistribution::
parent_if_shadowed_copy(TBaseElem* elem) const
{
	TBaseElem* parent = parent_if_copy(elem);
	if(parent && m_spSurfView->is_shadowed(parent))
		return parent;
	return NULL;
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__SURFACE_DOF_DISTRIBUTION__ */
