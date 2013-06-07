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
#include "managing_dof_distribution.h"
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

	///	defragments the index set
		void defragment();

	///	redistributes all dofs and resizes associated vectors afterwards.
		virtual void redistribute_dofs(bool bReinit);

	protected:
	///	initializes the indices
		void init();

	///	initializes the indices
		template <typename TBaseElem>
		void init();

	///	initializes the indices and permute values
		void reinit();

	///	initializes the indices and permute values
		template <typename TBaseElem>
		void reinit(std::vector<std::pair<size_t,size_t> >& vReplaced);

	///	removes holes in the index set
	/**
	 * This method removes holes in the index set such that the index set is
	 * contiguous. Therefore, free indices are replaced by those at the end
	 * of the index set. The replacement is stored in the vReplaced vector (and
	 * may be used to adjust associated data, e.g. a grid vector).
	 *
	 * \param[in,out]	vReplaced	vector with all pairs of replacements
	 */
		template <typename TBaseElem>
		bool defragment(std::vector<std::pair<size_t,size_t> >& vReplaced);

	/**
	 * Iterate over all elements and adds those whose
	 * dof-entry has not yet been assigned (whose index equals NOT_YET_ASSIGNED)*/
		template <typename TBaseElem>
		void add_unassigned_elements();

	///	removes dof-indices from ghosts (those are no longer contained in the surface grid).
		template <typename TBaseElem>
		void remove_ghost_entries();

	///	called by base class when parallel redistribution is done.
		virtual void parallel_redistribution_ended();


	///	adds indices to created objects
	/**
	 * When an element is inserted into the grid, this function is called an
	 * adds needed indices to the grid object.
	 */
		template <typename TBaseElem>
		inline void obj_created(TBaseElem* obj, GeometricObject* pParent = NULL,
		                        bool replacesParent = false);

	///	adds side elements of the given element, if those havn't been added already
		template <typename TElem>
		void add_unassigned_sides(TElem* e);

	///	removes indices, when a grid element is removed
	/**
	 * When a grid element is removed from the grid, this function is called
	 * and takes care about the indices. All indices associated with the element
	 * are removed and stored in a free index container, counters are adjusted.
	 * In general the removal of a grid element will lead to holes in the index
	 * set. Those can be removed by calling defragment.
	 *
	 * \param[in]		obj		grid object that will be removed
	 */
		template <typename TBaseElem>
		inline void obj_to_be_erased(TBaseElem* obj, TBaseElem* replacedBy = NULL);

	public:
		/// grid observer callbacks
		/// \{
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL, bool replacesParent = false);
		virtual void edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent = NULL, bool replacesParent = false);
		virtual void face_created(Grid* grid, Face* f, GeometricObject* pParent = NULL, bool replacesParent = false);
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL, bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy = NULL);
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e, EdgeBase* replacedBy = NULL);
		virtual void face_to_be_erased(Grid* grid, Face* f, Face* replacedBy = NULL);
		virtual void volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy = NULL);
		/// \}


	/// return the number of dofs distributed
		size_t num_indices() const {return lev_info().numIndex;}

	/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return lev_info().vNumIndexOnSubset[si];}

	///	returns adjacency graph if available
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	///	renames the indices
		void permute_indices(const std::vector<size_t>& vIndNew);
		
	/// set redistribution	
		void set_redistribution(bool bRedistribute){m_bRedistribute = bRedistribute;}
		
	///	returns parent != NULL, if it is a copy and shadowed by elem
		template <typename TBaseElem>
		TBaseElem* parent_if_shadowed_copy(TBaseElem* elem) const;

	protected:
		/// permutes the indices for an base element type
		template <typename TBaseElem>
		void permute_indices(const std::vector<size_t>& vNewInd);

		template <typename TBaseElem>
		void get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

		LevInfo<std::set<size_t> >& lev_info() {return m_levInfo;}
		const LevInfo<std::set<size_t> >& lev_info() const {return m_levInfo;}

	protected:
	///	MultiGrid Subset Handler
		SmartPtr<SurfaceView> m_spSurfView;

	///	DoF Info
		LevInfo<std::set<size_t> > m_levInfo;

	///	level
		int m_level;

	///	set of invalid indices, still contained in index set
		std::set<size_t> m_sFreeIndex;
		
		bool m_bRedistribute;

	protected:
#ifdef UG_PARALLEL
		void create_layouts_and_communicator();

		void create_index_layout(IndexLayout& layout, int keyType);

		template <typename TBaseElem>
		void add_indices_from_layouts(IndexLayout& indexLayout, int keyType);

	protected:
		DistributedGridManager* m_pDistGridMgr;
#endif

#ifdef UG_PARALLEL
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
