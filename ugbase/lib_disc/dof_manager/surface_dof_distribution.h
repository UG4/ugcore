/*
 * surface_dof_distribution.h
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#ifndef SURFACE_DOF_DISTRIBUTION_H_
#define SURFACE_DOF_DISTRIBUTION_H_

#include "mg_dof_distribution.h"
#include "managing_dof_distribution.h"
#include "lib_disc/domain_traits.h"
#include "lib_disc/function_spaces/local_transfer_interface.h"

namespace ug{

class SurfaceDoFDistribution : public MGDoFDistribution, public ManagingDoFDistribution
{
	typedef MGDoFDistribution base_type;

	public:
	///	constructor
		SurfaceDoFDistribution(SmartPtr<MultiGrid> spMG,
		                       SmartPtr<MGSubsetHandler> spMGSH,
		                       FunctionPattern& fctPatt,
		                       SmartPtr<SurfaceLevelView> spSurfLevelView,
		                       int level, bool bGrouped);

	///	defragments the index set
		void defragment();

	///	redistributes all dofs and resizes associated vectors afterwards.
		virtual void redistribute_dofs();

	///	returns all indices of the element
	///	\{
		void indices(GeometricObject* elem, LocalIndices& ind, bool bHang = false) const{base_type::indices(elem, ind, bHang);}
		void indices(VertexBase* elem, LocalIndices& ind, bool bHang = false) const{base_type::indices(elem, ind, bHang);}
		void indices(EdgeBase* elem, LocalIndices& ind, bool bHang = false) const{base_type::indices(elem, ind, bHang);}
		void indices(Face* elem, LocalIndices& ind, bool bHang = false) const{base_type::indices(elem, ind, bHang);}
		void indices(Volume* elem, LocalIndices& ind, bool bHang = false) const{base_type::indices(elem, ind, bHang);}
	/// \}

	/// get multi indices (Element + Closure of Element)
	/// \{
		size_t multi_indices(GeometricObject* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::multi_indices(elem, fct, ind);}
		size_t multi_indices(VertexBase* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::multi_indices(elem, fct, ind);}
		size_t multi_indices(EdgeBase* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::multi_indices(elem, fct, ind);}
		size_t multi_indices(Face* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::multi_indices(elem, fct, ind);}
		size_t multi_indices(Volume* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::multi_indices(elem, fct, ind);}
	/// \}

	/// get multi indices (only inner part of Element)
	/// \{
		size_t inner_multi_indices(VertexBase* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::inner_multi_indices(elem, fct, ind);}
		size_t inner_multi_indices(EdgeBase* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::inner_multi_indices(elem, fct, ind);}
		size_t inner_multi_indices(Face* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::inner_multi_indices(elem, fct, ind);}
		size_t inner_multi_indices(Volume* elem, size_t fct, std::vector<multi_index_type>& ind) const{return base_type::inner_multi_indices(elem, fct, ind);}
	/// \}

	/// get algebra indices (Element + Closure of Element)
	/// \{
		size_t algebra_indices(VertexBase* elem, std::vector<size_t>& ind) const{return base_type::algebra_indices(elem, ind);}
		size_t algebra_indices(EdgeBase* elem, std::vector<size_t>& ind) const{return base_type::algebra_indices(elem, ind);}
		size_t algebra_indices(Face* elem, std::vector<size_t>& ind) const{return base_type::algebra_indices(elem, ind);}
		size_t algebra_indices(Volume* elem, std::vector<size_t>& ind) const{return base_type::algebra_indices(elem, ind);}
	/// \}

	/// get algebra indices (only inner part of Element)
	/// \{
		size_t inner_algebra_indices(VertexBase* elem, std::vector<size_t>& ind) const{return base_type::inner_algebra_indices(elem,ind);}
		size_t inner_algebra_indices(EdgeBase* elem, std::vector<size_t>& ind) const{return base_type::inner_algebra_indices(elem,ind);}
		size_t inner_algebra_indices(Face* elem, std::vector<size_t>& ind) const{return base_type::inner_algebra_indices(elem,ind);}
		size_t inner_algebra_indices(Volume* elem, std::vector<size_t>& ind) const{return base_type::inner_algebra_indices(elem,ind);}
	/// \}

	protected:
	///	initializes the indices
		void init();

	///	initializes the indices
		template <typename TBaseElem>
		void init();

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
		void defragment(std::vector<std::pair<size_t,size_t> >& vReplaced);

	///	adds indices to created objects
	/**
	 * When an element is inserted into the grid, this function is called an
	 * adds needed indices to the grid object.
	 */
		template <typename TBaseElem>
		inline void obj_created(TBaseElem* obj, GeometricObject* pParent = NULL,
		                        bool replacesParent = false);

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

	public:
		template <typename TElem>
		struct traits
		{
			typedef TElem geometric_object;
			typedef typename SurfaceLevelView::traits<TElem>::iterator iterator;
			typedef typename SurfaceLevelView::traits<TElem>::const_iterator const_iterator;
		};

		template <int dim>
		struct dim_traits
		{
			typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;
			typedef typename SurfaceLevelView::traits<geometric_base_object>::iterator iterator;
			typedef typename SurfaceLevelView::traits<geometric_base_object>::const_iterator const_iterator;
		};

	public:
	///	returns grid level
		GridLevel grid_level() const {return GridLevel(m_level, GridLevel::SURFACE);}

		///////////////////////////////////////
		// Element Access
		///////////////////////////////////////

	/// iterator for elements where dofs are defined
	/// \{
		template <typename TElem>
		typename traits<TElem>::iterator begin() {return m_spSurfLevelView->begin<TElem>();}

		template <typename TElem>
		typename traits<TElem>::iterator end() {return m_spSurfLevelView->end<TElem>();}

		template <typename TElem>
		typename traits<TElem>::const_iterator begin() const {return m_spSurfLevelView->begin<TElem>();}

		template <typename TElem>
		typename traits<TElem>::const_iterator end() const {return m_spSurfLevelView->end<TElem>();}
	///	\}

	/// number of subsets where dofs are defined
		int num_subsets() const {return m_spSurfLevelView->num_subsets();}

	/// iterator for elements where dofs are defined
	/// \{
		template <typename TElem>
		typename traits<TElem>::iterator begin(int si) {return m_spSurfLevelView->begin<TElem>(si);}

		template <typename TElem>
		typename traits<TElem>::iterator end(int si) {return m_spSurfLevelView->end<TElem>(si);}

		template <typename TElem>
		typename traits<TElem>::const_iterator begin(int si) const {return m_spSurfLevelView->begin<TElem>(si);}

		template <typename TElem>
		typename traits<TElem>::const_iterator end(int si) const {return m_spSurfLevelView->end<TElem>(si);}
	///	\}

		///////////////////////////////////////
		// Index Access
		///////////////////////////////////////

	/// return the number of dofs distributed
		size_t num_indices() const {return lev_info().numIndex;}

	/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return lev_info().vNumIndexOnSubset[si];}

	///	returns adjacency graph if available
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	///	renames the indices
		void permute_indices(const std::vector<size_t>& vIndNew);

	///	add a transfer callback
		void add_transfer(SmartPtr<ILocalTransfer> transfer);

	///	add a transfer callback
		void remove_transfer(SmartPtr<ILocalTransfer> transfer);

	///	add a transfer callback
		void clear_transfers();

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
		SmartPtr<SurfaceLevelView> m_spSurfLevelView;

	///	DoF Info
		LevInfo<std::set<size_t> > m_levInfo;

	///	level
		int m_level;

	///	set of invalid indices, still contained in index set
		std::set<size_t> m_sFreeIndex;

	///	list of prolongations
		std::vector<SmartPtr<ILocalTransfer> >m_vProlongation[NUM_GEOMETRIC_BASE_OBJECTS];
		std::vector<SmartPtr<ILocalTransfer> >m_vRestriction[NUM_GEOMETRIC_BASE_OBJECTS];

	protected:
#ifdef UG_PARALLEL
		void create_layouts_and_communicator();

		void create_index_layout(IndexLayout& layout, int keyType);

		template <typename TBaseElem>
		void add_indices_from_layouts(IndexLayout& indexLayout, int keyType);

	protected:
		DistributedGridManager* m_pDistGridMgr;

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		const IndexLayout& master_layout() const {return lev_info().masterLayout;}
		const IndexLayout& slave_layout() const {return lev_info().slaveLayout;}
	/// \}

	/// returns the vertical slave/master index layout
	/// \{
		const IndexLayout& vertical_master_layout() const {return lev_info().verticalMasterLayout;}
		const IndexLayout& vertical_slave_layout() const {return lev_info().verticalSlaveLayout;}
	/// \}

	///	returns communicator
	/// \{
		const pcl::InterfaceCommunicator<IndexLayout>& communicator() const  {return lev_info().communicator;}
		const pcl::ProcessCommunicator& process_communicator() const {return lev_info().processCommunicator;}
	/// \}

	// \TODO: Non-const access should be private or be removed
	public:
	/// returns the horizontal slave/master index layout
	/// \{
		IndexLayout& master_layout(){return lev_info().masterLayout;}
		IndexLayout& slave_layout()	{return lev_info().slaveLayout;}
	/// \}

	/// returns the vertical slave/master index layout
	/// \{
		IndexLayout& vertical_master_layout() {return lev_info().verticalMasterLayout;}
		IndexLayout& vertical_slave_layout()  {return lev_info().verticalSlaveLayout;}
	/// \}

	///	returns communicator
	/// \{
		pcl::InterfaceCommunicator<IndexLayout>& communicator() {return lev_info().communicator;}
		pcl::ProcessCommunicator& process_communicator()	{return lev_info().processCommunicator;}
	/// \}
#endif
};


} // end namespace ug

#endif /* SURFACE_DOF_DISTRIBUTION_H_ */
