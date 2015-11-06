/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__

#include "lib_grid/tools/surface_view.h"
#include "lib_disc/domain_traits.h"
#include "lib_disc/common/local_algebra.h"
#include "dof_index_storage.h"
#include "dof_count.h"

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/algebra_layouts.h"
#endif

namespace ug{

class IGridFunction;

class DoFDistribution : public DoFDistributionInfoProvider
{
	public:
		///	constructor
		DoFDistribution(SmartPtr<MultiGrid> spMG,
		                SmartPtr<MGSubsetHandler> spMGSH,
		                ConstSmartPtr<DoFDistributionInfo> spDDInfo,
		                SmartPtr<SurfaceView> spSurfView,
		                const GridLevel& level, bool bGrouped,
		                SmartPtr<DoFIndexStorage> spDoFIndexStorage = SPNULL);

		/// destructor
		~DoFDistribution();

		///	returns the surface view
		ConstSmartPtr<SurfaceView> surface_view() const {return m_spSurfView;}

		///	returns the multigrid
		/// \{
		SmartPtr<MultiGrid> multi_grid() {return m_spMG;}
		ConstSmartPtr<MultiGrid> multi_grid() const {return m_spMG;}
		/// \}

		///	returns grid level
		const GridLevel& grid_level() const {return m_gridLevel;}

	public:
		template <typename TElem>
		struct traits
		{
			typedef TElem grid_object;
			typedef typename SurfaceView::traits<TElem>::iterator iterator;
			typedef typename SurfaceView::traits<TElem>::const_iterator const_iterator;
		};

		template <int dim>
		struct dim_traits
		{
			typedef typename domain_traits<dim>::grid_base_object grid_base_object;
			typedef typename SurfaceView::traits<grid_base_object>::iterator iterator;
			typedef typename SurfaceView::traits<grid_base_object>::const_iterator const_iterator;
		};

		/// iterator for elements where dofs are defined
		/// \{
		template <typename TElem>
		typename traits<TElem>::iterator begin()
		{return m_spSurfView->begin<TElem>(m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::iterator
		begin(SurfaceView::SurfaceConstants validStates)
		{return m_spSurfView->begin<TElem>(m_gridLevel, validStates);}

		template <typename TElem>
		typename traits<TElem>::iterator end()
		{return m_spSurfView->end<TElem>(m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::iterator
		end(SurfaceView::SurfaceConstants validStates)
		{return m_spSurfView->end<TElem>(m_gridLevel, validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator begin() const
		{return m_spSurfView->begin<TElem>(m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		begin(SurfaceView::SurfaceConstants validStates) const
		{return m_spSurfView->begin<TElem>(m_gridLevel, validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator end() const
		{return m_spSurfView->end<TElem>(m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		end(SurfaceView::SurfaceConstants validStates) const
		{return m_spSurfView->end<TElem>(m_gridLevel, validStates);}
		///	\}

		/// iterator for elements where dofs are defined
		/// \{
		template <typename TElem>
		typename traits<TElem>::iterator begin(int si)
		{return m_spSurfView->begin<TElem>(si, m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::iterator
		begin(int si, SurfaceView::SurfaceConstants validStates)
		{return m_spSurfView->begin<TElem>(si, m_gridLevel, validStates);}

		template <typename TElem>
		typename traits<TElem>::iterator end(int si)
		{return m_spSurfView->end<TElem>(si, m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::iterator
		end(int si, SurfaceView::SurfaceConstants validStates)
		{return m_spSurfView->end<TElem>(si, m_gridLevel, validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator begin(int si) const
		{return m_spSurfView->begin<TElem>(si, m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		begin(int si, SurfaceView::SurfaceConstants validStates) const
		{return m_spSurfView->begin<TElem>(si, m_gridLevel, validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator end(int si) const
		{return m_spSurfView->end<TElem>(si, m_gridLevel, defaultValidSurfState());}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		end(int si, SurfaceView::SurfaceConstants validStates) const
		{return m_spSurfView->end<TElem>(si, m_gridLevel, validStates);}
		///	\}

		/// returns the default valid surface state
		SurfaceView::SurfaceConstants defaultValidSurfState() const;

		///	returns the adjacent elements
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
								TElem* elem, bool clearContainer = true) const{
//			if(grid_level().is_level())
//				CollectAssociated(vAssElem, *m_pMG, elem, clearContainer);
//			else{
				m_spSurfView->collect_associated(vAssElem, elem, grid_level(), clearContainer);
//			}
		}

		/// returns if the grid object is part of the dof distribution
		template <class TGeomObj>
		bool is_contained(TGeomObj* obj) const{
			return m_spSurfView->is_contained(obj, grid_level());
		}


	public:
		///	returns if dofs are grouped
		bool grouped() const {return m_bGrouped;}

		/// return the number of dofs distributed
		size_t num_indices() const {return m_numIndex;}

		/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return m_vNumIndexOnSubset[si];}

	public:
		/// extracts all indices of the element (sorted)
		/**
		 * All Indices of the element (including the subelements) are extracted
		 * and stored in the LocalIndices structure. The order of the indices
		 * is sorted, i.e. the dofs are provided as specified in the local
		 * dof set of the local finite element trial space.
		 * If bHang is set to true, also the DoFs on the Constained Objects
		 * belonging to the constraining Subelements are extracted and added at
		 * the end of the indices.
		 *
		 * \param[in]		elem		the element
		 * \param[out]		ind			Local indices
		 * \param[in]		bHang		flag if extracting of constrained dofs required
		 */
		/// \{
		void indices(GridObject* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(Vertex* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(Edge* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(Face* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(Volume* elem, LocalIndices& ind, bool bHang = false) const;
		/// \}

		/// extracts all multiindices for a function (sorted)
		/**
		 * All Multi-Indices of a function living on the element (including the
		 * subelements) are extracted and stored in a std::vector. The order of
		 * the indices is sorted, i.e. the dofs are provided as specified in the
		 * local dof set of the local finite element trial space.
		 * If bHang is set to true, also the DoFs on the Constrained Objects
		 * belonging to the constraining Subelements are extracted and added at
		 * the end of the indices.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[in]		fct			the function
		 * \param[out]		ind			vector of multi indices
		 * \param[in]		bHang		flag if extracting of constrained dofs required
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		/// \{
		size_t dof_indices(GridObject* elem, size_t fct, std::vector<DoFIndex>& ind,
		                   bool bHang = false, bool bClear = true) const;
		size_t dof_indices(Vertex* elem, size_t fct, std::vector<DoFIndex>& ind,
		                   bool bHang = false, bool bClear = true) const;
		size_t dof_indices(Edge* elem, size_t fct, std::vector<DoFIndex>& ind,
		                   bool bHang = false, bool bClear = true) const;
		size_t dof_indices(Face* elem, size_t fct, std::vector<DoFIndex>& ind,
		                   bool bHang = false, bool bClear = true) const;
		size_t dof_indices(Volume* elem, size_t fct, std::vector<DoFIndex>& ind,
		                   bool bHang = false, bool bClear = true) const;
		/// \}

		/// extracts all multiindices of a function in the inner (sorted)
		/**
		 * All Multi-Indices of a function living on the element (including the
		 * subelements) are extracted and stored in a std::vector. The order of
		 * the indices is sorted, i.e. the dofs are provided as specified in the
		 * local dof set of the local finite element trial space.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[in]		fct			the function
		 * \param[out]		ind			vector of multi indices
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		/// \{
		size_t inner_dof_indices(GridObject* elem, size_t fct, std::vector<DoFIndex>& ind,
		                         bool bClear = true) const;
		size_t inner_dof_indices(Vertex* elem, size_t fct, std::vector<DoFIndex>& ind,
		                         bool bClear = true) const;
		size_t inner_dof_indices(Edge* elem, size_t fct, std::vector<DoFIndex>& ind,
		                         bool bClear = true) const;
		size_t inner_dof_indices(Face* elem, size_t fct, std::vector<DoFIndex>& ind,
		                         bool bClear = true) const;
		size_t inner_dof_indices(Volume* elem, size_t fct, std::vector<DoFIndex>& ind,
		                         bool bClear = true) const;
		/// \}

		/// extracts all algebra indices of an element (not sorted)
		/**
		 * All Algebra-Indices of the element (including the subelements) are
		 * extracted and stored in a std::vector. The order of the indices is
		 * not sorted, no constrained DoFs are extracted.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[out]		ind			vector of algebra indices
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		/// \{
		size_t algebra_indices(GridObject* elem,	std::vector<size_t>& ind,
		                       bool bClear = true) const;
		size_t algebra_indices(Vertex* elem, std::vector<size_t>& ind,
		                       bool bClear = true) const;
		size_t algebra_indices(Edge* elem, std::vector<size_t>& ind,
		                       bool bClear = true) const;
		size_t algebra_indices(Face* elem, std::vector<size_t>& ind,
		                       bool bClear = true) const;
		size_t algebra_indices(Volume* elem, std::vector<size_t>& ind,
		                       bool bClear = true) const;
		/// \}

		size_t inner_algebra_indices_for_fct(GridObject* elem, std::vector<size_t>& ind,
							bool bClear, int fct) const;

		/// extracts all algebra indices in the inner of the element (not sorted)
		/**
		 * All Algebra-Indices of the element (excluding the subelements) are
		 * extracted and stored in a std::vector. The order of the indices is
		 * not sorted, no constrained DoFs are extracted.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[out]		ind			vector of algebra indices
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		/// \{
		size_t inner_algebra_indices(GridObject* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;
		size_t inner_algebra_indices(Vertex* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;
		size_t inner_algebra_indices(Edge* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;
		size_t inner_algebra_indices(Face* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;
		size_t inner_algebra_indices(Volume* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;
		/// \}

	protected:
		template <typename TBaseElem>
		void _indices(TBaseElem* elem, LocalIndices& ind, bool bHang = false) const;

		template<typename TBaseElem>
		size_t _dof_indices(TBaseElem* elem, size_t fct,
		                     std::vector<DoFIndex>& ind,
		                     bool bHang = false, bool bClear = true) const;

		template<typename TBaseElem>
		size_t _inner_dof_indices(TBaseElem* elem, size_t fct,
		                           std::vector<DoFIndex>& ind,
		                           bool bClear = true) const;

		template<typename TBaseElem>
		size_t _algebra_indices(TBaseElem* elem,	std::vector<size_t>& ind,
		                       bool bClear = true) const;

		template<typename TBaseElem>
		size_t _inner_algebra_indices(TBaseElem* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;

	///	returns indices, that can be changed on the element
		template <typename TBaseElem>
		void changable_indices(std::vector<size_t>& vIndex,
		                       const std::vector<TBaseElem*>& vElem) const;

		///	extracts the indices of the vertices
		template<typename TBaseElem>
		void indices_on_vertex(TBaseElem* elem, const ReferenceObjectID roid,
		                       LocalIndices& ind,
		                       const Grid::SecureVertexContainer& vElem) const;

		///	extract dofs on constrained objects
		template <typename TConstraining, typename TConstrained, typename TBaseElem>
		void constrained_vertex_indices(LocalIndices& ind,
		                         const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const;

		template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
		void constrained_edge_indices(TBaseElem* elem,LocalIndices& ind,
		                         const typename Grid::traits<TSubElem>::secure_container& vSubElem) const;

		template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
		void constrained_face_indices(TBaseElem* elem,LocalIndices& ind,
		                         const typename Grid::traits<TSubElem>::secure_container& vSubElem) const;

		// sorts indices on constrained edges
		template <typename TBaseElem,typename TConstraining, typename TConstrained>
		void sort_constrained_edges(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const;

		// sorts indices on constrained faces
		template <typename TBaseElem,typename TConstraining, typename TConstrained>
		void sort_constrained_faces(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const;

		/// extracts the indices of the subelement of an element
		template<typename TBaseElem, typename TSubBaseElem>
		void indices(TBaseElem* elem, const ReferenceObjectID roid,
		             LocalIndices& ind,
		             const typename Grid::traits<TSubBaseElem>::secure_container& vElem) const;

		/// extracts the indices of a subelement of an element
		template<typename TBaseElem, typename TSubBaseElem>
		void dof_indices(TBaseElem* elem, const ReferenceObjectID roid,
		                   size_t fct, std::vector<DoFIndex>& ind,
		                   const typename Grid::traits<TSubBaseElem>::secure_container& vElem) const;

		/// multi indices on constrained vertices
		template <typename TConstraining, typename TConstrained, typename TBaseElem>
		void constrained_vertex_dof_indices(size_t fct,std::vector<DoFIndex>& ind,
									const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const;

		/// multi indices on constrained edges
		template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
		void constrained_edge_dof_indices(TBaseElem* elem,size_t fct,std::vector<DoFIndex>& ind,
									const typename Grid::traits<TSubElem>::secure_container& vSubElem) const;

		/// multi indices on constrained faces
		template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
		void constrained_face_dof_indices(TBaseElem* elem,size_t fct,std::vector<DoFIndex>& ind,
									const typename Grid::traits<TSubElem>::secure_container& vSubElem) const;

		/// adds all algebra indices of an geom object to the LocalIndices
		template <typename TBaseElem>
		size_t extract_inner_algebra_indices(TBaseElem* elem,
		                                     std::vector<size_t>& ind) const;

		///	adds all algebra indices of a set of geometric objects
		template<typename TBaseElem>
		void extract_inner_algebra_indices(const typename Grid::traits<TBaseElem>::secure_container& vElem,
		                                   std::vector<size_t>& ind) const;

	protected:
		/// adds indices to a geometric object
		template <typename TBaseObject>
		void add(TBaseObject* obj, const ReferenceObjectID roid, const int si);

		///	checks that subset assignment is ok
		void check_subsets();

		/// returns the index assigned to a grid object
		/// \{
		template <typename TElem>
		inline size_t& obj_index(TElem* obj) {return m_spDoFIndexStorage->obj_index(obj);}

		template <typename TElem>
		inline const size_t& obj_index(TElem* obj) const {return m_spDoFIndexStorage->obj_index(obj);}
		/// \}

	protected:
		///	grouping
		bool m_bGrouped;

		///	Multi Grid
		SmartPtr<MultiGrid> m_spMG;
		MultiGrid* m_pMG;

		/// Subset Handler
		SmartPtr<MGSubsetHandler> m_spMGSH;

		///	MultiGrid Subset Handler
		SmartPtr<SurfaceView> m_spSurfView;

		///	Grid level
		GridLevel m_gridLevel;

		/// DoF-Index Memory Storage
		SmartPtr<DoFIndexStorage> m_spDoFIndexStorage;

	protected:
		/// number of distributed indices on whole domain
		size_t m_numIndex;

		/// number of distributed indices on each subset
		std::vector<size_t> m_vNumIndexOnSubset;

	public:
		/// returns the connections
		void get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

		///	renames the indices
		void permute_indices(const std::vector<size_t>& vIndNew);

		///	initializes the indices
		void reinit();

	protected:
		///	initializes the indices
		template <typename TBaseElem>
		void reinit();

		template <typename TBaseElem>
		void permute_indices(const std::vector<size_t>& vNewInd);

		template <typename TBaseElem>
		void get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	public:
	///	registers a grid function for adaptation management
		void manage_grid_function(IGridFunction& gridFct);

	///	unregisters a grid function for adaptation management
		void unmanage_grid_function(IGridFunction& gridFct);

	public:
	///	permutes values in managed functions, if indices permuted
		void permute_values(const std::vector<size_t>& vIndNew);

	///	swaps values in managed functions, if indices swapped
		void copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                                         bool bDisjunct);

	///	changes values in managed functions, number of indices changed
		void resize_values(size_t newSize);

	protected:
	///	managed grid functions
	/**
	 * This vector holds a pointer to all grid functions, that should be
	 * managed (i.e. adapted), when the dof distribution is changed.
	 * NOTE:	No SmartPtr is used here on purpose. The GridFunction stores a
	 * 			SmartPtr to the DoFDistribution. If we would use a SmartPtr
	 * 			here those objects would reference each other and would never
	 * 			be deleted.
	 */
		std::vector<IGridFunction*> m_vpGridFunction;

#ifdef UG_PARALLEL
		public:
		///	returns algebra layouts
		///	\{
		ConstSmartPtr<AlgebraLayouts> layouts() const 	{return spAlgebraLayouts;}
		SmartPtr<AlgebraLayouts> layouts() 				{return spAlgebraLayouts;}
		///	\}

	protected:
		///	algebra layouts
		SmartPtr<AlgebraLayouts> spAlgebraLayouts;

		void reinit_layouts_and_communicator();

		void reinit_index_layout(IndexLayout& layout, int keyType);

		template <typename TBaseElem>
		void add_indices_from_layouts(IndexLayout& indexLayout, int keyType);
#endif

	public:
		DoFCount dof_count() const;

	protected:
		template <typename TBaseElem>
		void sum_dof_count(DoFCount& cnt) const;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__ */
