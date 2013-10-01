/*
 * level_dof_distribution.h
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__

#include "lib_grid/tools/surface_view.h"
#include "lib_disc/domain_traits.h"
#include "mg_dof_distribution.h"

namespace ug{

class IGridFunction;

class DoFDistribution : virtual public DoFDistributionInfoProvider
{
	public:
	///	constructor
		DoFDistribution(MGDoFDistribution& rMGDD,
		                SmartPtr<SurfaceView> spSurfView,
		                const GridLevel& level)
			: DoFDistributionInfoProvider(rMGDD.dof_distribution_info()),
			  m_rMGDD(rMGDD), m_spSurfView(spSurfView), m_level(level)
		{}

	/// destructor
		virtual ~DoFDistribution() {}

	///	returns the surface view
		ConstSmartPtr<SurfaceView> surface_view() const {return m_spSurfView;}

	///	returns the multi grid
		const MultiGrid& multi_grid() const {return *m_spSurfView->subset_handler()->multi_grid();}

	///	returns grid level
		const GridLevel& grid_level() const {return m_level;}

		///	returns if dofs are grouped
		bool grouped() const {return m_rMGDD.grouped();}

		///	returns blocksize
		std::string blocksize() const {return m_rMGDD.blocksize();}

		/// return the number of dofs distributed
		virtual size_t num_indices() const = 0;

		/// return the number of dofs distributed on subset si
		virtual size_t num_indices(int si) const = 0;

		/// returns the connections
		virtual bool get_connections(std::vector<std::vector<size_t> >& vvConnection) const = 0;

		///	renames the indices
		virtual void permute_indices(const std::vector<size_t>& vIndNew) = 0;

	public:
		template <typename TElem>
		struct traits
		{
			typedef TElem geometric_object;
			typedef typename SurfaceView::traits<TElem>::iterator iterator;
			typedef typename SurfaceView::traits<TElem>::const_iterator const_iterator;
		};

		template <int dim>
		struct dim_traits
		{
			typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;
			typedef typename SurfaceView::traits<geometric_base_object>::iterator iterator;
			typedef typename SurfaceView::traits<geometric_base_object>::const_iterator const_iterator;
		};

	/// iterator for elements where dofs are defined
	/// \{
		template <typename TElem>
		typename traits<TElem>::iterator
		begin(SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING)
		{return m_spSurfView->begin<TElem>(m_level, validSurfStates);}

		template <typename TElem>
		typename traits<TElem>::iterator
		end(SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING)
		{return m_spSurfView->end<TElem>(m_level, validSurfStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		begin(SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING) const
		{return m_spSurfView->begin<TElem>(m_level, validSurfStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		end(SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING) const
		{return m_spSurfView->end<TElem>(m_level, validSurfStates);}
	///	\}

	/// iterator for elements where dofs are defined
	/// \{
		template <typename TElem>
		typename traits<TElem>::iterator
		begin(int si, SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING)
		{return m_spSurfView->begin<TElem>(si, m_level, validSurfStates);}

		template <typename TElem>
		typename traits<TElem>::iterator
		end(int si, SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING)
		{return m_spSurfView->end<TElem>(si, m_level, validSurfStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		begin(int si, SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING) const
		{return m_spSurfView->begin<TElem>(si, m_level, validSurfStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		end(int si, SurfaceView::SurfaceConstants validSurfStates = SurfaceView::SURFACE_AND_SHADOWING) const
		{return m_spSurfView->end<TElem>(si, m_level, validSurfStates);}
	///	\}

	///	returns the adjacend elements
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
								TElem* elem, bool clearContainer = true) const{
			if(m_level.type() == GridLevel::LEVEL)
				CollectAssociated(vAssElem, *m_rMGDD.multi_grid(), elem, clearContainer);
			else{
				m_spSurfView->collect_associated(vAssElem, elem, clearContainer);
			}
		}

	public:
	///	returns all indices of the element
	///	\{
		void indices(GeometricObject* elem, LocalIndices& ind, bool bHang = false) const{m_rMGDD.indices(elem,ind, bHang);}
		void indices(VertexBase* elem, LocalIndices& ind, bool bHang = false) const{m_rMGDD.indices(elem,ind, bHang);}
		void indices(EdgeBase* elem, LocalIndices& ind, bool bHang = false) const{m_rMGDD.indices(elem,ind, bHang);}
		void indices(Face* elem, LocalIndices& ind, bool bHang = false) const{m_rMGDD.indices(elem,ind, bHang);}
		void indices(Volume* elem, LocalIndices& ind, bool bHang = false) const{m_rMGDD.indices(elem,ind, bHang);}
	/// \}

	/// get multi indices (Element + Closure of Element)
	/// \{
		size_t dof_indices(GeometricObject* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const{return m_rMGDD.dof_indices(elem,fct,ind,bHang,bClear);}
		size_t dof_indices(VertexBase* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const{return m_rMGDD.dof_indices(elem,fct,ind,bHang,bClear);}
		size_t dof_indices(EdgeBase* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const{return m_rMGDD.dof_indices(elem,fct,ind,bHang,bClear);}
		size_t dof_indices(Face* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const{return m_rMGDD.dof_indices(elem,fct,ind,bHang,bClear);}
		size_t dof_indices(Volume* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const{return m_rMGDD.dof_indices(elem,fct,ind,bHang,bClear);}
	/// \}

	/// get multi indices (only inner part of Element)
	/// \{
		size_t inner_dof_indices(GeometricObject* elem, size_t fct, std::vector<DoFIndex>& ind, bool bClear = true) const{return m_rMGDD.inner_dof_indices(elem,fct,ind,bClear);}
		size_t inner_dof_indices(VertexBase* elem, size_t fct, std::vector<DoFIndex>& ind, bool bClear = true) const{return m_rMGDD.inner_dof_indices(elem,fct,ind,bClear);}
		size_t inner_dof_indices(EdgeBase* elem, size_t fct, std::vector<DoFIndex>& ind, bool bClear = true) const{return m_rMGDD.inner_dof_indices(elem,fct,ind,bClear);}
		size_t inner_dof_indices(Face* elem, size_t fct, std::vector<DoFIndex>& ind, bool bClear = true) const{return m_rMGDD.inner_dof_indices(elem,fct,ind,bClear);}
		size_t inner_dof_indices(Volume* elem, size_t fct, std::vector<DoFIndex>& ind, bool bClear = true) const{return m_rMGDD.inner_dof_indices(elem,fct,ind,bClear);}
	/// \}

	/// get algebra indices (Element + Closure of Element)
	/// \{
		size_t algebra_indices(GeometricObject* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.algebra_indices(elem,ind,bClear);}
		size_t algebra_indices(VertexBase* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.algebra_indices(elem,ind,bClear);}
		size_t algebra_indices(EdgeBase* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.algebra_indices(elem,ind,bClear);}
		size_t algebra_indices(Face* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.algebra_indices(elem,ind,bClear);}
		size_t algebra_indices(Volume* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.algebra_indices(elem,ind,bClear);}
	/// \}

	/// get algebra indices (only inner part of Element)
	/// \{
		size_t inner_algebra_indices(GeometricObject* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.inner_algebra_indices(elem,ind,bClear);}
		size_t inner_algebra_indices(VertexBase* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.inner_algebra_indices(elem,ind,bClear);}
		size_t inner_algebra_indices(EdgeBase* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.inner_algebra_indices(elem,ind,bClear);}
		size_t inner_algebra_indices(Face* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.inner_algebra_indices(elem,ind,bClear);}
		size_t inner_algebra_indices(Volume* elem, std::vector<size_t>& ind, bool bClear = true) const{return m_rMGDD.inner_algebra_indices(elem,ind,bClear);}
	/// \}

#ifdef UG_PARALLEL
	public:
	///	returns the algebra layouts
		virtual ConstSmartPtr<AlgebraLayouts> layouts() const = 0;

	public:
	///	returns the algebra layouts
		virtual SmartPtr<AlgebraLayouts> layouts() = 0;
#endif

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

	protected:
	///	MultiGrid Level DoF Distribution
		MGDoFDistribution& m_rMGDD;

	///	MultiGrid Subset Handler
		SmartPtr<SurfaceView> m_spSurfView;

	///	Grid level
		GridLevel m_level;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__ */
