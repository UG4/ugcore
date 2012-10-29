// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)

#ifndef __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW__
#define __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW__

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif

#include <boost/iterator/iterator_facade.hpp>
#include "lib_grid/multi_grid.h"
#include "subset_handler_multi_grid.h"
#include "bool_marker.h"

namespace ug{

/** \ingroup lib_grid_tools
 *  \{ */

///	Represents the surface view of a multi-grid hierarchy.
/**	The surface of a multi-grid hierarchy consists of all elements, which
 * do not have children.
 * The surface view allows to iterate over the subsets of such a surface.
 * The surface_level_begin(...) and surface_level_end(...) methods even allow
 * to ignore all levels above a specified level, resulting in a new surface-view
 * in which only surface elements up to the specified level and the elements in
 * the specified level are regarded as surface-view-elements.
 */
class SurfaceView
{
	public:
		template <class TElem> class SurfaceViewElementIterator;
		template <class TElem> class ConstSurfaceViewElementIterator;

	public:
		template <class TElem>
		struct traits{
			typedef SurfaceViewElementIterator<TElem>		iterator;
			typedef ConstSurfaceViewElementIterator<TElem>	const_iterator;
		};

	public:
		SurfaceView(SmartPtr<MGSubsetHandler> spMGSH,
#ifdef UG_PARALLEL
		            DistributedGridManager* pDistGridMgr,
#endif
		            bool adaptiveMG = true);

	///	returns underlying subset handler
		inline SmartPtr<MGSubsetHandler> subset_handler();

	///	returns underlying subset handler
		inline ConstSmartPtr<MGSubsetHandler> subset_handler() const;

	///	returns if multigrid is adaptive
		inline bool is_adaptive() const;

	///	number of subsets
		int num_subsets() const {return m_spMGSH->num_subsets();}

	///	iterators over whole surface grid of toplevel
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator	begin();

		template <class TElem>
		typename traits<TElem>::iterator	end();

		template <class TElem>
		typename traits<TElem>::const_iterator	begin() const;

		template <class TElem>
		typename traits<TElem>::const_iterator	end() const;
	///	\}

	///	iterators over subset of surface grid of toplevel
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator	begin(int si);

		template <class TElem>
		typename traits<TElem>::iterator	end(int si);

		template <class TElem>
		typename traits<TElem>::const_iterator	begin(int si) const;

		template <class TElem>
		typename traits<TElem>::const_iterator	end(int si) const;
	/// \}

	///	iterators over whole surface grid w.r.t to level 'lvl'
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator	surface_level_begin(int lvl);

		template <class TElem>
		typename traits<TElem>::iterator	surface_level_end(int lvl);

		template <class TElem>
		typename traits<TElem>::const_iterator	surface_level_begin(int lvl) const;

		template <class TElem>
		typename traits<TElem>::const_iterator	surface_level_end(int lvl) const;
	///	\}

	///	iterators over subset of surface grid w.r.t to level 'lvl
	///	 \{
		template <class TElem>
		typename traits<TElem>::iterator	surface_level_begin(int si, int lvl);

		template <class TElem>
		typename traits<TElem>::iterator	surface_level_end(int si, int lvl);

		template <class TElem>
		typename traits<TElem>::const_iterator	surface_level_begin(int si, int lvl) const;

		template <class TElem>
		typename traits<TElem>::const_iterator	surface_level_end(int si, int lvl) const;
	///	\}

	///	returns the level in grid hierarchy of an element in the surface
		template <class TGeomObj>
		inline int get_level(TGeomObj* obj) const;

	///	returns if the element is contained in the surface view
	/**	Retruns true e.g. for unshadowed constrained (hanging) vertices
	 * \sa SurfaceView::is_shadowed*/
		template <class TGeomObj>
		inline bool is_contained(TGeomObj* obj) const;

	///	returns if the element is shadowed and thus not contained in the surface view
	/**	A shadowed element has at least one child, which is also a member of the
	 * surface view.
	 * \sa SurfaceView::is_contained*/
		template <class TGeomObj>
		inline bool is_shadowed(TGeomObj* obj) const;

	///	returns if the element is ghost and thus not contained in the surface view
	/**	Vertical master elements which are not contained in a horizontal interface
	 * are called ghost elements, since no operations should be executed on them.
	 * Those operations are instead executed on associated vertical slave elements.
	 * \sa SurfaceView::is_contained*/
		template <class TGeomObj>
		inline bool is_ghost(TGeomObj* obj) const;

	///	returns if an element has children
		template <typename TBaseElem>
		inline bool has_children(TBaseElem* elem) const;

	///	returns parent != NULL if copy
		template <typename TBaseElem>
		inline TBaseElem* parent_if_copy(TBaseElem* elem) const;

	///	returns parent != NULL if of same base object type
		template <typename TBaseElem>
		inline TBaseElem* parent_if_same_type(TBaseElem* elem) const;

	///	returns child != NULL if copy
		template <typename TBaseElem>
		inline TBaseElem* child_if_copy(TBaseElem* elem) const;

	///	marks the shadows, must be called after a grid change
		void mark_shadows();

	public:
	///	Iterator to traverse the surface of a multi-grid hierarchy
		template <class TElem>
		class SurfaceViewElementIterator
//			: public boost::iterator_facade<SurfaceViewElementIterator<TElem>, TElem*,
//											boost::forward_traversal_tag>
		{
			public:
				SurfaceViewElementIterator();

			private:
				typedef SurfaceViewElementIterator<TElem> this_type;
				typedef TElem* TValue;

				friend class SurfaceView;
				friend class ConstSurfaceViewElementIterator<TElem>;

				SurfaceViewElementIterator(SurfaceView* surfView,
				                           int fromSubset, int toSubset,
				                           int startLvl, int topLvl,
				                           typename geometry_traits<TElem>::iterator elemIter);

			public:
				this_type operator ++()	{increment(); return *this;}
				this_type operator ++(int unused)	{this_type i = *this; increment(); return i;}

				bool operator ==(const this_type& iter) const {return equal(iter);}
				bool operator !=(const this_type& iter) const {return !equal(iter);}

				TValue operator *()	{return dereference();}

			private:
				friend class boost::iterator_core_access;

				inline bool equal(SurfaceViewElementIterator<TElem> const& other) const;

			///	returns next valid iterator
				void increment();

			///	returns begin-iterator of next non-empty section, returns false if not available
				bool increment_section();

			///	dereference
				inline TValue dereference() const;

			private:
			// \todo: Maybe it would be convenient to drop some of those variabes
			//		  below to save copy-time when an iterator is copied.
				int m_fromSI;
				int m_toSI;
				int m_si;
				int m_lvl;
				int m_topLvl;
				typename geometry_traits<TElem>::iterator m_elemIter;
				typename geometry_traits<TElem>::iterator m_iterEndSection;
				SurfaceView* m_surfView;
		};

	///	Const iterator to traverse the surface of a multi-grid hierarchy
		template <class TElem>
		class ConstSurfaceViewElementIterator
//			: public boost::iterator_facade<ConstSurfaceViewElementIterator<TElem>, /*const*/ TElem*,
//											boost::forward_traversal_tag>
		{
			public:
				ConstSurfaceViewElementIterator();

				ConstSurfaceViewElementIterator(const SurfaceViewElementIterator<TElem>& iter);

			private:
				typedef ConstSurfaceViewElementIterator<TElem> this_type;
				// \todo: should return const TElem*
				typedef TElem* TValue;
				//typedef const TElem* TValue;

				friend class SurfaceView;

				ConstSurfaceViewElementIterator(const SurfaceView* surfView,
				                                int fromSubset, int toSubset,
				                                int startLvl, int topLvl,
				                                typename geometry_traits<TElem>::const_iterator elemIter);

			public:
				this_type operator ++()	{increment(); return *this;}
				this_type operator ++(int unused)	{this_type i = *this; increment(); return i;}

				bool operator ==(const this_type& iter) const {return equal(iter);}
				bool operator !=(const this_type& iter) const {return !equal(iter);}

				TValue operator *()	{return dereference();}

			private:
				friend class boost::iterator_core_access;

				inline bool equal(ConstSurfaceViewElementIterator<TElem> const& other) const;

			///	returns next valid iterator
				void increment();

			///	returns begin-iterator of next non-empty section, returns false if not available
				bool increment_section();

			///	dereference
				inline TValue dereference() const;

			private:
				int m_fromSI;
				int m_toSI;
				int m_si;
				int m_lvl;
				int m_topLvl;
				typename geometry_traits<TElem>::const_iterator m_elemIter;
				typename geometry_traits<TElem>::const_iterator m_iterEndSection;
				const SurfaceView* m_surfView;
		};

	private:
		SmartPtr<MGSubsetHandler> 	m_spMGSH;
		bool						m_adaptiveMG;
		MultiGrid*					m_pMG;
		BoolMarker					m_Marker;
#ifdef UG_PARALLEL
		DistributedGridManager* 	m_pDistGridMgr;
#endif
};


///	Represents a surface-level view of a given multi-grid hierarchy
/**	The SurfaceLevelView allows to iterate over all subsets of a given
 * surface level. A surface level therby consists exactly of all surface elements
 * up to the specified level and of all the elements in the specified level.
 *
 * A SurfaceLevelView is always constructed for a given SurfaceView. If you
 * pass the constant SLV_TOPLEVEL instead of a concrete level during construction,
 * then the SurfaceLevelView will always represent the topmost surface-level of
 * the associated multi-grid hierarchy.
 */
class SurfaceLevelView
{
	public:
		template <class TElem>
		struct traits{
			typedef SurfaceView::SurfaceViewElementIterator<TElem>		iterator;
			typedef SurfaceView::ConstSurfaceViewElementIterator<TElem>	const_iterator;
		};

		enum{SLV_TOPLEVEL = -1};

	public:
		SurfaceLevelView(SmartPtr<SurfaceView> spSV, int topLvl = SLV_TOPLEVEL);

	///	number of subsets
		int num_subsets() const {return m_spSV->num_subsets();}

	///	retruns surface view
		ConstSmartPtr<SurfaceView> surface_view() const {return m_spSV;}

	///	returns the level in grid hierarchy of an element in the surface
		template <class TGeomObj>
		inline int get_level(TGeomObj* obj) const {return m_spSV->get_level(obj);}

	///	returns the adjacend elements
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
		                        TElem* elem, bool clearContainer = true) const;

	///	iterator over whole grid
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator	begin();

		template <class TElem>
		typename traits<TElem>::iterator	end();

		template <class TElem>
		typename traits<TElem>::const_iterator	begin() const;

		template <class TElem>
		typename traits<TElem>::const_iterator	end() const;
	///	\}

	///	iterator over subset of grid
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator	begin(int si);

		template <class TElem>
		typename traits<TElem>::iterator	end(int si);

		template <class TElem>
		typename traits<TElem>::const_iterator	begin(int si) const;

		template <class TElem>
		typename traits<TElem>::const_iterator	end(int si) const;
	/// \}

	private:
		SmartPtr<SurfaceView>	m_spSV;
		int						m_topLvl;
};

/** \} */

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "surface_view_impl.hpp"

#endif
