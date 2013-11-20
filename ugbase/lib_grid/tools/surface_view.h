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
#include "lib_grid/tools/grid_level.h"
#include "subset_handler_multi_grid.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "common/util/flags.h"

namespace ug{

/** \ingroup lib_grid_tools
 *  \{ */

///	Represents the surface view of a multi-grid hierarchy.
/**	The surface of a multi-grid hierarchy consists of all elements, which
 * do not have children.
 * The surface view allows to iterate over the subsets of such a surface.
 * The surface_begin(...) and surface_end(...) methods even allow
 * to ignore all levels above a specified level, resulting in a new surface-view
 * in which only surface elements up to the specified level and the elements in
 * the specified level are regarded as surface-view-elements.
 */
class SurfaceView
{
	public:
		enum{TOPLEVEL = -1};

		template <class TElem> class SurfaceViewElementIterator;
		template <class TElem> class ConstSurfaceViewElementIterator;

		template <class TElem>
		struct traits{
			typedef SurfaceViewElementIterator<TElem>		iterator;
			typedef ConstSurfaceViewElementIterator<TElem>	const_iterator;
		};

		enum SurfaceConstants{
			UNASSIGNED = 1,
			PURE_SURFACE = 1 << 1,
			SHADOWING = 1 << 2,
			SHADOW_COPY = 1 << 3,
			SHADOW_NONCOPY = 1 << 4,

			SHADOW = SHADOW_COPY | SHADOW_NONCOPY,
			ALL = PURE_SURFACE | SHADOWING | SHADOW,
			SURFACE_AND_SHADOWING = PURE_SURFACE | SHADOWING,
			SURFACE_NONCOPY = PURE_SURFACE | SHADOWING | SHADOW_NONCOPY
		};
		typedef Flag<SurfaceConstants, byte, SS_NONE>	SurfaceState;
		typedef Attachment<SurfaceState>				ASurfaceState;

	public:
		SurfaceView(SmartPtr<MGSubsetHandler> spMGSH,
		            bool adaptiveMG = true);

		~SurfaceView();

	///	returns underlying subset handler
		inline SmartPtr<MGSubsetHandler> subset_handler();

	///	returns underlying subset handler
		inline ConstSmartPtr<MGSubsetHandler> subset_handler() const;

	///	returns if multigrid is adaptive
		inline bool is_adaptive() const;

	///	number of subsets
		int num_subsets() const {return m_spMGSH->num_subsets();}

	///	iterators of grid level
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator
		begin(int si, const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING);

		template <class TElem>
		typename traits<TElem>::iterator
		end(int si, const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING);

		template <class TElem>
		typename traits<TElem>::const_iterator
		begin(int si, const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING) const;

		template <class TElem>
		typename traits<TElem>::const_iterator
		end(int si, const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING) const;

		template <class TElem>
		typename traits<TElem>::iterator
		begin(const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING);

		template <class TElem>
		typename traits<TElem>::iterator
		end(const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING);

		template <class TElem>
		typename traits<TElem>::const_iterator
		begin(const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING) const;

		template <class TElem>
		typename traits<TElem>::const_iterator
		end(const GridLevel& gl, SurfaceConstants validSurfStates = SURFACE_AND_SHADOWING) const;
	///	\}

	public:
	///	returns the level in grid hierarchy of an element in the surface
		template <class TGeomObj>
		inline int get_level(TGeomObj* obj) const;

	///	returns if the element is contained in the surface view
	/**	Retruns true e.g. for unshadowed constrained (hanging) vertices.
	 * A top-level may optionally be specified. True is returned for all elements
	 * in the top level, as long as they aren't ghosts. If you specify a topLevel < 0,
	 * the topLevel parameter will be ignored and the method behaves as if no topLevel
	 * was specified.
	 * \sa SurfaceView::is_shadowed, SurfaceView::is_shadowing
	 * \{ */
		template <class TGeomObj>
		inline bool is_surface_element(TGeomObj* obj) const;

		template <class TGeomObj>
		inline bool is_surface_element(TGeomObj* obj, int topLevel) const;
	/** \} */

	///	returns if the element is ghost
	/**	ghost elements are vertical masters that are in no other interfaces.*/
		template <class TGeomObj>
		inline bool is_ghost(TGeomObj* obj) const;

	///	returns if the element is shadowed and thus not contained in the surface view
	/**	A shadowed element has a child and at least one adjacent element which is
	 * a surface element
	 * \sa SurfaceView::is_surface_element, SurfaceView::is_shadowing*/
		template <class TGeomObj>
		inline bool is_shadowed(TGeomObj* obj) const;

	///	returns true if the given element is a shadowing element
	/**	An element is considered to be shadowing if it is the child of a
	 * shadowed element. Not that in a distributed grid, this can be true even
	 * if the element has no parent on the local process.
	 * \sa SurfaceView::is_shadowed, SurfaceView::is_surface_element*/
		template <class TGeomObj>
		inline bool is_shadowing(TGeomObj* obj) const;

	///	returns parent != NULL if copy
		template <typename TBaseElem>
		inline TBaseElem* parent_if_copy(TBaseElem* elem) const;

	///	returns parent != NULL if of same base object type
		template <typename TBaseElem>
		inline TBaseElem* parent_if_same_type(TBaseElem* elem) const;

	///	returns child != NULL if copy
		template <typename TBaseElem>
		inline TBaseElem* child_if_copy(TBaseElem* elem) const;

	///	refresh_surface_states must be called after a grid change
		void refresh_surface_states();

	///	returns an or combination of current surface states
	/**	Please use the methods is_surface_element, is_shadowed and is_shadowing
	 * instead of this method.
	 * The only reason it is publicly accessible is for debugging reasons.
	 *
	 * \note a protected non-const version exists, which returns a reference to the state.*/
		template <class TElem>
		SurfaceState get_surface_state(TElem* elem) const		{return m_aaSurfState[elem];}

	///	returns the adjacend elements w.r.t. the surface view
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
								TElem* elem, bool clearContainer = true) const;

	public:
	///	Iterator to traverse the surface of a multi-grid hierarchy
		template <class TElem>
		class SurfaceViewElementIterator
		{
			public:
				SurfaceViewElementIterator();

			private:
				typedef SurfaceViewElementIterator<TElem> this_type;
				typedef TElem* TValue;

				friend class SurfaceView;
				friend class ConstSurfaceViewElementIterator<TElem>;

				SurfaceViewElementIterator(bool start,
				                           SurfaceView* surfView,
				                           const GridLevel& gl,
				                           SurfaceState validStates,
				                           int si = -1);

			public:
				this_type operator ++()	{increment(); return *this;}
				this_type operator ++(int unused)	{this_type i = *this; increment(); return i;}

				bool operator ==(const this_type& iter) const {return equal(iter);}
				bool operator !=(const this_type& iter) const {return !equal(iter);}

				TValue operator *()	{return dereference();}

			private:
				friend class boost::iterator_core_access;

				inline bool equal(SurfaceViewElementIterator<TElem> const& other) const;

			///	returns if valid element, i.e. contained in iterator loop
				template <class TGeomObj>
				inline bool is_contained(TGeomObj* obj) const;

			///	returns next valid iterator
				void increment();

			///	returns begin-iterator of next non-empty section, returns false if not available
				bool increment_section();

			///	dereference
				inline TValue dereference() const;

			private:
				SurfaceView* m_pSurfView;
				GridLevel m_gl;
				SurfaceState m_validSurfStates;
				int m_fromSI;
				int m_toSI;
				int m_si;
				int m_topLvl;
				int m_lvl;
				typename geometry_traits<TElem>::iterator m_elemIter;
				typename geometry_traits<TElem>::iterator m_iterEndSection;
		};

	///	Const iterator to traverse the surface of a multi-grid hierarchy
		template <class TElem>
		class ConstSurfaceViewElementIterator
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

				ConstSurfaceViewElementIterator(bool start,
						                        const SurfaceView* surfView,
						                        const GridLevel& gl,
						                        SurfaceState validStates,
						                        int si = -1);

			public:
				this_type operator ++()	{increment(); return *this;}
				this_type operator ++(int unused)	{this_type i = *this; increment(); return i;}

				bool operator ==(const this_type& iter) const {return equal(iter);}
				bool operator !=(const this_type& iter) const {return !equal(iter);}

				TValue operator *()	{return dereference();}

			private:
				friend class boost::iterator_core_access;

				inline bool equal(ConstSurfaceViewElementIterator<TElem> const& other) const;

			///	returns if valid element, i.e. contained in iterator loop
				template <class TGeomObj>
				inline bool is_contained(TGeomObj* obj) const;

			///	returns next valid iterator
				void increment();

			///	returns begin-iterator of next non-empty section, returns false if not available
				bool increment_section();

			///	dereference
				inline TValue dereference() const;

			private:
				const SurfaceView* m_pSurfView;
				GridLevel m_gl;
				SurfaceState m_validSurfStates;
				int m_fromSI;
				int m_toSI;
				int m_si;
				int m_topLvl;
				int m_lvl;
				typename geometry_traits<TElem>::const_iterator m_elemIter;
				typename geometry_traits<TElem>::const_iterator m_iterEndSection;
			};

	private:
	///	returns true if the element is a surface element locally
	/**	This method disregards possible copies of the given element on other processes.
	 * The method is used during surface-state computation in refresh_surface_states.
	 * The method returns false if the given element has children or if it is
	 * contained in a vertical master interface.*/
		template <class TElem>
		bool is_local_surface_view_element(TElem* elem);

	///	only call for elements of highest dimension
		template <class TElem>
		void refresh_surface_states();

	///	recursively marks sides and sides of sides as surface or shadow
	/**	This method is used during refresh_surface_states to assign surface
	 * states to sides of surface elements.
	 * Make sure that all elements in lower levels have already been processed!*/
		template <class TElem, class TSide>
		void mark_sides_as_surface_or_shadow(TElem* elem,
											 byte surfaceState = PURE_SURFACE);

		template <class TElem>
		void mark_shadowing(bool markSides = false);

	///	adjusts surface states in a parallel environment
		template <class TElem>
		void adjust_parallel_surface_states();

		template <class TElem>
		SurfaceState surface_state(TElem* elem) const	{return m_aaSurfState[elem];}

		template <class TElem>
		SurfaceState& surface_state(TElem* elem)		{return m_aaSurfState[elem];}

//		template <class TElem>
//		void set_surface_state(TElem* elem, byte ss)	{m_aaElemSurfState[elem] = ss;}
//
//		template <class TElem>
//		bool has_surface_state(TElem* elem, byte s) const	{return (surface_state(elem) & s) == s;}

		template <class TElem>
		bool is_vmaster(TElem* elem) const;

	private:
		SmartPtr<MGSubsetHandler> 		m_spMGSH;
		bool							m_adaptiveMG;
		MultiGrid*						m_pMG;
		DistributedGridManager*			m_distGridMgr;
		ASurfaceState									m_aSurfState;
		MultiElementAttachmentAccessor<ASurfaceState>	m_aaSurfState;
};

/** \} */

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "surface_view_impl.hpp"

#endif
