/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW__
#define __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW__

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif

#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/grid_level.h"
#include "subset_handler_multi_grid.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "common/util/flags.h"
#include "common/types.h"

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
		/**
		 * Byte-Constants, that identify the SurfaceState of an grid-object.
		 *
		 * Every grid-object is in exactly one SurfaceState. In addition some
		 * combinations of the states are named for an easier usage.
		 *
		 * IMPORTANT: The order of the byte-flags is currently crucial. Do not
		 * 			  change them. See ComPol_GatherSurfaceStates.
		 */
		enum SurfaceConstants{
			// each grid-object has exactly on of these states (begin)
			MG_SHADOW_PURE = 1,               ///< full-covered (inner)
			MG_SURFACE_PURE = 1 << 1,         ///< surface, i.e., without children (inner)
			MG_SURFACE_RIM = 1 << 2,          ///< surface, i.e., without children (at rim)
			MG_SHADOW_RIM_COPY = 1 << 3,      ///< covered (at rim) with identical child
			MG_SHADOW_RIM_NONCOPY = 1 << 4,   ///< covered (at rim) with non-identical child(ren)
			// each grid-object has exactly on of these states (end)

			// combo-states with flags as in multi-grid (begin)
			MG_SHADOW_RIM = MG_SHADOW_RIM_COPY | MG_SHADOW_RIM_NONCOPY,    //!< all rim-shadows
			MG_SHADOW = MG_SHADOW_RIM | MG_SHADOW_PURE,                    //!< all shadows
			MG_SURFACE = MG_SURFACE_PURE | MG_SURFACE_RIM,                 //!< all surface
			MG_ALL_BUT_SHADOW_COPY = MG_SURFACE | MG_SHADOW_RIM_NONCOPY,   //!< surface + rim-non-copy-shadows
			MG_ALL = MG_SURFACE | MG_SHADOW_RIM,                           //!< all (except pure-shadow)
			// combo-states with flags as in multi-grid (end)

			TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE = 1 << 5,  //! pseudo-state

			// combo-states with flags as in level-view (begin)
			SHADOW_PURE = MG_SHADOW_PURE               | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			SURFACE_PURE = MG_SURFACE_PURE             | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			SURFACE_RIM = MG_SURFACE_RIM               | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			SHADOW_RIM_COPY = MG_SHADOW_RIM_COPY       | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			SHADOW_RIM_NONCOPY = MG_SHADOW_RIM_NONCOPY | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,

			SHADOW_RIM = MG_SHADOW_RIM | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			SHADOW = MG_SHADOW         | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			SURFACE = MG_SURFACE       | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			ALL_BUT_SHADOW_COPY = MG_ALL_BUT_SHADOW_COPY | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE,
			ALL = MG_ALL               | TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE
			// combo-states with flags as in level-view (end)
		};
		typedef Flag<SurfaceConstants, byte_t, SS_NONE>	SurfaceState;
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

	///	returns if the element is contained in the surface view
		template <class TGeomObj>
		inline bool is_contained(TGeomObj* obj, const GridLevel& gl,
		                         SurfaceState validStates = ALL) const;

	///	returns the surface states, when considered as part of grid level
		template <class TElem>
		SurfaceState surface_state(TElem* elem, const GridLevel& gl) const;

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

	///	refresh_surface_states must be called after a grid change
		void refresh_surface_states();

	///	returns an or combination of current surface states
	/**	Please use the methods is_surface_element, is_shadowed and is_shadowing
	 * instead of this method.
	 * The only reason it is publicly accessible is for debugging reasons.
	 *
	 * \note a protected non-const version exists, which returns a reference to the state.*/
		template <class TElem>
		SurfaceState surface_state(TElem* elem) const {return m_aaSurfState[elem];}

	///	returns the adjacent elements w.r.t. the grid level
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
								TElem* elem,
								const GridLevel& gl,
								bool clearContainer = true) const;

	public:
		template <class TElem> class SurfaceViewElementIterator;
		template <class TElem> class ConstSurfaceViewElementIterator;

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
			///	returns if this iterator equals another
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
				SurfaceState m_validStates;
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
			///	returns if this iterator equals another
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
				SurfaceState m_validStates;
				int m_fromSI;
				int m_toSI;
				int m_si;
				int m_topLvl;
				int m_lvl;
				typename geometry_traits<TElem>::const_iterator m_elemIter;
				typename geometry_traits<TElem>::const_iterator m_iterEndSection;
			};

	public:
		template <class TElem>
		struct traits{
			typedef SurfaceViewElementIterator<TElem>		iterator;
			typedef ConstSurfaceViewElementIterator<TElem>	const_iterator;
		};

	///	iterators of grid level
	///	\{
		template <class TElem>
		typename traits<TElem>::iterator
		begin(int si, const GridLevel& gl, SurfaceState validStates);

		template <class TElem>
		typename traits<TElem>::iterator
		end(int si, const GridLevel& gl, SurfaceState validStates);

		template <class TElem>
		typename traits<TElem>::const_iterator
		begin(int si, const GridLevel& gl, SurfaceState validStates) const;

		template <class TElem>
		typename traits<TElem>::const_iterator
		end(int si, const GridLevel& gl, SurfaceState validStates) const;

		template <class TElem>
		typename traits<TElem>::iterator
		begin(const GridLevel& gl, SurfaceState validStates);

		template <class TElem>
		typename traits<TElem>::iterator
		end(const GridLevel& gl, SurfaceState validStates);

		template <class TElem>
		typename traits<TElem>::const_iterator
		begin(const GridLevel& gl, SurfaceState validStates) const;

		template <class TElem>
		typename traits<TElem>::const_iterator
		end(const GridLevel& gl, SurfaceState validStates) const;
	///	\}

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
											 byte_t surfaceState = MG_SURFACE_PURE);

		template <class TElem>
		void mark_shadowing(bool markSides = false);

	///	adjusts surface states in a parallel environment
		template <class TElem>
		void adjust_parallel_surface_states();

		template <class TElem>
		SurfaceState& surface_state(TElem* elem)		{return m_aaSurfState[elem];}

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
