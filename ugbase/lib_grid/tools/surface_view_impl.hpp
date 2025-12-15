/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW_IMPL__
#define __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW_IMPL__

#include "surface_view.h"

#include "lib_grid/parallelization/distributed_grid.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
//	implementation of SurfaceViewElementIterator
////////////////////////////////////////////////////////////////////////////////

template <typename TElem>
SurfaceView::SurfaceViewElementIterator<TElem>::
SurfaceViewElementIterator(bool start,
                           SurfaceView* sv,
                           const GridLevel& gl,
                           SurfaceState validStates,
                           int si) :
	m_pSurfView(sv),
	m_gl(gl),
	m_validStates(validStates),

	m_fromSI( (si >= 0) ? si : 0 ),
	m_toSI( (si >= 0) ? si : (sv->subset_handler()->num_subsets() - 1) ),
	m_si(start ? m_fromSI : m_toSI),

	m_topLvl(gl.top() ? (sv->subset_handler()->num_levels()-1) : gl.level()),
	m_lvl((start && gl.is_surface() && sv->is_adaptive()) ? 0 : m_topLvl),

	m_elemIter(start ? sv->subset_handler()->begin<TElem>(m_si, m_lvl)
					 : sv->subset_handler()->end<TElem>(m_toSI, m_topLvl)),
	m_iterEndSection(sv->subset_handler()->end<TElem>(m_si, m_lvl))

{
	UG_ASSERT(m_topLvl >= 0 && m_topLvl < (int)sv->subset_handler()->num_levels(),
	          "Invalid level: "<<m_topLvl<<" [min: 0, max: "<<sv->subset_handler()->num_levels()<<"]");
	UG_ASSERT(m_lvl >= 0 && m_lvl < (int)sv->subset_handler()->num_levels(),
	          "Invalid level: "<<m_lvl<<" [min: 0, max: "<<sv->subset_handler()->num_levels()<<"]");
	UG_ASSERT(m_si >= 0 && m_si < sv->subset_handler()->num_subsets(),
	          "Invalid subset: "<<m_si<<" [min: 0, max: "<<sv->subset_handler()->num_subsets()<<"]");
	UG_ASSERT(m_toSI >= 0 && m_toSI < sv->subset_handler()->num_subsets(),
	          "Invalid subset: "<<m_toSI<<" [min: 0, max: "<<sv->subset_handler()->num_subsets()<<"]");

//	if at end of section -> increase until next non-empty section
	if(m_elemIter == m_iterEndSection)
		if(!increment_section())
			return;

//	m_elemIter has to point to a valid surface view element
	if(!is_contained(*m_elemIter)){increment(); return;}
}

template <typename TElem>
SurfaceView::SurfaceViewElementIterator<TElem>::
SurfaceViewElementIterator() :
	m_pSurfView(nullptr),
	m_validStates(SurfaceConstants::MG_UNDEFINED),
	m_fromSI(0),
	m_toSI(0),
	m_si(0),
	m_topLvl(0),
	m_lvl(0),
	m_elemIter(),
	m_iterEndSection()
{}

template <typename TElem>
bool SurfaceView::SurfaceViewElementIterator<TElem>::
equal(SurfaceViewElementIterator const& other) const
{
	return (m_elemIter == other.m_elemIter);
}

template <typename TElem>
bool SurfaceView::SurfaceViewElementIterator<TElem>::
increment_section()
{
//	check if end of section reached
	do
	{
	//	a) if still subsets left to loop on level
		if(m_si < m_toSI)
		{
		//	increase subset, set new section iterators
			++m_si;
			m_elemIter = m_pSurfView->subset_handler()->begin<TElem>(m_si, m_lvl);
			m_iterEndSection = m_pSurfView->subset_handler()->end<TElem>(m_si, m_lvl);
		}
	//	b) if still levels left to be looped
		else if(m_lvl < m_topLvl)
		{
		//	increase level, reset subset to fromSubset, set new section iterators
			++m_lvl;
			m_si = m_fromSI;
			m_elemIter = m_pSurfView->subset_handler()->begin<TElem>(m_si, m_lvl);
			m_iterEndSection = m_pSurfView->subset_handler()->end<TElem>(m_si, m_lvl);
		}
	//	c) no section left, we're done (m_elemIter is end iterator now)
		else {
			return false;
		}
	}
	while(m_elemIter == m_iterEndSection);
	return true;
}

template <typename TElem>
void SurfaceView::SurfaceViewElementIterator<TElem>::
increment()
{
//	we search the next non-shadowed element
	do
	{
	//	increase iterator
		++m_elemIter;

	//	check if end of section reached
		if(m_elemIter == m_iterEndSection){
			if(!increment_section())
				return;
		}

	}while(!is_contained(*m_elemIter));
}

template <typename TElem>
typename SurfaceView::SurfaceViewElementIterator<TElem>::TValue
SurfaceView::SurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}

template <typename TElem>
template <typename TGeomObj>
bool SurfaceView::SurfaceViewElementIterator<TElem>::
is_contained(TGeomObj* obj) const
{
	#ifdef UG_PARALLEL
	if(m_pSurfView->is_ghost(obj)){
		if(m_gl.ghosts() && (m_lvl == m_topLvl)) return true;
		else return false;
	}
	#endif

	SurfaceState oss = m_pSurfView->surface_state(obj);

	if( m_validStates.contains(SurfaceConstants::TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE)
		&& (m_lvl == m_topLvl)
		&& oss.partially_contains(SurfaceConstants::MG_SHADOW))
	{
		oss = SurfaceConstants::MG_SURFACE_PURE;
	}

	return m_validStates.contains(oss);
}


////////////////////////////////////////////////////////////////////////////////
//	implementation of ConstSurfaceViewElementIterator
////////////////////////////////////////////////////////////////////////////////

template <typename TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator(const SurfaceViewElementIterator<TElem>& iter)
{
	m_pSurfView = iter.m_pSurfView;
	m_gl = iter.m_gl;
	m_validStates = iter.m_validStates;
	m_fromSI = iter.m_fromSI;
	m_toSI = iter.m_toSI;
	m_si = iter.m_si;
	m_topLvl = iter.m_topLvl;
	m_lvl = iter.m_lvl;
	m_elemIter = typename geometry_traits<TElem>::const_iterator(iter.m_elemIter);
	m_iterEndSection = typename geometry_traits<TElem>::const_iterator(iter.m_iterEndSection);
}

template <typename TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator() :
	m_pSurfView(nullptr),
	m_validStates(SurfaceConstants::MG_UNDEFINED),
	m_fromSI(0),
	m_toSI(0),
	m_si(0),
	m_topLvl(0),
	m_lvl(0),
	m_elemIter(),
	m_iterEndSection()
{}

template <typename TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator(bool start,
                                const SurfaceView* sv,
                                const GridLevel& gl,
                                SurfaceState validStates,
                                int si) :
	m_pSurfView(sv),
	m_gl(gl),
	m_validStates(validStates),

	m_fromSI( (si >= 0) ? si : 0 ),
	m_toSI( (si >= 0) ? si : (sv->subset_handler()->num_subsets() - 1) ),
	m_si(start ? m_fromSI : m_toSI),

	m_topLvl(gl.top() ? (sv->subset_handler()->num_levels()-1) : gl.level()),
	m_lvl((start && gl.is_surface() && sv->is_adaptive()) ? 0 : m_topLvl),

	m_elemIter(start ? sv->subset_handler()->begin<TElem>(m_si, m_lvl)
					 : sv->subset_handler()->end<TElem>(m_toSI, m_topLvl)),
	m_iterEndSection(sv->subset_handler()->end<TElem>(m_si, m_lvl))

{
	UG_ASSERT(m_topLvl >= 0 && m_topLvl < (int)sv->subset_handler()->num_levels(),
			  "Invalid level: "<<m_topLvl<<" [min: 0, max: "<<sv->subset_handler()->num_levels()<<"]");
	UG_ASSERT(m_lvl >= 0 && m_lvl < (int)sv->subset_handler()->num_levels(),
			  "Invalid level: "<<m_lvl<<" [min: 0, max: "<<sv->subset_handler()->num_levels()<<"]");
	UG_ASSERT(m_si >= 0 && m_si < sv->subset_handler()->num_subsets(),
			  "Invalid subset: "<<m_si<<" [min: 0, max: "<<sv->subset_handler()->num_subsets()<<"]");
	UG_ASSERT(m_toSI >= 0 && m_toSI < sv->subset_handler()->num_subsets(),
			  "Invalid subset: "<<m_toSI<<" [min: 0, max: "<<sv->subset_handler()->num_subsets()<<"]");

//	if at end of section -> increase until next non-empty section
	if(m_elemIter == m_iterEndSection)
		if(!increment_section())
			return;

//	m_elemIter has to point to a valid surface view element
	if(!is_contained(*m_elemIter)){increment(); return;}
}

template <typename TElem>
bool SurfaceView::ConstSurfaceViewElementIterator<TElem>::
equal(ConstSurfaceViewElementIterator const& other) const
{
	return (m_elemIter == other.m_elemIter);
}

template <typename TElem>
bool SurfaceView::ConstSurfaceViewElementIterator<TElem>::
increment_section()
{
//	check if end of section reached
	do
	{
	//	a) if still subsets left to loop on level
		if(m_si < m_toSI)
		{
		//	increase subset, set new section iterators
			++m_si;
			m_elemIter = m_pSurfView->subset_handler()->begin<TElem>(m_si, m_lvl);
			m_iterEndSection = m_pSurfView->subset_handler()->end<TElem>(m_si, m_lvl);
		}
	//	b) if still levels left to be looped
		else if(m_lvl < m_topLvl)
		{
		//	increase level, reset subset to fromSubset, set new section iterators
			++m_lvl;
			m_si = m_fromSI;
			m_elemIter = m_pSurfView->subset_handler()->begin<TElem>(m_si, m_lvl);
			m_iterEndSection = m_pSurfView->subset_handler()->end<TElem>(m_si, m_lvl);
		}
	//	c) no section left, we're done (m_elemIter is end iterator now)
		else {
			return false;
		}
	}
	while(m_elemIter == m_iterEndSection);
	return true;
}

template <typename TElem>
void SurfaceView::ConstSurfaceViewElementIterator<TElem>::
increment()
{
//	this iterator should only be enabled for optimization work
	// PROFILE_FUNC_GROUP("SurfaceView::iterator")
//	we search the next non-shadowed element
	do
	{
	//	increase iterator
		++m_elemIter;

	//	check if end of section reached
		if(m_elemIter == m_iterEndSection){
			if(!increment_section())
				return;
		}
	}while(!is_contained(*m_elemIter));
}

template <typename TElem>
typename SurfaceView::ConstSurfaceViewElementIterator<TElem>::TValue
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}

template <typename TElem>
template <typename TGeomObj>
bool SurfaceView::ConstSurfaceViewElementIterator<TElem>::
is_contained(TGeomObj* obj) const
{
	#ifdef UG_PARALLEL
	if(m_pSurfView->is_ghost(obj)){
		if(m_gl.ghosts() && (m_lvl == m_topLvl)) return true;
		else return false;
	}
	#endif

	SurfaceState oss = m_pSurfView->surface_state(obj);

	if( m_validStates.contains(SurfaceConstants::TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE)
		&& (m_lvl == m_topLvl)
		&& oss.partially_contains(SurfaceConstants::MG_SHADOW))
	{
		oss = SurfaceConstants::MG_SURFACE_PURE;
	}

	return m_validStates.contains(oss);
}

////////////////////////////////////////////////////////////////////////////////
//	grid level iterators
////////////////////////////////////////////////////////////////////////////////

template <typename TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
begin(int si, const GridLevel& gl, SurfaceState validStates)
{
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::iterator(true, this, gl, validStates, si);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
end(int si, const GridLevel& gl, SurfaceState validStates)
{
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::iterator(false, this, gl, validStates, si);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
begin(int si, const GridLevel& gl, SurfaceState validStates) const
{
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::const_iterator(true, this, gl, validStates, si);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
end(int si, const GridLevel& gl, SurfaceState validStates) const
{
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::const_iterator(false, this, gl, validStates, si);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
begin(const GridLevel& gl, SurfaceState validStates)
{
	return typename traits<TElem>::iterator(true, this, gl, validStates);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
end(const GridLevel& gl, SurfaceState validStates)
{
	return typename traits<TElem>::iterator(false, this, gl, validStates);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
begin(const GridLevel& gl, SurfaceState validStates) const
{
	return typename traits<TElem>::const_iterator(true, this, gl, validStates);
}

template <typename TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
end(const GridLevel& gl, SurfaceState validStates) const
{
	return typename traits<TElem>::const_iterator(false, this, gl, validStates);
}

////////////////////////////////////////////////////////////////////////////////
//	implementation of SurfaceView
////////////////////////////////////////////////////////////////////////////////

SmartPtr<MGSubsetHandler> SurfaceView::subset_handler()
{
	return m_spMGSH;
}

ConstSmartPtr<MGSubsetHandler> SurfaceView::subset_handler() const
{
	return m_spMGSH;
}

bool SurfaceView::is_adaptive() const
{
	return m_adaptiveMG;
}

template <typename TGeomObj>
bool SurfaceView::is_contained(TGeomObj* obj, const GridLevel& gl,
                               SurfaceState validStates) const
{
	const int lvl = m_pMG->get_level(obj);
	const int topLvl = (gl.top() ? (subset_handler()->num_levels()-1) : gl.level());
	if(lvl > topLvl) return false;

	#ifdef UG_PARALLEL
	if(is_ghost(obj)){
		if(gl.ghosts() && (lvl == topLvl)) return true;
		else return false;
	}
	#endif

	SurfaceState oss = surface_state(obj);

	if( validStates.contains(SurfaceConstants::TREAT_TOP_LVL_SHADOWS_AS_SURFACE_PURE)
		&& (lvl == topLvl)
		&& oss.partially_contains(SurfaceConstants::MG_SHADOW))
	{
		oss = SurfaceConstants::MG_SURFACE_PURE;
	}

	return validStates.contains(oss);
}

template <typename TElem>
SurfaceView::SurfaceState SurfaceView::surface_state(TElem* elem, const GridLevel& gl) const
{
	const int lvl = m_pMG->get_level(elem);
	const int topLvl = (gl.top() ? (subset_handler()->num_levels()-1) : gl.level());
	if(lvl > topLvl)
		UG_THROW("SurfaceView::surface_state: Call only on objects contained "
				"in the grid level. (Else result is undefined)");

	#ifdef UG_PARALLEL
	if(is_ghost(elem)){
		if(gl.ghosts() && (lvl == topLvl)) return SurfaceConstants::MG_SURFACE_PURE;
		else {
			UG_THROW("SurfaceView::surface_state: Call only on objects contained "
					"in the grid level. (Else result is undefined)");
		}
	}
	#endif

	SurfaceState oss = surface_state(elem);

	if( (lvl == topLvl)
		&& oss.partially_contains(SurfaceConstants::MG_SHADOW))
	{
		oss = SurfaceConstants::MG_SURFACE_PURE;
	}

	return oss;
}

template <typename TGeomObj>
bool SurfaceView::is_ghost(TGeomObj* obj) const
{
#ifdef UG_PARALLEL
	return m_distGridMgr->is_ghost(obj);
#else
	return false;
#endif
}

template <typename TGeomObj>
bool SurfaceView::is_shadowed(TGeomObj* obj) const
{
	return surface_state(obj).partially_contains(SurfaceConstants::MG_SHADOW_RIM);
}

template <typename TGeomObj>
bool SurfaceView::is_shadowing(TGeomObj* obj) const
{
	return surface_state(obj).contains(SurfaceConstants::MG_SURFACE_RIM);
}

template <typename TElem>
bool SurfaceView::is_vmaster(TElem* elem) const
{
	#ifdef UG_PARALLEL
		return m_distGridMgr->contains_status(elem, ElementStatusTypes::ES_V_MASTER);
	#else
		return false;
	#endif
}

template <typename TElem, typename TBaseElem>
void SurfaceView::
collect_associated(std::vector<TBaseElem*>& vAssElem,
                   TElem* elem,
                   const GridLevel& gl,
                   bool clearContainer) const
{
	if(clearContainer) vAssElem.clear();

//	collect associated on this level
	if(is_contained(elem, gl, SurfaceConstants::SURFACE)){
		std::vector<TBaseElem*> vCoarseElem;
		CollectAssociated(vCoarseElem, *m_pMG, elem, true);
		for(size_t i = 0; i < vCoarseElem.size(); ++i)
		{
			if(!is_contained(vCoarseElem[i], gl, SurfaceConstants::ALL)) continue;
			vAssElem.push_back(vCoarseElem[i]);
		}
	}

//	if at border of a level grid, there may be connections of the "shadow" element
//	to surface elements on the coarser level. These must be taken into account.
	if(is_contained(elem, gl, SurfaceConstants::SURFACE_RIM))
	{
	//	get parent if copy
		GridObject* pParent = m_pMG->get_parent(elem);
		TElem* parent = dynamic_cast<TElem*>(pParent);
		if(parent == nullptr) return;
		if(m_pMG->num_children<TBaseElem>(parent) != 1) return;

	//	Get connected elements
		std::vector<TBaseElem*> vCoarseElem;
		CollectAssociated(vCoarseElem, *m_pMG, parent, true);

		for(size_t i = 0; i < vCoarseElem.size(); ++i)
		{
			if(!is_contained(vCoarseElem[i], gl, SurfaceConstants::SURFACE)) continue;
			vAssElem.push_back(vCoarseElem[i]);
		}
	}

// alternative implementation taking into account possible SHADOW_RIM_COPY elem.
#if 0
	// collect associated on this level
	if (is_contained(elem, gl, ALL))
	{
		std::vector<TBaseElem*> vCoarseElem;
		CollectAssociated(vCoarseElem, *m_pMG, elem, true);
		for (size_t i = 0; i < vCoarseElem.size(); ++i)
			if (is_contained(vCoarseElem[i], gl, ALL))
				vAssElem.push_back(vCoarseElem[i]);
	}

	// if at border of a grid level, there may be connections of a "shadow-rim-copy" element
	// to surface elements on the finer level. These must be taken into account.
	if (is_contained(elem, gl, SHADOW_RIM_COPY))
	{
		if (m_pMG->num_children<TElem>(elem) > 0)
		{
			TElem* child = m_pMG->get_child<TElem>(elem, 0);
			if (is_contained(child, gl, SURFACE_RIM))
			{
				// get connected elements
				std::vector<TBaseElem*> vFineElem;
				CollectAssociated(vFineElem, *m_pMG, child, true);
				for (size_t i = 0; i < vFineElem.size(); ++i)
					if (is_contained(vFineElem[i], gl, ALL))
						vAssElem.push_back(vFineElem[i]);
			}
		}
	}
#endif
}


}//	end of namespace

#endif
