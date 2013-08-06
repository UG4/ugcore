// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.11.2011 (m,d,y)

#ifndef __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW_IMPL__
#define __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW_IMPL__

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	implementation of SurfaceViewElementIterator
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
SurfaceView::SurfaceViewElementIterator<TElem>::
SurfaceViewElementIterator(SurfaceView* surfView,
                           int fromSubset, int toSubset,
                           int startLvl, int topLvl,
                           bool start, bool withGhosts,
						   SurfaceState validStates) :
	m_fromSI(fromSubset),
	m_toSI(toSubset),
	m_si(m_fromSI),
	m_lvl(startLvl),
	m_topLvl(topLvl),
	m_elemIter(start ? surfView->subset_handler()->begin<TElem>(m_si, m_lvl)
					 : surfView->subset_handler()->end<TElem>(m_toSI, m_topLvl)),
	m_iterEndSection(surfView->subset_handler()->end<TElem>(m_si, m_lvl)),
	m_surfView(surfView),
	m_bWithGhosts(withGhosts),
	m_validSurfStates(validStates)
{
//	if at end of section -> increase until next non-empty section
	if(m_elemIter == m_iterEndSection)
		if(!increment_section())
			return;

//	m_elemIter has to point to a valid surface view element
	if(!is_contained(*m_elemIter)){increment(); return;}
}

template <class TElem>
SurfaceView::SurfaceViewElementIterator<TElem>::
SurfaceViewElementIterator() :
	m_fromSI(0),
	m_toSI(0),
	m_si(0),
	m_lvl(0),
	m_topLvl(0),
	m_elemIter(),
	m_iterEndSection(),
	m_surfView(NULL),
	m_validSurfStates(SurfaceView::SURFACE_AND_SHADOWING)
{}

template <class TElem>
bool SurfaceView::SurfaceViewElementIterator<TElem>::
equal(SurfaceView::SurfaceViewElementIterator<TElem> const& other) const
{
	return (m_elemIter == other.m_elemIter);
}

template <class TElem>
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
			MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
			m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
			m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
		}
	//	b) if still levels left to be looped
		else if(m_lvl < m_topLvl)
		{
		//	increase level, reset subset to fromSubset, set new section iterators
			++m_lvl;
			m_si = m_fromSI;
			MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
			m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
			m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
		}
	//	c) no section left, we're done (m_elemIter is end iterator now)
		else {
			return false;
		}
	}
	while(m_elemIter == m_iterEndSection);
	return true;
}

template <class TElem>
void SurfaceView::SurfaceViewElementIterator<TElem>::
increment()
{
//	we search the next non-shadowed element
	do
	{
	//	increase iterator
		++m_elemIter;

	//	check if end of section reached
		while(m_elemIter == m_iterEndSection)
		{
		//	a) if still subsets left to loop on level
			if(m_si < m_toSI)
			{
			//	increase subset, set new section iterators
				++m_si;
				MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
				m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
				m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
			}
		//	b) if still levels left to be looped
			else if(m_lvl < m_topLvl)
			{
			//	increase level, reset subset to fromSubset, set new section iterators
				++m_lvl;
				m_si = m_fromSI;
				MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
				m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
				m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
			}
		//	c) no section left, we're done (m_elemIter is end iterator now)
			else {
				return;
			}
		}
	}while(!is_contained(*m_elemIter));
}

template <class TElem>
typename SurfaceView::SurfaceViewElementIterator<TElem>::TValue
SurfaceView::SurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}

template <class TElem>
template <class TGeomObj>
bool SurfaceView::SurfaceViewElementIterator<TElem>::
is_contained(TGeomObj* obj) const
{
	UG_ASSERT(m_surfView->get_level(obj) == m_lvl, "Wrong level");

//	if(m_surfView->is_surface_element(obj))
//		return true;
	if(m_validSurfStates.contains(m_surfView->surface_state(obj)))
		return true;

	if(m_lvl == m_topLvl){
#ifdef UG_PARALLEL
		if(m_bWithGhosts) return true;
		else return !m_surfView->is_ghost(obj);
#else
		return true;
#endif
	}

	return false;
}


////////////////////////////////////////////////////////////////////////////////
//	implementation of ConstSurfaceViewElementIterator
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator(const SurfaceView::SurfaceViewElementIterator<TElem>& iter)
{
	m_fromSI = iter.m_fromSI;
	m_toSI = iter.m_toSI;
	m_si = iter.m_si;
	m_lvl = iter.m_lvl;
	m_topLvl = iter.m_topLvl;
	m_elemIter = iter.m_elemIter;
	m_iterEndSection = iter.m_iterEndSection;
	m_surfView = iter.m_surfView;
	m_bWithGhosts = iter.m_bWithGhosts;
	m_validSurfStates = iter.m_validSurfStates;
}

template <class TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator() :
	m_fromSI(0),
	m_toSI(0),
	m_si(0),
	m_lvl(0),
	m_topLvl(0),
	m_elemIter(),
	m_iterEndSection(),
	m_surfView(NULL),
	m_bWithGhosts(false),
	m_validSurfStates(SurfaceView::SURFACE_AND_SHADOWING)
{}

template <class TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator(const SurfaceView* surfView,
                                int fromSubset, int toSubset,
                                int startLvl, int topLvl,
                                bool start, bool withGhosts,
                                SurfaceState validStates) :
	m_fromSI(fromSubset),
	m_toSI(toSubset),
	m_si(m_fromSI),
	m_lvl(startLvl),
	m_topLvl(topLvl),
	m_elemIter(start ? surfView->subset_handler()->begin<TElem>(m_si, m_lvl)
					 : surfView->subset_handler()->end<TElem>(m_toSI, m_topLvl)),
	m_iterEndSection(surfView->subset_handler()->end<TElem>(m_si, m_lvl)),
	m_surfView(surfView),
	m_bWithGhosts(withGhosts),
	m_validSurfStates(validStates)
{
//	if at end of section -> increase until next non-empty section
	if(m_elemIter == m_iterEndSection)
		if(!increment_section())
			return;

//	m_elemIter has to point to a valid surface view element
	if(!is_contained(*m_elemIter)){increment(); return;}
}

template <class TElem>
bool SurfaceView::ConstSurfaceViewElementIterator<TElem>::
equal(SurfaceView::ConstSurfaceViewElementIterator<TElem> const& other) const
{
	return (m_elemIter == other.m_elemIter);
}

template <class TElem>
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
			const MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
			m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
			m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
		}
	//	b) if still levels left to be looped
		else if(m_lvl < m_topLvl)
		{
		//	increase level, reset subset to fromSubset, set new section iterators
			++m_lvl;
			m_si = m_fromSI;
			const MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
			m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
			m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
		}
	//	c) no section left, we're done (m_elemIter is end iterator now)
		else {
			return false;
		}
	}
	while(m_elemIter == m_iterEndSection);
	return true;
}

template <class TElem>
void SurfaceView::ConstSurfaceViewElementIterator<TElem>::
increment()
{
//	we search the next non-shadowed element
	do
	{
	//	increase iterator
		++m_elemIter;

	//	check if end of section reached
		while(m_elemIter == m_iterEndSection)
		{
		//	a) if still subsets left to loop on level
			if(m_si < m_toSI)
			{
			//	increase subset, set new section iterators
				++m_si;
				const MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
				m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
				m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
			}
		//	b) if still levels left to be looped
			else if(m_lvl < m_topLvl)
			{
			//	increase level, reset subset to fromSubset, set new section iterators
				++m_lvl;
				m_si = m_fromSI;
				const MGSubsetHandler& mgsh = *(m_surfView->subset_handler());
				m_elemIter = mgsh.begin<TElem>(m_si, m_lvl);
				m_iterEndSection = mgsh.end<TElem>(m_si, m_lvl);
			}
		//	c) no section left, we're done (m_elemIter is end iterator now)
			else {
				return;
			}
		}
	}while(!is_contained(*m_elemIter));
}

template <class TElem>
typename SurfaceView::ConstSurfaceViewElementIterator<TElem>::TValue
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}

template <class TElem>
template <class TGeomObj>
bool SurfaceView::ConstSurfaceViewElementIterator<TElem>::
is_contained(TGeomObj* obj) const
{
	UG_ASSERT(m_surfView->get_level(obj) == m_lvl, "Wrong level");

//	if(m_surfView->is_surface_element(obj))
//		return true;
	if(m_validSurfStates.contains(m_surfView->surface_state(obj)))
		return true;

	if(m_lvl == m_topLvl){
#ifdef UG_PARALLEL
		if(m_bWithGhosts) return true;
		else return !m_surfView->is_ghost(obj);
#else
		return true;
#endif
	}

	return false;
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

////////////////////////////////////////////////////////////////////////////////
//	surface iterators
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_begin(int si, int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	if(is_adaptive())
		return typename traits<TElem>::iterator(this, si, si, 0, lvl, true, withGhosts, validSurfStates);
	else
		return typename traits<TElem>::iterator(this, si, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_end(int si, int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	return typename traits<TElem>::iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_begin(int si, int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	if(is_adaptive())
		return typename traits<TElem>::const_iterator(this, si, si, 0, lvl, true, withGhosts, validSurfStates);
	else
		return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_end(int si, int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_begin(int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	const int si =  m_spMGSH->num_subsets()-1;
	if(is_adaptive())
		return typename traits<TElem>::iterator(this, 0, si, 0, lvl, true, withGhosts, validSurfStates);
	else
		return typename traits<TElem>::iterator(this, 0, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_end(int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_begin(int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	const int si =  m_spMGSH->num_subsets()-1;
	if(is_adaptive())
		return typename traits<TElem>::const_iterator(this, 0, si, 0, lvl, true, withGhosts, validSurfStates);
	else
		return typename traits<TElem>::const_iterator(this, 0, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_end(int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	if(lvl == TOPLEVEL)
		lvl = m_spMGSH->num_levels() - 1;

	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

////////////////////////////////////////////////////////////////////////////////
//	level iterators
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
level_begin(int si, int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::iterator(this, si, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
level_end(int si, int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
level_begin(int si, int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
level_end(int si, int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	UG_ASSERT(si >= 0 && si < m_spMGSH->num_subsets(), "Invalid subset: "<<si);
	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
level_begin(int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::iterator(this, 0, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
level_end(int lvl, bool withGhosts, SurfaceState validSurfStates)
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
level_begin(int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::const_iterator(this, 0, si, lvl, lvl, true, withGhosts, validSurfStates);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
level_end(int lvl, bool withGhosts, SurfaceState validSurfStates) const
{
	UG_ASSERT(lvl >= 0 && lvl < (int)m_spMGSH->num_levels(), "Invalid level: "<<lvl);
	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl, false, withGhosts, validSurfStates);
}

////////////////////////////////////////////////////////////////////////////////
//	grid level iterators
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
begin(int si, const GridLevel& gl, SurfaceState validSurfStates)
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_begin<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_begin<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
end(int si, const GridLevel& gl, SurfaceState validSurfStates)
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_end<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_end<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
begin(int si, const GridLevel& gl, SurfaceState validSurfStates) const
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_begin<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_begin<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
end(int si, const GridLevel& gl, SurfaceState validSurfStates) const
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_end<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_end<TElem>(si, gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
begin(const GridLevel& gl, SurfaceState validSurfStates)
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_begin<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_begin<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
end(const GridLevel& gl, SurfaceState validSurfStates)
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_end<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_end<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
begin(const GridLevel& gl, SurfaceState validSurfStates) const
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_begin<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_begin<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
end(const GridLevel& gl, SurfaceState validSurfStates) const
{
	if(gl.type() == GridLevel::SURFACE){
		return surface_end<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}else if(gl.type() == GridLevel::LEVEL){
		return level_end<TElem>(gl.level(), gl.with_ghosts(), validSurfStates);
	}
	else UG_THROW("Cannot find level type.");
}

////////////////////////////////////////////////////////////////////////////////
//	util fct
////////////////////////////////////////////////////////////////////////////////

template <class TGeomObj>
bool SurfaceView::is_surface_element(TGeomObj* obj) const
{
	SurfaceState surfState = surface_state(obj);
	return (surfState.contains(PURE_SURFACE)
			||	surfState.contains(SHADOWING))
			&& (!is_shadowed(obj));
}

template <class TGeomObj>
bool SurfaceView::is_surface_element(TGeomObj* obj, int topLevel) const
{
	int lvl = get_level(obj);
	if((topLevel > -1) && (lvl > topLevel))
		return false;

	if(is_surface_element(obj))
		return true;

	if(lvl == topLevel){
		#ifdef UG_PARALLEL
			return !m_distGridMgr->is_ghost(obj);
		#else
			return true;
		#endif
	}
}

template <class TGeomObj>
bool SurfaceView::is_ghost(TGeomObj* obj) const
{
#ifdef UG_PARALLEL
	return m_distGridMgr->is_ghost(obj);
#else
	return false;
#endif
}

template <class TGeomObj>
bool SurfaceView::is_shadowed(TGeomObj* obj) const
{
	return (surface_state(obj) & SHADOW) != 0;
}

template <class TGeomObj>
bool SurfaceView::is_shadowing(TGeomObj* obj) const
{
	return surface_state(obj).contains(SHADOWING);
}

template <class TGeomObj>
int SurfaceView::get_level(TGeomObj* obj) const
{
	return m_pMG->get_level(obj);
}

template <class TElem>
bool SurfaceView::is_vmaster(TElem* elem) const
{
	#ifdef UG_PARALLEL
		return m_distGridMgr->contains_status(elem, ES_V_MASTER);
	#else
		return false;
	#endif
}

template <typename TBaseElem>
TBaseElem* SurfaceView::parent_if_copy(TBaseElem* elem) const
{
	GeometricObject* pParent = m_pMG->get_parent(elem);
	TBaseElem* parent = dynamic_cast<TBaseElem*>(pParent);
	if(parent != NULL &&
		m_pMG->num_children<TBaseElem>(parent) == 1) return parent;
	else return NULL;
}

template <typename TBaseElem>
TBaseElem* SurfaceView::parent_if_same_type(TBaseElem* elem) const
{
	GeometricObject* pParent = m_pMG->get_parent(elem);
	return dynamic_cast<TBaseElem*>(pParent);
}

///	returns child != NULL if copy
template <typename TBaseElem>
TBaseElem* SurfaceView::child_if_copy(TBaseElem* elem) const
{
	if(m_pMG->num_children<TBaseElem>(elem) != 1) return NULL;
	return m_pMG->get_child<TBaseElem>(elem, 0);
}

template <typename TElem, typename TBaseElem>
void SurfaceView::
collect_associated(std::vector<TBaseElem*>& vAssElem,
                   TElem* elem, bool clearContainer) const
{
	if(clearContainer) vAssElem.clear();

//	collect associated on this level
	CollectAssociated(vAssElem, *m_pMG, elem, false);

//	if at border of a level grid, there may be connections of the "shadow" element
//	to surface elements on the coarser level. These must be taken into account.

	if(is_adaptive())
	{
	//	get parent
		TElem* parent = parent_if_copy(elem);

	//	Get connected elements
		if(parent != NULL)
		{
			std::vector<TBaseElem*> vCoarseElem;
			CollectAssociated(vCoarseElem, *m_pMG, parent, true);

			for(size_t i = 0; i < vCoarseElem.size(); ++i)
			{
			//	if shadowed, not in surface view of the requested level
			//todo:	it could make sense to pass the level of this SurfaceLevelView
			//		to the is_surface_element method.
				if(!is_surface_element(vCoarseElem[i])) continue;

			//	else this must be added to adjacent elements
				vAssElem.push_back(vCoarseElem[i]);
			}
		}
	}
}


}//	end of namespace

#endif
