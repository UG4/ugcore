// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.11.2011 (m,d,y)

#ifndef __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW_IMPL__
#define __H__UG__LIB_GRID__TOOLS__SURFACE_VIEW_IMPL__

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	implementation of SurfaceViewElementIterator
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class TElem>
SurfaceView::SurfaceViewElementIterator<TElem>::
SurfaceViewElementIterator(SurfaceView* surfView,
                           int fromSubset, int toSubset,
                           int startLvl, int topLvl,
                           typename geometry_traits<TElem>::iterator elemIter) :
	m_fromSI(fromSubset),
	m_toSI(toSubset),
	m_si(m_fromSI),
	m_lvl(startLvl),
	m_topLvl(topLvl),
	m_elemIter(elemIter),
	m_iterEndSection(surfView->subset_handler()->end<TElem>(m_si, m_lvl)),
	m_surfView(surfView)
{
//	if at end of section -> increase until next non-empty section
	if(m_elemIter == m_iterEndSection)
		if(!increment_section())
			return;

//	m_elemIter has to point to a valid surface view element
	if(!m_surfView->is_surface_element(*m_elemIter)){increment(); return;}
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
	m_surfView(NULL)
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

	//	if on top level no shadows can appear
		if(m_lvl == m_topLvl) return;

	}while(!m_surfView->is_surface_element(*m_elemIter));
}

template <class TElem>
typename SurfaceView::SurfaceViewElementIterator<TElem>::TValue
SurfaceView::SurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	implementation of ConstSurfaceViewElementIterator
////////////////////////////////////////////////////////////////////////////////
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
	m_surfView(NULL)
{}

template <class TElem>
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
ConstSurfaceViewElementIterator(const SurfaceView* surfView,
                                int fromSubset, int toSubset,
                                int startLvl, int topLvl,
                                typename geometry_traits<TElem>::const_iterator elemIter) :
	m_fromSI(fromSubset),
	m_toSI(toSubset),
	m_si(m_fromSI),
	m_lvl(startLvl),
	m_topLvl(topLvl),
	m_elemIter(elemIter),
	m_iterEndSection(surfView->subset_handler()->end<TElem>(m_si, m_lvl)),
	m_surfView(surfView)
{
//	if at end of section -> increase until next non-empty section
	if(m_elemIter == m_iterEndSection)
		if(!increment_section())
			return;

//	m_elemIter has to point to a valid surface view element
	if(!m_surfView->is_surface_element(*m_elemIter)){increment(); return;}
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

	//	if on top level no shadows can appear
		if(m_lvl == m_topLvl) return;

	}while(!m_surfView->is_surface_element(*m_elemIter));
}

template <class TElem>
typename SurfaceView::ConstSurfaceViewElementIterator<TElem>::TValue
SurfaceView::ConstSurfaceViewElementIterator<TElem>::
dereference() const
{
	return *m_elemIter;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	implementation of SurfaceView
////////////////////////////////////////////////////////////////////////////////
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

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
begin(int si)
{
	return surface_level_begin<TElem>(si, m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
end(int si)
{
	return surface_level_end<TElem>(si, m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
begin(int si) const
{
	return surface_level_begin<TElem>(si, m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
end(int si) const
{
	return surface_level_end<TElem>(si, m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
begin()
{
	return surface_level_begin<TElem>(m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
end()
{
	return surface_level_end<TElem>(m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
begin() const
{
	return surface_level_begin<TElem>(m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
end() const
{
	return surface_level_end<TElem>(m_spMGSH->num_levels() - 1);
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_level_begin(int si, int lvl)
{
	if(is_adaptive())
		return typename traits<TElem>::iterator(this, si, si, 0, lvl,
										m_spMGSH->begin<TElem>(si, 0));
	else
		return typename traits<TElem>::iterator(this, si, si, lvl, lvl,
											 	 m_spMGSH->begin<TElem>(si, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_level_end(int si, int lvl)
{
	return typename traits<TElem>::iterator(this, si, si, lvl, lvl,
										 	 m_spMGSH->end<TElem>(si, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_level_begin(int si, int lvl) const
{
	if(is_adaptive())
		return typename traits<TElem>::const_iterator(this, si, si, 0, lvl,
											subset_handler()->begin<TElem>(si, 0));
	else
		return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl,
											subset_handler()->begin<TElem>(si, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_level_end(int si, int lvl) const
{
	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl,
											subset_handler()->end<TElem>(si, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_level_begin(int lvl)
{
	if(is_adaptive())
		return typename traits<TElem>::iterator(this, 0, m_spMGSH->num_subsets()-1, 0, lvl,
										m_spMGSH->begin<TElem>(0, 0));
	else
		return typename traits<TElem>::iterator(this, 0, m_spMGSH->num_subsets()-1, lvl, lvl,
											 	 m_spMGSH->begin<TElem>(0, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::iterator SurfaceView::
surface_level_end(int lvl)
{
	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::iterator(this, si, si, lvl, lvl,
										 	 m_spMGSH->end<TElem>(si, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_level_begin(int lvl) const
{
	if(is_adaptive())
		return typename traits<TElem>::const_iterator(this, 0, m_spMGSH->num_subsets()-1, 0, lvl,
											subset_handler()->begin<TElem>(0, 0));
	else
		return typename traits<TElem>::const_iterator(this, 0, m_spMGSH->num_subsets()-1, lvl, lvl,
											subset_handler()->begin<TElem>(0, lvl));
}

template <class TElem>
typename SurfaceView::traits<TElem>::const_iterator SurfaceView::
surface_level_end(int lvl) const
{
	const int si =  m_spMGSH->num_subsets()-1;
	return typename traits<TElem>::const_iterator(this, si, si, lvl, lvl,
											subset_handler()->end<TElem>(si, lvl));
}


template <class TGeomObj>
bool SurfaceView::is_surface_element(TGeomObj* obj) const
{
	byte surfState = surface_state(obj);
	return (surfState & ESS_SURFACE) && !(surfState & ESS_HIDDEN);
}

template <class TGeomObj>
bool SurfaceView::is_shadowed(TGeomObj* obj) const
{
	return has_surface_state(obj, ESS_SHADOW);
}

template <class TGeomObj>
bool SurfaceView::is_shadowing(TGeomObj* obj) const
{
	return has_surface_state(obj, ESS_SHADOWING);
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	implementation of SurfaceLevelView
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class TElem>
typename SurfaceLevelView::traits<TElem>::iterator
SurfaceLevelView::
begin()
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_begin<TElem>(m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_begin<TElem>(m_topLvl);
}

template <class TElem>
typename SurfaceLevelView::traits<TElem>::iterator
SurfaceLevelView::
end()
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_end<TElem>(m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_end<TElem>(m_topLvl);
}

template <class TElem>
typename SurfaceLevelView::traits<TElem>::const_iterator
SurfaceLevelView::
begin() const
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_begin<TElem>(m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_begin<TElem>(m_topLvl);
}

template <class TElem>
typename SurfaceLevelView::traits<TElem>::const_iterator
SurfaceLevelView::
end() const
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_end<TElem>(m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_end<TElem>(m_topLvl);
}


template <class TElem>
typename SurfaceLevelView::traits<TElem>::iterator
SurfaceLevelView::
begin(int si)
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_begin<TElem>(si,
										m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_begin<TElem>(si, m_topLvl);
}

template <class TElem>
typename SurfaceLevelView::traits<TElem>::iterator
SurfaceLevelView::
end(int si)
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_end<TElem>(si,
										m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_end<TElem>(si, m_topLvl);
}

template <class TElem>
typename SurfaceLevelView::traits<TElem>::const_iterator
SurfaceLevelView::
begin(int si) const
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_begin<TElem>(si,
										m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_begin<TElem>(si, m_topLvl);
}

template <class TElem>
typename SurfaceLevelView::traits<TElem>::const_iterator
SurfaceLevelView::
end(int si) const
{
	if(m_topLvl == SLV_TOPLEVEL)
		return m_spSV->surface_level_end<TElem>(si,
										m_spSV->subset_handler()->num_levels() - 1);
	else
		return m_spSV->surface_level_end<TElem>(si, m_topLvl);
}


template <typename TElem, typename TBaseElem>
void SurfaceLevelView::
collect_associated(std::vector<TBaseElem*>& vAssElem,
                   TElem* elem, bool clearContainer) const
{
	if(clearContainer) vAssElem.clear();

	MultiGrid& multiGrid = *const_cast<MultiGrid*>(m_spSV->subset_handler()->multi_grid());

//	collect associated on this level
	CollectAssociated(vAssElem, multiGrid, elem, false);

//	if at border of a level grid, there may be connections of the "shadow" element
//	to surface elements on the coarser level. These must be taken into account.

	if(m_spSV->is_adaptive())
	{
	//	get parent
		TElem* parent = m_spSV->parent_if_copy(elem);

	//	Get connected elements
		if(parent != NULL)
		{
			std::vector<TBaseElem*> vCoarseElem;
			CollectAssociated(vCoarseElem, multiGrid, parent, true);

			for(size_t i = 0; i < vCoarseElem.size(); ++i)
			{
			//	if shadowed, not in surface view of the requested level
				if(!m_spSV->is_surface_element(vCoarseElem[i])) continue;

			//	else this must be added to adjacend elements
				vAssElem.push_back(vCoarseElem[i]);
			}
		}
	}
}

}//	end of namespace

#endif
