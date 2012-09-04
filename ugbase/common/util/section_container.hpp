//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d10

#ifndef __UTIL__SECTION_CONTAINER__IMPL__
#define __UTIL__SECTION_CONTAINER__IMPL__

#include <cassert>
#include "section_container.h"

/*
I began work on reverse iterators. All associated code is in comments.
The method erase has not yet been adjusted.
I stopped implementing those reverse-iterators, since they would require quite
some additional instructions during element insertion (Each time the rbegin
iterator of the current and rend iterator of the next section would
have to be adjusted).
If reverse iterators are required, one should think about constructing them
directly from normal iterators in the rbegin / rend calls. This is done
similar in the current implementation of back.
*/

namespace ug
{

template <class TValue, class TContainer>
SectionContainer<TValue, TContainer>::
SectionContainer() : m_numElements(0)
{
}

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
clear()
{
//	clear the container. correct all sections afterwards.
	m_container.clear();
	m_numElements = 0;
	for(typename SectionVec::iterator iter = m_vSections.begin();
		iter != m_vSections.end(); ++iter)
	{
		Section& sec = *iter;
		sec.m_elemsBegin = sec.m_elemsEnd = m_container.end();
//		sec.m_elemsRBegin = sec.m_elemsREnd = m_container.rend();
		sec.m_numElements = 0;
	}
}

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
clear_section(int sectionIndex)
{
	assert((sectionIndex >= 0) && (sectionIndex < num_sections()) &&
			"ERROR in SectionContainer::clear_section(): bad sectionIndex");

	iterator iterStart = section_begin(sectionIndex);
	iterator iterEnd = section_end(sectionIndex);
	if(iterStart != iterEnd)
	{
		m_container.erase(iterStart, iterEnd);
		m_vSections[sectionIndex].m_elemsBegin = m_vSections[sectionIndex].m_elemsEnd;
//		m_vSections[sectionIndex].m_elemsRBegin = m_vSections[sectionIndex].m_elemsREnd;
		m_numElements -= m_vSections[sectionIndex].m_numElements;
		m_vSections[sectionIndex].m_numElements = 0;

	//	iterators of previous sections have to be adjusted
		for(int i = sectionIndex - 1; i >= 0; --i){
			if(m_vSections[i].m_numElements > 0){
				m_vSections[i].m_elemsEnd = m_vSections[sectionIndex].m_elemsBegin;
				break;//	we don't have to correct any more sections
			}
		}

	//	same for iterators of following sections
//		for(int i = sectionIndex + 1; i < num_sections(); ++i){
//			if(m_vSections[i].m_numElements > 0){
//				m_vSections[i].m_elemsREnd = m_vSections[sectionIndex].m_elemsRBegin;
//				break;//	we don't have to correct any more sections
//			}
//		}
	}
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::iterator
SectionContainer<TValue, TContainer>::
section_begin(int sectionIndex)
{
	if(sectionIndex < 0)
		return m_container.begin();
	else if(sectionIndex >= num_sections())
		return m_container.end();

	return m_vSections[sectionIndex].m_elemsBegin;
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::const_iterator
SectionContainer<TValue, TContainer>::
section_begin(int sectionIndex) const
{
	if(sectionIndex < 0)
		return m_container.begin();
	else if(sectionIndex >= num_sections())
		return m_container.end();

	return m_vSections[sectionIndex].m_elemsBegin;
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::iterator
SectionContainer<TValue, TContainer>::
section_end(int sectionIndex)
{
	if(sectionIndex >= num_sections() || sectionIndex < 0)
		return m_container.end();

	return m_vSections[sectionIndex].m_elemsEnd;
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::const_iterator
SectionContainer<TValue, TContainer>::
section_end(int sectionIndex) const
{
	if(sectionIndex >= num_sections() || sectionIndex < 0)
		return m_container.end();

	return m_vSections[sectionIndex].m_elemsEnd;
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::value_type&
SectionContainer<TValue, TContainer>::
front(int secIndex)
{
	assert((secIndex >= -1) && (secIndex < num_sections()) && "Bad section index");
	if(secIndex == -1)
		return m_container.front();
	return *m_vSections[secIndex].m_elemsBegin;
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::value_type&
SectionContainer<TValue, TContainer>::
back(int secIndex)
{
	assert((secIndex >= -1) && (secIndex < num_sections()) && "Bad section index");
	if(secIndex == -1)
		return m_container.back();

//	things are a little more complicated here, since there is no rbegin
//	check whether the last element is the last element of the underlying
//	container, too. If so, we'll return the last element of this container.
	if(m_vSections[secIndex].m_elemsEnd == m_container.end())
		return m_container.back();
	
//	since it is not, it points to the element in m_container, which
//	directly follows the last element of this section.
//	we'll thus create a reverse iterator on m_elemsEnd, increase it
//	once, and return the element, to which the iterator now points.
	iterator titer(m_vSections[secIndex].m_elemsEnd);
	--titer;
	return *titer;
}

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
add_sections(int num)
{
	for(int i = 0; i < num; ++i)
		m_vSections.push_back(Section(m_container.end(), m_container.end(),
									  //m_container.rend(), m_container.rend(),
									  0));
}

template <class TValue, class TContainer>
uint
SectionContainer<TValue, TContainer>::
num_elements(int sectionIndex) const
{
	assert((sectionIndex >= 0) &&
			"ERROR in SectionContainer::num_elements(): bad sectionIndex");

	if(sectionIndex >= num_sections())
		return 0;
	return m_vSections[sectionIndex].m_numElements;
}

template <class TValue, class TContainer>
typename SectionContainer<TValue, TContainer>::iterator
SectionContainer<TValue, TContainer>::
insert(const TValue& val, int sectionIndex)
{
	assert((sectionIndex >= 0) &&
			"ERROR in SectionContainer::insert(): bad sectionIndex");

	iterator nHandle;

//	check if a new section has to be created.
	if(sectionIndex >= num_sections())
		add_sections(sectionIndex - num_sections() + 1);

//	if the current section is empty, we have to initialize it.
	if(num_elements(sectionIndex) == 0)
	{
	//	it's the first element in this section.
	//	iterate through the sections and find the next one which is not empty
		int nextValidSection = -1;
		for(int i = sectionIndex + 1; i < num_sections(); ++i){
			if(num_elements(i) > 0){
				nextValidSection = i;
				break;
			}
		}

	//	we have to find the previous valid section, too
		int prevValidSection = -1;
		for(int i = sectionIndex - 1; i >= 0; --i){
			if(num_elements(i) > 0){
				prevValidSection = i;
				break;
			}
		}

	//	get the iterator before which the element shall be inserted.
		if(nextValidSection != -1){
			m_vSections[sectionIndex].m_elemsEnd = m_vSections[nextValidSection].m_elemsBegin;
		}
		else
			m_vSections[sectionIndex].m_elemsEnd = m_container.end();

	//	insert the element
		nHandle = m_vSections[sectionIndex].m_elemsBegin = m_container.insert(m_vSections[sectionIndex].m_elemsEnd, val);
		m_vSections[sectionIndex].m_numElements++;
//		m_vSections[sectionIndex].m_elemsRBegin = reverse_iterator(nHandle);
		
	//	adjust rend iterator of the next section
//		if(nextValidSection != -1)
//			m_vSections[nextValidSection].m_elemsREnd = m_vSections[sectionIndex].m_elemsRBegin;
			
	//	adjust the end-iterator of the previous valid section.
		if(prevValidSection != -1){
			m_vSections[prevValidSection].m_elemsEnd = m_vSections[sectionIndex].m_elemsBegin;
//			m_vSections[sectionIndex].m_elemsREnd = m_vSections[prevValidSection].m_elemsRBegin;
		}
//		else
//			m_vSections[sectionIndex].m_elemsREnd = m_container.rend();
	}
	else
	{
	//	the section is not empty. Simply insert the element before the sections end-iterator.
		nHandle = m_container.insert(m_vSections[sectionIndex].m_elemsEnd, val);
		m_vSections[sectionIndex].m_numElements++;
	}

	m_numElements++;

	return nHandle;
}

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
erase(const typename ug::SectionContainer<TValue, TContainer>::iterator& elemHandle, int sectionIndex)
{
	assert((sectionIndex >= 0) && (sectionIndex < num_sections()) &&
			"ERROR in SectionContainer::erase(): bad sectionIndex");

	iterator hNext;

//	check if the element is the first in the section.
//	if this is the case, we have to update the begin-handle of this section
	if(elemHandle == m_vSections[sectionIndex].m_elemsBegin)
	{
	//	get the previous valid section
		int prevValidSection = -1;
		{
			for(int i = sectionIndex - 1; i >= 0; --i)
			{
				if(num_elements(i) > 0)
				{
					prevValidSection = i;
					break;
				}
			}
		}

	//	erase the element
		hNext = m_container.erase(elemHandle);
		m_vSections[sectionIndex].m_numElements--;

	//	update the current section and the previous valid one.
		m_vSections[sectionIndex].m_elemsBegin = hNext;

		if(prevValidSection != -1)
			m_vSections[prevValidSection].m_elemsEnd = hNext;
	}
	else
	{
	//	erase the element
		hNext = m_container.erase(elemHandle);
		m_vSections[sectionIndex].m_numElements--;
	}

	m_numElements--;

//	check if hNext points to the end of the container. If so we will adjust the
//	m_elemsEnd iterator of the section.
//	included for safety reasons. Shouldn't be needed...

	if(hNext == m_container.end())
		m_vSections[sectionIndex].m_elemsEnd = m_container.end();
}

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
append(const SectionContainer& c)
{
	for(int i = 0; i < c.num_sections(); ++i){
		for(const_iterator iter = c.section_begin(i); iter != c.section_end(i); ++iter){
			insert(*iter, i);
		}
	}
}

}

#endif
