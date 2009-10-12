//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d10

#ifndef __UTIL__SECTION_CONTAINER__IMPL__
#define __UTIL__SECTION_CONTAINER__IMPL__

#include <cassert>

namespace ug
{

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
clear()
{
//	clear the container. correct all sections afterwards.
	m_container.clear();
	m_numElements = 0;
	for(typename SectionVec::iterator iter = m_vSections.begin(); iter != m_vSections.end(); ++iter)
	{
		Section& sec = *iter;
		sec.m_elemsBegin = sec.m_elemsEnd = m_container.end();
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
		m_numElements -= m_vSections[sectionIndex].m_numElements;
		m_vSections[sectionIndex].m_numElements = 0;

	//	end-iterators of previous sections have to be corrected
		for(int i = sectionIndex - 1; i >= 0; --i)
		{
			if(m_vSections[i].m_numElements > 0)
			{
				m_vSections[i].m_elemsEnd = m_vSections[sectionIndex].m_elemsBegin;
				break;//	we don't have to correct any more sections
			}
		}
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
typename SectionContainer<TValue, TContainer>::iterator
SectionContainer<TValue, TContainer>::
section_end(int sectionIndex)
{
	if(sectionIndex >= num_sections() || sectionIndex < 0)
		return m_container.end();

	return m_vSections[sectionIndex].m_elemsEnd;
}

template <class TValue, class TContainer>
void
SectionContainer<TValue, TContainer>::
add_sections(int num)
{
	for(int i = 0; i < num; ++i)
		m_vSections.push_back(Section(m_container.end(), m_container.end(), 0));
}

template <class TValue, class TContainer>
uint
SectionContainer<TValue, TContainer>::
num_elements(int sectionIndex)
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
		{
			for(int i = sectionIndex + 1; i < num_sections(); ++i)
			{
				if(num_elements(i) > 0)
				{
					nextValidSection = i;
					break;
				}
			}
		}

	//	we have to find the previous valid section, too
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

	//	get the iterator before which the element shall be inserted.
		if(nextValidSection != -1)
			m_vSections[sectionIndex].m_elemsEnd = m_vSections[nextValidSection].m_elemsBegin;
		else
			m_vSections[sectionIndex].m_elemsEnd = m_container.end();

	//	insert the element
		nHandle = m_vSections[sectionIndex].m_elemsBegin = m_container.insert(m_vSections[sectionIndex].m_elemsEnd, val);
		m_vSections[sectionIndex].m_numElements++;

	//	adjust the end-iterator of the previous valid section.
		if(prevValidSection != -1)
			m_vSections[prevValidSection].m_elemsEnd = m_vSections[sectionIndex].m_elemsBegin;
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
erase(const typename SectionContainer<TValue, TContainer>::iterator& elemHandle, int sectionIndex)
{
	assert((sectionIndex >= 0) && (sectionIndex < num_sections()) &&
			"ERROR in SectionContainer::erase(): bad sectionIndex");

	iterator hNext;

//	check if the element is the first in the section.
//	if this is the case, we have to update the begin-handle of this section and the end-handle of the previous one.
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

}

#endif
