//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d10

#ifndef __UTIL__SECTION_CONTAINER__
#define __UTIL__SECTION_CONTAINER__

#include <vector>
#include "../types.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SectionContainer
///	A container that is divided into different sections.
/**
 * This container can be used to store values that have a common base,
 * but can be sorted in different sections. The SectionContainer supplies
 * you with interfaces to iterate over a certain section or to iterate
 * over the complete container.
 * The SectionContainer works on containers of the standard template library.
 * Be sure, that TContainer holds values of type TValue.
 * TContainer has to support the container operations of a std::list.
 */
template <class TValue, class TContainer>
class SectionContainer
{
	public:
		typedef TContainer							Container;
		typedef typename Container::iterator		iterator;
		typedef typename Container::const_iterator	const_iterator;

	public:
		SectionContainer();

		void clear();
		void clear_section(int sectionIndex);

		iterator insert(const TValue& val, int sectionIndex);
		void erase(const iterator& iter, int sectionIndex);

		inline iterator begin()	{return m_container.begin();}
		inline iterator end()	{return m_container.end();}

		const_iterator begin() const	{return m_container.begin();}
		const_iterator end() const		{return m_container.end();}

	/**if the section is empty section_begin and section_end return the same iterators.
	 * However, no assumptions on the positions of these iterators should be made.*/
		iterator section_begin(int sectionIndex);

		const_iterator section_begin(int sectionIndex) const;

	/**if the section is empty section_begin and section_end return the same iterators.
	   However, no assumptions on the positions of these iterators should be made.*/
		iterator section_end(int sectionIndex);

		const_iterator section_end(int sectionIndex) const;

		uint num_elements(int sectionIndex) const;
		inline uint num_elements() const	{return m_numElements;}
		inline int num_sections() const		{return m_vSections.size();}

	protected:
		void add_sections(int num);

	protected:
		struct Section
		{
			Section()	{}
			Section(const iterator& elemsBegin, const iterator& elemsEnd, int numElems) :
				m_elemsBegin(elemsBegin), m_elemsEnd(elemsEnd), m_numElements(numElems)
				{}

			iterator	m_elemsBegin;
			iterator	m_elemsEnd;
			uint		m_numElements;
		};

		typedef std::vector<Section>	SectionVec;

	protected:
		Container		m_container;
		SectionVec		m_vSections;
		uint			m_numElements;
};

}//	end of namespace

//	include implementation
#include "section_container.hpp"

#endif
