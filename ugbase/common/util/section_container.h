#ifndef __UTIL__SECTION_CONTAINER__
#define __UTIL__SECTION_CONTAINER__

#include <vector>
#include "../types.h"

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

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
		typedef TValue								value_type;
		typedef TContainer							Container;
		typedef typename Container::iterator		iterator;
		typedef typename Container::const_iterator	const_iterator;
		//typedef typename Container::reverse_iterator		reverse_iterator;
		//typedef typename Container::const_reverse_iterator	const_reverse_iterator;
		
	public:
		SectionContainer();

		void clear();
		void clear_section(int sectionIndex);

		iterator insert(const TValue& val, int sectionIndex);
		void erase(const iterator& iter, int sectionIndex);

		inline iterator begin()	{return m_container.begin();}
		inline iterator end()	{return m_container.end();}

		inline const_iterator begin() const	{return m_container.begin();}
		inline const_iterator end() const	{return m_container.end();}

	/**if the section is empty section_begin and section_end return the same iterators.
	 * However, no assumptions on the positions of these iterators should be made.*/
		iterator section_begin(int sectionIndex);

		const_iterator section_begin(int sectionIndex) const;

	/**if the section is empty section_begin and section_end return the same iterators.
	   However, no assumptions on the positions of these iterators should be made.*/
		iterator section_end(int sectionIndex);

		const_iterator section_end(int sectionIndex) const;

	///	returns the first entry in the given section.
	/**	use index -1 to get the first entry of the complete chain.
	 *	Make sure to only call this method if there are elements in the
	 *	given section at all.*/
		value_type& front(int secIndex = -1);

	///	returns the last entry in the given section.
	/**	use index -1 to get the last entry of the complete chain.
	 *	Make sure to only call this method if there are elements in the
	 *	given section at all.*/
		value_type& back(int secIndex = -1);
			
		uint num_elements(int sectionIndex) const;
		inline uint num_elements() const	{return m_numElements;}
		inline int num_sections() const		{return m_vSections.size();}

	///	returns the container for raw access.
	/**	Use this method with extreme care. Changes to the elements
	 * and to the layout of the container may most likely result in a
	 * corruption of the SectionContainer.
	 */
		inline Container& get_container()			{return m_container;}

	///	appends the elements of the given container to the current one
	/**	Note that the append operation is performed for each section separately.
	 * \warning	This method should not be used if the underlying container and
	 * 			the given one operate on the same AttachedElemList. Otherwise
	 * 			severe side effect will occur! Use transfer_elements instead!*/
		void append(const SectionContainer& c);

	///	takes all elements from the given section container and transfers them to this one.
		void transfer_elements(SectionContainer& c);

	protected:
		void add_sections(int num);

	protected:
		struct Section
		{
			Section()	{}
			Section(const iterator& elemsBegin, const iterator& elemsEnd,
					//const reverse_iterator& elemsRBegin,
					//const reverse_iterator& elemsREnd,
					int numElems) :
				m_elemsBegin(elemsBegin), m_elemsEnd(elemsEnd),
				//m_elemsRBegin(elemsRBegin), m_elemsREnd(elemsREnd),
				m_numElements(numElems)
				{}

			iterator	m_elemsBegin;
			iterator	m_elemsEnd;
			//reverse_iterator 	m_elemsRBegin;
			//reverse_iterator 	m_elemsREnd;
			uint		m_numElements;
		};

		typedef std::vector<Section>	SectionVec;

	protected:
		Container		m_container;
		SectionVec		m_vSections;
		uint			m_numElements;
};

// end group ugbase_common_util
/// \}

}//	end of namespace

//	include implementation
#include "section_container.hpp"

#endif
