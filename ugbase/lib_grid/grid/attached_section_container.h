// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 01.04.2011 (m,d,y)

#ifndef __H__UG__attached_section_container__
#define __H__UG__attached_section_container__

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	AttachedSectionContainer
///	A container that is divided into different sections, living in an attachment.
/**
 * This container can be used to store values that have a common base,
 * but can be sorted in different sections. The SectionContainer supplies
 * you with interfaces to iterate over a certain section or to iterate
 * over the complete container.
 * The section container stores its entries directly in an attachment.
 */
template <class TValue, class TAttachmentPipe>
class AttachedSectionContainer
{
	public:
		typedef TValue								value_type;

	public:
		AttachedSectionContainer();

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

}//	end of namespace

#endif
