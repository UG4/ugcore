/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__attached_list__
#define __H__UG__attached_list__

#include <iterator>
#include "attachment_pipe.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
//	predeclarations
template <class TAAEntry> class ConstAttachedElementListIterator;

////////////////////////////////////////////////////////////////////////////////
///	A special iterator which allows to iterate over elements in a AttachedElementList.
template <class TAAEntry>
class AttachedElementListIterator : public std::iterator<
											std::bidirectional_iterator_tag,
											typename TAAEntry::element>
{
	public:
		using element = typename TAAEntry::element;
		using iterator = AttachedElementListIterator;

		AttachedElementListIterator() : m_curElem(nullptr)	{}
		AttachedElementListIterator(element curElem, const TAAEntry& aaEntry) :
			m_aaEntry(aaEntry), m_curElem(curElem)	{}
		AttachedElementListIterator(const AttachedElementListIterator& cpy) :
			m_aaEntry(cpy.m_aaEntry), m_curElem(cpy.m_curElem)	{}

		const iterator& operator=(const iterator& iter)
		{
			m_aaEntry = iter.m_aaEntry;
			m_curElem = iter.m_curElem;
			return *this;
		}

	///	note that the * operator is read only.
		element operator * () const	{return m_curElem;}
		element* operator -> () const	{return &m_curElem;}

		iterator operator ++ ()		{m_curElem = m_aaEntry[m_curElem].next; return *this;}
		iterator operator ++ (int)
		{//	post-increment
			iterator iter = *this;
			m_curElem = m_aaEntry[m_curElem].next;
			return iter;
		}

		iterator operator -- ()		{m_curElem = m_aaEntry[m_curElem].prev; return *this;}
		iterator operator -- (int)
		{
			iterator iter = *this;
			m_curElem = m_aaEntry[m_curElem].prev;
			return iter;
		}

	//	note that each element may only be in the list once.
		bool operator == (const iterator& iter) const	{return m_curElem == iter.m_curElem;}
		bool operator != (const iterator& iter) const	{return m_curElem != iter.m_curElem;}

	private:
		TAAEntry	m_aaEntry;
		element		m_curElem;

	friend class ConstAttachedElementListIterator<TAAEntry>;
};

////////////////////////////////////////////////////////////////////////////////
///	A special iterator which allows to iterate over elements in a AttachedElementList.
template <class TAAEntry>
class ConstAttachedElementListIterator : public std::iterator<
											std::bidirectional_iterator_tag,
											const typename TAAEntry::element>
{ // ø todo iterator
	public:
		using element = typename TAAEntry::element;
		using iterator = ConstAttachedElementListIterator;

		ConstAttachedElementListIterator() : m_curElem(nullptr)	{}
		ConstAttachedElementListIterator(element curElem, const TAAEntry& aaEntry) :
			m_aaEntry(aaEntry), m_curElem(curElem)	{}
		ConstAttachedElementListIterator(const ConstAttachedElementListIterator& it) :
			m_aaEntry(it.m_aaEntry), m_curElem(it.m_curElem)	{}
		ConstAttachedElementListIterator(const AttachedElementListIterator<TAAEntry>& it) :
			m_aaEntry(it.m_aaEntry), m_curElem(it.m_curElem)	{}

		const iterator& operator=(const iterator& iter)
		{
			m_aaEntry = iter.m_aaEntry;
			m_curElem = iter.m_curElem;
			return *this;
		}

	///	note that the * operator is read only.
		element operator*() const {return m_curElem;}
		const element* operator->() const {return &m_curElem;}

		iterator operator ++ () {m_curElem = m_aaEntry[m_curElem].next; return *this;}
		iterator operator ++ (int)
		{//	post-increment
			iterator iter = *this;
			m_curElem = m_aaEntry[m_curElem].next;
			return iter;
		}

		iterator operator--() {m_curElem = m_aaEntry[m_curElem].prev; return *this;}
		iterator operator--(int)
		{
			iterator iter = *this;
			m_curElem = m_aaEntry[m_curElem].prev;
			return iter;
		}

	//	note that each element may only be in the list once.
		bool operator==(const iterator& iter) const	{return m_curElem == iter.m_curElem;}
		bool operator!=(const iterator& iter) const	{return m_curElem != iter.m_curElem;}

	private:
		TAAEntry	m_aaEntry;
		element		m_curElem;
};

////////////////////////////////////////////////////////////////////////////////
///	A linked list of elements living in an attachment
/**	This is a highly specialized linked list mainly (only!) used in the Grid,
 * SubsetHandler, Selector, ... classes to maintain a list of GridObjects.
 *
 * Only elements registered at the underlying attachment pipe may be added
 * to the list. Note that each element may only be added once.
 *
 * Make sure that the attachment pipe on which the list operates
 * exists at least as long as the list itself, or call set_pipe(nullptr),
 * if the pipe is erased before.
 *
 * This list assumes that TAttachmentPipe::element is a pointer type.
 *
 * Note that instances of this class are never lightweight, since memory
 * for a potentially full list is allocated in the attachment pipe.
 *
 * If the list operates on a shared attachment, the one who created the
 * attachment is responsible for releasing it.
 *
 * \warning	Special care has to be taken with this type of list. It is e.g. not
 * 			possible to insert the same element multiple times into the list.
 * 			Instead upon new insertion, the old entry will be erased.
 */
template <class TAttachmentPipe>
class AttachedElementList
{
	public:
		using element = typename TAttachmentPipe::element;

		struct Entry{
			Entry() : prev(nullptr), next(nullptr)	{}
			Entry(const element& p, const element& n) : prev(p), next(n) {}
			element prev;
			element next;
		};

		using AEntry = Attachment<Entry>;
		using ElemHandler = typename TAttachmentPipe::ElementHandler;
		using AAEntry = AttachmentAccessor<element, AEntry, ElemHandler>;

		using iterator = AttachedElementListIterator<AAEntry>;
		using const_iterator = ConstAttachedElementListIterator<AAEntry>;

	public:
		AttachedElementList(TAttachmentPipe* pipe = nullptr) :
			m_pipe(nullptr),
			m_aEntry("AttachedElementList_Entry", false),
			m_front(nullptr),
			m_back(nullptr),
			m_bSharedAttachment(false)
		{
			if(pipe)
				set_pipe(pipe);
		}

	///	Note that auto-copy on aEntry has to be disabled.
		AttachedElementList(AEntry aEntry) :
			m_pipe(nullptr),
			m_aEntry(aEntry),
			m_front(nullptr),
			m_back(nullptr),
			m_bSharedAttachment(false)
		{}

	///	Note that auto-copy on aEntry has to be disabled.
		AttachedElementList(TAttachmentPipe* pipe, AEntry aEntry) :
			m_pipe(nullptr),
			m_aEntry(aEntry),
			m_front(nullptr),
			m_back(nullptr),
			m_bSharedAttachment(true)
		{
			if(pipe)
				set_pipe(pipe);
		}

		AttachedElementList(const AttachedElementList& ael) : m_pipe(nullptr)
		{
			m_bSharedAttachment = ael.m_bSharedAttachment;
			if(m_bSharedAttachment)
				m_aEntry = ael.m_aEntry;

			set_pipe(ael.m_pipe);
			if(ael.m_front){
				iterator iter(ael.m_front, m_aaEntry);
				while(*iter){
					push_back(*iter);
					++iter;
				}
			}
		}

		~AttachedElementList()
		{
			set_pipe(nullptr);
		}

		const AttachedElementList& operator=(const AttachedElementList& ael)
		{
			clear();

			if(ael.m_bSharedAttachment)
				set_pipe(ael.m_pipe, ael.m_aEntry);
			else
				set_pipe(ael.m_pipe);

			if(ael.m_front){
				iterator iter(ael.m_front, m_aaEntry);
				while(*iter){
					push_back(*iter);
					++iter;
				}
			}

			return *this;
		}

	///	set the attachment pipe on which the list shall operate
		void set_pipe(TAttachmentPipe* pipe)
		{
			if(!m_bSharedAttachment && m_pipe)
				m_pipe->detach(m_aEntry);
			m_aaEntry.invalidate();

			m_pipe = pipe;
			if(m_pipe){
				if(!m_pipe->has_attachment(m_aEntry))
					m_pipe->attach(m_aEntry);
				m_aaEntry.access(*m_pipe, m_aEntry);
			}

			m_front = m_back = nullptr;
		}

	///	Sets the pipe and a shared entry-attachment on which the list will operate
	/** Note that auto-copy on aEntry has to be disabled.*/
		void set_pipe(TAttachmentPipe* pipe, AEntry aEntry)
		{
			if(!m_bSharedAttachment && m_pipe)
				m_pipe->detach(m_aEntry);
			m_aaEntry.invalidate();

			m_pipe = pipe;
			m_aEntry = aEntry;
			m_bSharedAttachment = true;

			if(m_pipe){
			//	since we use a shared attachment in this case, it may already be
			//	attached to the pipe.
				if(!m_pipe->has_attachment(m_aEntry))
					m_pipe->attach(m_aEntry);
				m_aaEntry.access(*m_pipe, m_aEntry);
			}

			m_front = m_back = nullptr;
		}

	///	clears the list. begin() == end() afterwards.
		void clear()
		{
			if(m_pipe){
				iterator iter = begin();
				while(iter != end()){
					element elem = *iter;
					iter++;
					m_aaEntry[elem] = Entry();
				}

				m_front = m_back = nullptr;
			}
		}


	///	retunrs true if the list is empty
		bool empty() const				{return m_front == nullptr;}

	///	returns the first element in the list
		element front()					{return m_front;}
	///	returns the last element in the list
		element back()					{return m_back;}

	///	returns the first element in the list (const)
		const element front() const		{return m_front;}
	///	returns the last element in the list (const)
		const element back() const		{return m_back;}

	///	pushes an element to the end of the list
		void push_back(const element& elem)
		{
			m_aaEntry[elem] = Entry(m_back, nullptr);
			if(empty())
				m_front = elem;
			else
				m_aaEntry[m_back].next = elem;
			m_back = elem;
		}

	///	inserts an element right before the specified position.
	/**	returns the iterator of the newly inserted element.
	 */
		iterator insert(iterator position, const element& elem)
		{
			if(empty() || !(*position))
				push_back(elem);
			else{
				Entry& entry = m_aaEntry[*position];

				if(entry.prev){
					m_aaEntry[entry.prev].next = elem;
					m_aaEntry[elem].prev = entry.prev;
				}
				else{
					m_aaEntry[elem].prev = nullptr;
					m_front = elem;
				}

				m_aaEntry[elem].next = *position;
				entry.prev = elem;
			}

			return get_iterator(elem);
		}

	///	erases the element at the given position
	/**	returns an iterator to the element directly behind position.
	 */
		iterator erase(iterator position)
		{
			Entry& entry = m_aaEntry[*position];
			if(entry.prev)
				m_aaEntry[entry.prev].next = entry.next;
			else
				m_front = entry.next;

			if(entry.next)
				m_aaEntry[entry.next].prev = entry.prev;
			else
				m_back = entry.prev;

		//	get next element and reset entry.
			element nextElem = entry.next;
			entry = Entry();
			return iterator(nextElem, m_aaEntry);
		}

	///	erases a sequence of entries
		iterator erase(iterator begin, iterator end)
		{
			iterator iter = begin;
			while(iter != end){
				iterator titer = iter++;
				erase(titer);
			}
			return iter;
		}

	///	returns the iterator to the given element
	/**	Note that the return value is undefined if element is not
	 * a member of the list!
	 */
		iterator get_iterator(const element& elem)
		{
			return iterator(elem, m_aaEntry);
		}

	///	returns a const pointer to an element.
	/**	This pointer is valid until the content of the list is changed.*/
		element const* get_pointer_to_element(const element& elem)
		{
			if(elem == m_front)
				return &m_front;
			return &m_aaEntry[m_aaEntry[elem].prev].next;
		}
		
	///	returns true if the element is in the list
		bool is_in_list(const element& elem)
		{
			return (m_front == elem) || (m_aaEntry[elem].prev != nullptr);
		}
	///	returns an iterator to the beginning of the sequence.
		iterator begin()				{return iterator(m_front, m_aaEntry);}

	///	returns an iterator to the end of the sequence.
	/**	Note that this is the iterator to the element behind the last one.*/
		iterator end()					{return iterator(nullptr, m_aaEntry);}

	///	returns an iterator to the beginning of the sequence.
		const_iterator begin() const	{return const_iterator(m_front, m_aaEntry);}

	///	returns an iterator to the end of the sequence.
	/**	Note that this is the iterator to the element behind the last one.*/
		const_iterator end() const		{return const_iterator(nullptr, m_aaEntry);}

	private:
	//	the attachment pipe on which we'll operate
		TAttachmentPipe*	m_pipe;

	//	The entry attachment
		AEntry		m_aEntry;

	//	the attachment accessor with which the entries are accessed
		AAEntry		m_aaEntry;

	//	front and back elements
		element		m_front;
		element		m_back;

	//	tells whether the entry attachment is shared with other lists
		bool		m_bSharedAttachment;
};



}//	end of namespace

#endif
