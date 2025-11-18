/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __UTIL__ATTACHMENT_PIPE__
#define __UTIL__ATTACHMENT_PIPE__

#include <list>
#include <vector>
#include <stack>
#include <cassert>
#include "common/static_assert.h"
#include "common/types.h"
#include "common/util/uid.h"
#include "common/util/hash.h"
#include "common/ug_config.h"


namespace ug
{

// PREDECLARATIONS
class IAttachment;
class IAttachmentDataContainer;


// CONSTANTS
enum ATTACHMENT_CONSTANTS
{
	INVALID_ATTACHMENT_INDEX = 0xFFFFFFFF
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	IAttachedDataContainer
///	the interface for an attachment-data-container.
/** In order to use an attachment-data-container you have to supply several type defintins.
* Take a look at the derived generic class AttachmentDataContainer<T> to see which
* type defintion and operations have to be supplied in addition to the interface-methods.
* if possible you should always use the derived generic class AttachmentDataContainer<T>
* instead creating inheriting your own version of IAttachedDataContainer.
*/
class UG_API IAttachmentDataContainer
{
	public:
		virtual ~IAttachmentDataContainer()	= default;

		virtual void resize(size_t iSize) = 0;///< resize the data array
		virtual size_t size() = 0;///< returns the size of the data-array.
		virtual void copy_data(size_t indFrom, size_t indTo) = 0;///< copy data from entry indFrom to entry indTo.
		virtual void reset_entry(size_t index) = 0;///< resets the entry to its default value.
		
	/**	copies entries from the this-container to the specified target container.
	 * For the i-th entry in the destination-container, pIndexMap has to contain
	 * the index of the associated source entry in this container.
	 * num specifies the number of entries to be copied.
	 * Make sure, that pDest can hold 'num' elements.
	 *
	 * pDestCon has to have the same or a derived dynamic type as the container
	 * on which this method is called.*/
		virtual void copy_to_container(IAttachmentDataContainer* pDestCon,
									   int* indexMap, int num) const = 0;
		
	/** defragment should clear the containers data from unused entries.
	 * pNewIndices should be an array of indices, wich holds a new index for each
	 * entry in IAttachedDataContainer. pNewIndices has thus to hold as many indices
	 * as there are entries in IAttachedDataContainer. If an entry shall not appear
	 * in the defragmented container, its new index has to be set to INVALID_ATTACHMENT_INDEX.
	 * numValidElements has to specify the size of the defragmented container -
	 * it thus has to equal the number of valid indices in pNewIndices.*/
		virtual void defragment(size_t* pNewIndices, size_t numValidElements) = 0;

	///	returns the size in bytes, which the container occupies
	/** Mainly for debugging purposes.*/
		virtual size_t occupied_memory() = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////
///	define reference and const reference types for attachment values.
/**	The default traits should be fine in all cases. A specialization for bool exists.
 */
template <class TValue>
struct attachment_value_traits{
	using reference = TValue&;
	using const_reference = const TValue&;
};

///	specialization of attachment_value_traits for the bool type
template <>
struct attachment_value_traits<bool>{
	using reference = std::vector<bool>::reference;
	using const_reference = const std::vector<bool>::reference;
};

/* THOUGHTS
*	AttachmentDataContainer<T> should probably be defined in another header, since it is somehow specialised for libGrid.
*	same goes for Attachment<T>
*/
////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachedDataContainer
///	A generic specialization of IAttachedDataContainer.
/** This template-class not only simplifies the creation of custom containers,
 * it also defines some types, operators and values, which are essential to use an AttachmentDataContainer with libGrid.
 * In particular libGrids AttachmentAccessors require these definitions.
 */
template <class T> class UG_API AttachmentDataContainer : public IAttachmentDataContainer
{
	private:
		using ClassType = AttachmentDataContainer<T>;
		using DataContainer = std::vector<T>;
		using TRef = typename attachment_value_traits<T>::reference;
		using TConstRef = typename attachment_value_traits<T>::const_reference;

	public:
		using ValueType = T;

		AttachmentDataContainer(const T& defaultValue = T())	: m_defaultValue(defaultValue)	{}

		~AttachmentDataContainer() override {m_vData.clear();}

		void resize(size_t iSize) override {
				if(iSize > 0)
					m_vData.resize(iSize, m_defaultValue);
				else
					m_vData.clear();
			}

		size_t size() override {return m_vData.size();}

		void copy_data(size_t indFrom, size_t indTo) override {m_vData[indTo] = m_vData[indFrom];}

		void reset_entry(size_t index) override {
				m_vData[index] = m_defaultValue;
			}

		void defragment(size_t* pNewIndices, size_t numValidElements) override {
				size_t numOldElems = size();
				DataContainer vDataOld = m_vData;
				m_vData.resize(numValidElements);
				for(size_t i = 0; i < numOldElems; ++i)
				{
					size_t nInd = pNewIndices[i];
					if(nInd != INVALID_ATTACHMENT_INDEX)
						m_vData[nInd] = vDataOld[i];
				}
			}
	
	/**	copies entries from the this-container to the container
	 *	specified by pDestCon.
	 *	For the i-th entry in the target container, pIndexMap has to contain
	 *	the index of the associated source entry in this container.
	 *	num specifies the number of entries to be copied.
	 *	Make sure, that pDest can hold 'num' elements.
	 *
	 *	pDestCon has to have the same or a derived dynamic type as the container
	 *	on which this method is called.*/
		void copy_to_container(IAttachmentDataContainer* pDestCon,
		                       int* indexMap, int num) const override {
				ClassType* destConT = dynamic_cast<ClassType*>(pDestCon);
				assert(destConT && "Type of pDestBuf has to be the same as the"
						"type of this buffer");

				if(!destConT)
					return;

				DataContainer& dest = destConT->get_data_container();
				for(int i = 0; i < num; ++i)
					dest[i] = m_vData[indexMap[i]];
			}


	///	returns the memory occupied by the container
		size_t occupied_memory() override {
			return m_vData.capacity() * sizeof(T);
		}

		inline TConstRef get_elem(size_t index) const		{return m_vData[index];}
		inline TRef get_elem(size_t index)					{return m_vData[index];}
		inline TConstRef operator[] (size_t index) const	{return m_vData[index];}
		inline TRef operator[] (size_t index)				{return m_vData[index];}

	///	swaps the buffer content of associated data
		void swap(AttachmentDataContainer<T>& container) noexcept {m_vData.swap(container.m_vData);}

	protected:
		DataContainer& get_data_container()			{return m_vData;}
		//inline const T* get_ptr() const			{return &m_vData.front();}
		//inline T* get_ptr()						{return &m_vData.front();}
		
	protected:
		DataContainer	m_vData;
		T				m_defaultValue;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	IAttachment
///	the interface for attachments.
/** Attachments can be attached to an AttachmentPipe and thus enhance the pipes elements by data,
 * whose type, container and behavior is defined by the Attachment itself.
 * In order to use an Attachment with libGrid (in particular with libGrids AttachmentAccessors),
 * derivatives of IAttachment have to feature some special type defintion (see Attachment<T> for more information).
 * Whenever possible you should use the template-derivative Attachment<T> instead of IAttachment.
 */
class UG_API IAttachment : public UID
{
	public:
		IAttachment() : m_name("undefined")   {}
		IAttachment(const char* name) : m_name(name)
			{assert(m_name);}

		~IAttachment() override = default;
		virtual IAttachment* clone() = 0;
		virtual IAttachmentDataContainer*	create_container() = 0;
		virtual bool default_pass_on_behaviour() const = 0;

		const char* get_name()  {return m_name;}    ///< should only be used for debug purposes.

	protected:
		const char* m_name; //only for debug
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	Attachment
/// A generic specialization of IAttachment
/** This class is intended to simplify the process of Attachment creation.
 * Note that there are type definitions, which are required by libGrids AttachmentAccessors.
 */
template <class T> class UG_API Attachment : public IAttachment
{
	public:
		using ContainerType = AttachmentDataContainer<T>;
		using ValueType = T;

		Attachment() : m_passOnBehaviour(false)    {}
		Attachment(bool passOnBehaviour) : m_passOnBehaviour(passOnBehaviour)    {}
		Attachment(const char* name) : IAttachment(name), m_passOnBehaviour(false)   		{}
		Attachment(const char* name, bool passOnBehaviour) : IAttachment(name), m_passOnBehaviour(passOnBehaviour)	{}

		~Attachment() override = default;

		IAttachment* clone() override {IAttachment* pA = new Attachment<T>; *pA = *this; return pA;}
		IAttachmentDataContainer* create_container() override {return new ContainerType;}
		bool default_pass_on_behaviour() const override {return m_passOnBehaviour;}
		IAttachmentDataContainer* create_container(const T& defaultValue)	{return new ContainerType(defaultValue);}

	protected:
		bool m_passOnBehaviour;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachmentEntry
///	This struct is used by AttachmentPipe in order to manage its attachments
struct AttachmentEntry
{
	AttachmentEntry() : m_pAttachment(nullptr), m_pContainer(nullptr), m_userData(0)	{}
	AttachmentEntry(IAttachment* pAttachment, IAttachmentDataContainer* pContainer, uint userData = 0) :
			m_pAttachment(pAttachment), m_pContainer(pContainer), m_userData(userData)	{}

	IAttachment*				m_pAttachment;
	IAttachmentDataContainer*	m_pContainer;
	uint						m_userData;
};


////////////////////////////////////////////////////////////////////////////////////////////////
///	define the interface that enables you to use your own types as element-types in an AttachmentPipe.
/** By creating a template specialization for your own element types, you can
 * use arbitrary types as element-types in an AttachmentPipe.
 */
template<class TElem, class TElemHandler>
class attachment_traits
{
	public:
		using ElemRef = TElem&;
		using ElemPtr = void*;
		using ConstElemPtr = const void*;
		using ElemHandlerPtr = void*;
		using ConstElemHandlerPtr = const void*;

		using element_iterator = void;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)					{}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)					{}
		static inline uint get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem)			{return INVALID_ATTACHMENT_INDEX;/*STATIC_ASSERT(0, INVALID_ATTACHMENT_TRAITS);*/}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, size_t index)	{/*STATIC_ASSERT(0, INVALID_ATTACHMENT_TRAITS);*/}
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachmentPipe
///	Handles data which has been attached to the pipe using callbacks for the element.
/** The AttachmentPipe can be used to attach data to a collection of elements.
 * Elements have to be registered at the AttachmentPipe. Using the methods
 * defined in the attachment_traits template class, registered elements are
 * associated with their data-entries.
 *
 * TODO:
 * enum Options : CopyAllowed
 *
 * - copy_values(TElem from, TElem to)
 * - swap_entry_indices()?!?
 */
template<class TElem, class TElemHandler>
class UG_API AttachmentPipe
{
	public:
		using element = TElem;
		using ElementHandler = TElemHandler;
		using AttachmentEntryContainer = std::list<AttachmentEntry>;
		using AttachmentEntryIterator = AttachmentEntryContainer::iterator;
		using ConstAttachmentEntryIterator = AttachmentEntryContainer::const_iterator;
		using AttachmentEntryIteratorHash = Hash<uint, AttachmentEntryIterator>;
		using atraits = attachment_traits<TElem, TElemHandler>;

	public:
		AttachmentPipe();
		AttachmentPipe(typename atraits::ElemHandlerPtr pHandler);
		~AttachmentPipe();
		
		inline typename atraits::ElemHandlerPtr get_elem_handler()		{return m_pHandler;}
		
	///	calls both clear_elements and clear_attachments.
		void clear();
	///	clears elements and associated data but keeps registered attachments
		void clear_elements();
	///	clears registered attachments and associated data but keeps existing elements.
		void clear_attachments();

	///	Reserves memory for element- and data-entries.
	/**	This method reserves memory for data and elements.
	 * Note that numElems specifies the total amount of elements
	 * for which data shall be reserved. If numElems is smaller than
	 * num_elements, then nothing will be done.
	 * Note that a call to reserve does not change num_elements nor
	 * num_data_entries.
	 *
	 * If you want to reserve space for additional N elements, please
	 * call reserve with numElems = num_elements() + N.
	 *
	 * Note that reserve is optional and only serves performance benefits.
	 */
		void reserve(size_t numElems);

	///	Registers a new element at the attachment pipe.
	/**	Registers the element and reserves memory for all registered
	 * attachments. Note that this method may be faster if memory has
	 * previously been reserved using reserve_new.
	 */
		void register_element(TElem elem);

	///	Unregisters the given element.
	/**	Unregisters the given element but does not yet clear associated data.
	 * Indeed the memory is kept until defragment is called or new elements
	 * reuse the memory.
	 */
		void unregister_element(const TElem& elem);

	/**	Aligns data with elements and removes unused data-memory.*/
		void defragment();

	/**\brief attaches a new data-array to the pipe.
	 *
	 * Attachs a new attachment and creates a container which holds the
	 * data associated with the elements through this attachment.
	 * Several overloads exist: You may pass an attachment together with
	 * a default value and with an option constant or simply the attachment
	 * together with an optional option constant. Not that in the first
	 * case the concrete type of the attachment has to be known, whereas
	 * in the second case only the interface IAttachment has to be specified.
	 *
	 * The option constant is not used by the attachment system but may be
	 * used by a user to store a constant with each attachment.
	 * \{
	 */
		template <class TAttachment>
		void attach(TAttachment& attachment,
					const typename TAttachment::ValueType& defaultValue,
					uint options);

		void attach(IAttachment& attachment, uint options = 0);
	/**	\}	*/

	///	Removes the data associated with the given attachment from the pipe.
		void detach(IAttachment& attachment);

	///	Returns true if the given attachment is currently attached to the pipe.
		bool has_attachment(IAttachment& attachment) const;

	///	Lets you access the raw data array associated with the given attachment.
	/**	If you access several arrays through this method, it is guaranteed that
	 * the data associated with one specific object are found at the same indices
	 * in those arrays.
	 *
	 * Note that if you access the data arrays after a call to defragment, then
	 * the i-th data-entry corresponds to the i-th element.
	 *
	 * Please not that the pointer may be invalidated through the following operations:
	 * 		- register_element
	 * 		- defragment
	 * 		- clear, clear_elements, clear_attachments
	 */
		template <class TAttachment>
		typename TAttachment::ValueType*
		get_data_array(TAttachment& attachment);
		
	/**	\brief Returns the data container managing the data array for the given attachment.
	 * \{ */
		IAttachmentDataContainer* get_data_container(IAttachment& attachment) const;

		template <class TAttachment>
		typename TAttachment::ContainerType*
		get_data_container(TAttachment& attachment);
	/**	\} */

		inline ConstAttachmentEntryIterator attachments_begin() const	{return m_attachmentEntryContainer.begin();}
		inline ConstAttachmentEntryIterator attachments_end() const		{return m_attachmentEntryContainer.end();}

	///	Returns the number of registered elements
		inline size_t num_elements() const		{return m_numElements;}
	///	Returns the size of the associated data arrays
	/**	Note: If the pipe is fragmented, then num_elements and num_data_entries
	 * differ. If the pipe however isn't fragmented, then both values are the same.
	 */
		inline size_t num_data_entries() const	{return m_numDataEntries;}

	///	Returns whether the attachment pipe is fragmented.
	/**	The pipe gets fragmented whenever elements are erased.
	 * Through the creation of new elements the pipe may be automatically
	 * defragmented.
	 * Through a call to defragment() the pipe can be defragmented at any time.
	 */
		inline bool is_fragmented() const		{return m_numElements != m_numDataEntries;}

		void reset_values(size_t dataIndex);///< fills the attached data at the given index with the default values.

	protected:
		void resize_attachment_containers(size_t newSize);
		void grow_attachment_containers(size_t newMinSize);
		inline size_t get_container_size();

	protected:
		using UINTStack = std::stack<size_t>;

	protected:
		AttachmentEntryContainer	m_attachmentEntryContainer;
		AttachmentEntryIteratorHash	m_attachmentEntryIteratorHash;

		UINTStack		m_stackFreeEntries;	///< holds indices to free entries.

		size_t			m_numElements;
		size_t			m_numDataEntries;
		size_t			m_containerSize; ///< total size of containers.
		
		typename atraits::ElemHandlerPtr	m_pHandler;
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachmentAccessor
///	Used to access data that has been attached to an attachment pipe.
/** Once initialized, an AttachmentAccessor can be used to access the data stored
 * in the given Attachment at the given AttachmentPipe.
 * The reference type of the associated value is taken from attachment_value_traits.
 * By default this is the standard reference type.
 *
 * To initialize an AttachmentAccessor, you may either use its constructor or
 * its access method, which returns false if the specified attachment is not present
 * in the specified attachment pipe.
 *
 * Note that data-access using an attachment accessor is cheap. Setting up
 * a new attachment accessor however involves some work. While this is generally
 * fast, too, it would introduce an unnecessary overhead inside inner loops
 * or frequently called methods. You should thus try to minimize calls to access
 * or to the accessors constructor.
 */
template <class TElem, class TAttachment, class TElemHandler>
class UG_API AttachmentAccessor
{
	public:
		using attachment = TAttachment;
		using element = TElem;
		using ValueType = typename TAttachment::ValueType;
		using ContainerType = typename TAttachment::ContainerType;
		using ElemHandler = TElemHandler;
		using atraits = attachment_traits<TElem, TElemHandler>;
		using attachment_pipe = AttachmentPipe<TElem, TElemHandler>;

	public:
		AttachmentAccessor();
		AttachmentAccessor(const AttachmentAccessor& aa);
		AttachmentAccessor(AttachmentPipe<TElem, TElemHandler>& attachmentPipe, TAttachment& attachment);

		bool access(attachment_pipe& attachmentPipe, TAttachment& attachment);

		inline typename attachment_value_traits<ValueType>::reference
		operator[](typename atraits::ConstElemPtr elem)
			{
				assert((attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem) != INVALID_ATTACHMENT_INDEX) &&
						"ERROR in AttachmentAccessor::operator[]: accessing element with invalid attachment index!");
				assert(m_pContainer && "ERROR in AttachmentAccessor::operator[]: no AttachmentPipe assigned.");
				return m_pContainer->get_elem(attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem));
			}
			
		inline typename attachment_value_traits<ValueType>::const_reference
		operator[](typename atraits::ConstElemPtr elem) const
			{
				assert((attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem) != INVALID_ATTACHMENT_INDEX) &&
						"ERROR in AttachmentAccessor::operator[]: accessing element with invalid attachment index!");
				assert(m_pContainer && "ERROR in AttachmentAccessor::operator[]: no AttachmentPipe assigned.");
				return m_pContainer->get_elem(attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem));
			}
			
/*
		inline ValueType&
		operator[](int index)
			{
				assert(m_pContainer && "ERROR in AttachmentAccessor::operator[]: no AttachmentPipe assigned.");
				return m_pContainer->get_elem(index);
			}
*/
		inline bool valid() const
			{return m_pContainer != nullptr;}

		inline void invalidate()
			{m_pContainer = nullptr;}
			
	///	calls swap on associated containers
		void swap(AttachmentAccessor<TElem, TAttachment, TElemHandler>& acc) noexcept {
			m_pContainer->swap(*acc.m_pContainer);
		}

	///	returns the raw pointer to the data of the associated container
	/**	ATTENTION: Use this method with extreme care!
	 * Returns nullptr if no container was associated.
	 */
		ValueType* raw_data()
		{
			if(m_pContainer){
				if(m_pContainer->size() > 0)
					return &(*m_pContainer)[0];
			}
			return nullptr;
		}
		
	///	returns the data index of the given element regarding the associated container.
	/**	Note that the data index does not stay constant all the time. If the
	 * associated container is, e.g., defragmented, the data index may change.*/
		size_t element_data_index(typename atraits::ConstElemPtr elem)
		{
			return attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem);
		}
		
	protected:
		ContainerType*	m_pContainer;
		TElemHandler*	m_pHandler;
};


}//	end of namespace

//	include implementation
#include "attachment_pipe.hpp"

#endif
