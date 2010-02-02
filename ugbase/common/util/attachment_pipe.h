//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d05

#ifndef __UTIL__ATTACHMENT_PIPE__
#define __UTIL__ATTACHMENT_PIPE__

#include <list>
#include <vector>
#include <stack>
#include "../static_assert.h"
#include "../types.h"
#include "uid.h"
#include "hash.h"

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
/**
*	In order to use an attachment-data-container you have to supply several Typedefs.
*	Take a look at the template-derivate AttachmentDataContainer<T> to see wich typedefs
	and operations have to be supplied in addition to the interface-methods.
*	if possible you should use the derivate-class AttachmentDataContainer<T> instead creating
	your own derivate of IAttachedDataContainer.
*/
class IAttachmentDataContainer
{
	public:
		virtual ~IAttachmentDataContainer()		{}

		virtual void resize(uint iSize) = 0;///< resize the data array
		virtual uint size() = 0;///< returns the size of the data-array.
		virtual void copy_data(uint indFrom, uint indTo) = 0;///< copy data from entry indFrom to entry indTo.
		virtual void reset_entry(uint index) = 0;///< resets the entry to its default value.
		
	///	returns a pointer to the data-buffer. Only use this method if not avoidable.
		virtual void* get_data_buffer() = 0;

	/**	copies entrys from the this-container to the buffer
	 *	specified by pDest.
	 *	For the i-th entry in the buffer, pIndexMap has to contain
	 *	the index of the source entry in this container.
	 *	num specifies the number of entries to be copied.
	 *	Nake sure, that pDest can hold 'num' elementes of the
	 *	correct size.*/
		virtual void copy_to_buffer(void* pDest, int* indexMap, int num) const = 0;
		
	/**
	*	defragment should clear the containers data from unused entries.
	*	pNewIndices should be an array of indices, wich holds a new index for each entry in IAttachedDataContainer.
	*	pNewIndices has thus to hold as many indices as there are entries in IAttachedDataContainer.
	*	If an entry shall not appear in the defragmented container, its new index has to be set to INVALID_ATTACHMENT_INDEX.
	*	numValidElements has to specify the size of the defragmented container - it thus has to equal the number of valid indices in pNewIndices.
	*/
		virtual void defragment(uint* pNewIndices, uint numValidElements) = 0;
};



/* THOUGHTS
*	AttachmentDataContainer<T> should probably be defined in another header, since it is somehow specialised for libGrid.
*	same goes for Attachment<T>
*/
////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachedDataContainer
///	A generic specialization of IAttachedDataContainer.
/**
*	This template-class not only simplifies the creation of custom containers,
	it also defines some types, operators and values, which are essential to use an AttachmentDataContainer with libGrid.
	In particular libGrids AttachmentAccessors require these definitions.
*/
template <class T> class AttachmentDataContainer : public IAttachmentDataContainer
{
	public:
		typedef T	ValueType;

		AttachmentDataContainer()	:	m_bDefaultValueSet(false)	{}
		AttachmentDataContainer(const T& defaultValue)	: m_bDefaultValueSet(true), m_defaultValue(defaultValue)	{}

		virtual ~AttachmentDataContainer()			{m_vData.clear();}

		virtual void resize(uint iSize)
			{
				if(iSize > 0)
				{
					if(m_bDefaultValueSet)
						m_vData.resize(iSize, m_defaultValue);
					else
						m_vData.resize(iSize);
				}
				else
					m_vData.clear();
			}

		virtual uint size()	{return m_vData.size();}

		virtual void copy_data(uint indFrom, uint indTo)    {m_vData[indTo] = m_vData[indFrom];}

		virtual void reset_entry(uint index)
			{
				if(m_bDefaultValueSet)
					m_vData[index] = m_defaultValue;
				else
					m_vData[index] = ValueType();
			}

		virtual void defragment(uint* pNewIndices, uint numValidElements)
			{
				uint numOldElems = size();
				std::vector<T> vDataOld = m_vData;
				m_vData.resize(numValidElements);
				for(uint i = 0; i < numOldElems; ++i)
				{
					uint nInd = pNewIndices[i];
					if(nInd != INVALID_ATTACHMENT_INDEX)
						m_vData[nInd] = vDataOld[i];
				}
			}
	
	/**	copies entrys from the this-container to the buffer
	 *	specified by pDest.
	 *	For the i-th entry in the buffer, pIndexMap has to contain
	 *	the index of the source entry in this container.
	 *	num specifies the number of entries to be copied.
	 *	Nake sure, that pDest can hold 'num' elementes of the
	 *	correct size.*/
		virtual void copy_to_buffer(void* pDest, int* indexMap, int num) const
			{
				ValueType* dest = static_cast<ValueType*>(pDest);
				for(int i = 0; i < num; ++i)
					dest[i] = m_vData[indexMap[i]];
			}


		inline const T& get_elem(uint index) const      {return m_vData[index];}
		inline T& get_elem(uint index)                  {return m_vData[index];}
		inline const T& operator[] (uint index) const	{return m_vData[index];}
		inline T& operator[] (uint index)				{return m_vData[index];}

		inline const T* get_ptr() const			{return &m_vData.front();}
		inline T* get_ptr()						{return &m_vData.front();}

		virtual void* get_data_buffer()			{return get_ptr();}

	protected:
		std::vector<T>	m_vData;
		bool			m_bDefaultValueSet;
		T				m_defaultValue;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	IAttachment
///	the interface for attachments.
/**
*	Attachments can be attached to an AttachmentPipe and thus enhance the pipes elements by data,
	whose type, container and behavior is defined by the Attachment itself.
*	In order to use an Attachment with libGrid (in particular with libGrids AttachmentAccessors),
	derivatives of IAttachment have to feature some special typedefs (see Attachment<T> for more information).
*	Whenever possible you should use the template-derivative Attachment<T> instead of IAttachment.
*/
class IAttachment : public UID
{
	public:
		IAttachment() : m_name("")   {}
		IAttachment(const char* name) : m_name(name)   {}

		virtual ~IAttachment()  {}
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
/**
*	This class is intended to simplify the process of Attachment creation.
*	Note that there are typedefs, which are required by libGrids AttachmentAccessors.
*/
template <class T> class Attachment : public IAttachment
{
	public:
		typedef AttachmentDataContainer<T>	ContainerType;
		typedef T							ValueType;

		Attachment() : IAttachment(), m_passOnBehaviour(false)    {}
		Attachment(bool passOnBehaviour) : IAttachment(), m_passOnBehaviour(passOnBehaviour)    {}
		Attachment(const char* name) : IAttachment(name), m_passOnBehaviour(false)   		{}
		Attachment(const char* name, bool passOnBehaviour) : IAttachment(name), m_passOnBehaviour(passOnBehaviour)	{}

		virtual ~Attachment()	{}
		virtual IAttachment* clone()							{IAttachment* pA = new Attachment<T>; *pA = *this; return pA;}
		virtual IAttachmentDataContainer* create_container()	{return new ContainerType;}
		virtual bool default_pass_on_behaviour() const			{return m_passOnBehaviour;}
		IAttachmentDataContainer* create_container(const T& defaultValue)	{return new ContainerType(defaultValue);}

	protected:
		bool m_passOnBehaviour;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachmentEntry
///	This struct is used by AttachmentPipe in order to manage its attachments
struct AttachmentEntry
{
	AttachmentEntry() : m_pAttachment(NULL), m_pContainer(NULL), m_userData(0)	{}
	AttachmentEntry(IAttachment* pAttachment, IAttachmentDataContainer* pContainer, uint userData = 0) :
			m_pAttachment(pAttachment), m_pContainer(pContainer), m_userData(userData)	{}

	IAttachment*				m_pAttachment;
	IAttachmentDataContainer*	m_pContainer;
	uint						m_userData;
};


////////////////////////////////////////////////////////////////////////////////////////////////
///	attachment_traits define the interface that enables you to use your own data-types with the AttachmentPipe.
/**
 * Perform template-specialization for your own data-types and their respective handlers.
 */
template<class TElem, class TElemHandler>
class attachment_traits
{
	public:
		typedef TElem&		ElemRef;
		typedef void*		ElemPtr;
		typedef const void* ConstElemPtr;
		typedef void*		ElemHandlerPtr;
		typedef const void*	ConstElemHandlerPtr;

	///	mark the element as invalid.
	/**	You may do something like elem = NULL, if TElem is a pointer type.*/
		static inline void invalidate_entry(ElemHandlerPtr pHandler, ElemRef elem)				{/*STATIC_ASSERT(0, INVALID_ATTACHMENT_TRAITS);*/}
		static inline bool entry_is_invalid(ElemHandlerPtr pHandler, ElemRef elem)				{return true;/*STATIC_ASSERT(0, INVALID_ATTACHMENT_TRAITS);*/}
		static inline uint get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem)			{return INVALID_ATTACHMENT_INDEX;/*STATIC_ASSERT(0, INVALID_ATTACHMENT_TRAITS);*/}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index)	{/*STATIC_ASSERT(0, INVALID_ATTACHMENT_TRAITS);*/}
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachmentPipe
///	Handles data which has been attached to the pipe using callbacks for the element.
/**
 * The AttachmentPipe can be used to attach data to a collection of elements.
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
class AttachmentPipe
{
	public:
		typedef std::list<AttachmentEntry>			AttachmentEntryContainer;
		typedef AttachmentEntryContainer::iterator	AttachmentEntryIterator;
		typedef AttachmentEntryContainer::const_iterator	ConstAttachmentEntryIterator;
		typedef Hash<AttachmentEntryIterator, uint>		AttachmentEntryIteratorHash;
		typedef attachment_traits<TElem, TElemHandler>	atraits;

	public:
		AttachmentPipe() : m_pHandler(NULL)	{}
		AttachmentPipe(typename atraits::ElemHandlerPtr pHandler) : m_pHandler(pHandler)	{}
		
		inline typename atraits::ElemHandlerPtr get_elem_handler()		{return m_pHandler;}
		
		void clear();
		void clear_elements();
		void clear_attachments();

		void register_element(TElem elem);
		void unregister_element(const TElem& elem);

		void defragment();

		template <class TAttachment>
		void attach(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, uint options);

		void attach(IAttachment& attachment, uint options = 0);
		void detach(IAttachment& attachment);
		bool has_attachment(IAttachment& attachment) const;

		template <class TAttachment>
		typename TAttachment::ValueType*
		get_data_array(TAttachment& attachment);
		
		IAttachmentDataContainer* get_data_container(IAttachment& attachment) const;
		
		template <class TAttachment>
		typename TAttachment::ContainerType*
		get_data_container(TAttachment& attachment);

		inline ConstAttachmentEntryIterator attachments_begin() const	{return m_attachmentEntryContainer.begin();}
		inline ConstAttachmentEntryIterator attachments_end() const		{return m_attachmentEntryContainer.end();}

		inline uint num_elements() const		{return m_numElements;}
		inline uint num_data_entries() const	{return m_vEntries.size();}
		inline bool is_fragmented() const		{return m_numElements != m_vEntries.size();}

		void reset_values(uint dataIndex);///< fills the attached data at the given index with the default values.

	protected:
		void resize_attachment_containers(uint newSize);
		void grow_attachment_containers(uint newMinSize);
		uint get_container_size();

	protected:
		typedef std::vector<TElem>		ElemEntryVec;
		typedef std::stack<uint>		UINTStack;

	protected:
		AttachmentEntryContainer	m_attachmentEntryContainer;
		AttachmentEntryIteratorHash	m_attachmentEntryIteratorHash;

		ElemEntryVec	m_vEntries;	///< same size as attachment containers.
		UINTStack		m_stackFreeEntries;	///< holds indices to free entries.

		uint			m_numElements;
		
		typename atraits::ElemHandlerPtr	m_pHandler;
};


////////////////////////////////////////////////////////////////////////////////////////////////
//	AttachmentAccessor
///	Used to access data that has been attached to an attachment pipe.
/**
 * Once initialized (use the constructor), an AttachmentAccessor can be used to access
 * the data stored in the given AttachmentPipe
 */
template <class TElem, class TAttachment, class TElemHandler>
class AttachmentAccessor
{
	public:
		typedef typename TAttachment::ValueType		ValueType;
		typedef typename TAttachment::ContainerType			ContainerType;
		typedef attachment_traits<TElem, TElemHandler>	atraits;

	public:
		AttachmentAccessor();
		AttachmentAccessor(const AttachmentAccessor& aa);
		AttachmentAccessor(AttachmentPipe<TElem, TElemHandler>& attachmentPipe, TAttachment& attachment);

		void access(AttachmentPipe<TElem, TElemHandler>& attachmentPipe, TAttachment& attachment);

		inline ValueType&
		operator[](typename atraits::ConstElemPtr elem)
			{
				assert((attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem) != INVALID_ATTACHMENT_INDEX) &&
						"ERROR in AttachmentAccessor::operator[]: accessing element with invalid attachment index!");
				assert(m_pContainer && "ERROR in AttachmentAccessor::operator[]: no AttachmentPipe assigned.");
				return m_pContainer->get_elem(attachment_traits<TElem, TElemHandler>::get_data_index(m_pHandler, elem));
			}
			
		inline const ValueType&
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
		inline bool valid()
			{return m_pContainer != NULL;}

		inline void invalidate()
			{m_pContainer = NULL;}

	protected:
		ContainerType*	m_pContainer;
		TElemHandler*	m_pHandler;
};


}//	end of namespace

//	include implementation
#include "attachment_pipe.hpp"

#endif
