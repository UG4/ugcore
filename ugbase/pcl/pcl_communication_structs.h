//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_COMMUNICATION_STRUCTS__
#define __H__PCL__PCL_COMMUNICATION_STRUCTS__

#include <iterator>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include "common/util/smart_pointer.h"
#include "common/util/binary_buffer.h"
#include "common/error.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
//	type-traits
///	associate internally used types with an external typename
/**
 * By default it is assumed, that the external typename and the element
 * type used in interfaces match.
 *
 * You may specialize type_traits for your own types through template
 * specialization. Be sure to specialize them before you use any
 * pcl-classes with your type.
 */
template <class TType>
struct type_traits
{
	typedef TType Elem;	///	Type of interface elements
};


////////////////////////////////////////////////////////////////////////
//	Interface Categories
///	Interface tags allow to differentiate between interfaces with different features
namespace interface_tags
{
/**	Every interface has to feature this tag (at least indirectly - by tag-derivation).
 *	This tag simply says that one may iterate over the elements of the
 *	interface.
 *
 *	Additionally you may receive the local_src_id. This id represents the
 *	sender during local communication (communication on one process only).
 *	Values < 0 mark the sender as invalid (during local communication).
 *	The src-id is ignored during parallel communication. Instead pcl::GetProcRank
 *	is used.
 *
 *	typedefs that have to be featured in such interfaces:
 *		- category_tag
 *		- iterator
 *		- const_iterator
 *		- Type
 *		- Element
 *
 *	Methods that have to be featured in such interfaces:
 *		- iterator begin();
 *		- iterator end();
 *		- const_iterator begin() const;
 *		- const_iterator end() const;
 *		- Element& get_element(iterator iter);
 *		- const Element& get_element(iterator iter) const;
 *		- int get_target_proc();
 *		- void swap(Interface& interface);
 */
class basic_interface_tag									{};

/**	The ordered_interface_tag derives from the basic_interface_tag.
 *	Thus all classes and methods that may operate on interfaces
 *	with the basic_interface_tag will operate on interfaces with
 *	this tag, too.
 *
 *	Interfaces with this category have to associated a local id with
 *	each entry.
 *
 *	Methods that have to be featured in such interfaces:
 *		- get_local_id(iterator iter);
 */
class ordered_interface_tag : public basic_interface_tag	{};
}//	end of namespace interface_tags


////////////////////////////////////////////////////////////////////////
//	BasicInterface
///	You may add elements to this interface and iterate over them
/**
 * This interface type features a minimal set of methods that is
 * required to actually work with it.
 * You may add new elements, erase old ones and iterate through them.
 *
 * In order to access the associated element of an iterator you should
 * use get_element (not the * operator of the iterator). This increases
 * the flexibility of your code.
 *
 * You may specify a stl-container compatible type that is used to store.
 * the elements. Supported types are std::vector and std::list.
 *
 * TContainer defaults to std::vector. This is the best option if your
 * interface is mainly used static or grows, but is considerably slower
 * than std::list if you want to erase elements often.
 * For dynamic interfaces std::list may be the better option.
 */
template <class TType,
		  template <class T, class Alloc> class TContainer = std::vector,
		  template <class T> class TAlloc = std::allocator>
class BasicInterface
{
	protected:
		typedef typename type_traits<TType>::Elem		TElem;
		typedef TContainer<TElem, TAlloc<TElem> >		ElemContainer;

	public:
		typedef interface_tags::basic_interface_tag		category_tag;

		typedef TType									Type;
		typedef typename type_traits<TType>::Elem		Element;
		typedef typename ElemContainer::iterator		iterator;
		typedef typename ElemContainer::const_iterator	const_iterator;

	public:
		BasicInterface(int targetProc = -1)	: m_size(0), m_targetProc(targetProc)	{};

		inline iterator push_back(const Element& elem)	{++m_size; return m_elements.insert(m_elements.end(), elem);}
		inline iterator erase(iterator iter)				{--m_size; return m_elements.erase(iter);}

		inline iterator begin()		{return m_elements.begin();}
		inline iterator end()		{return m_elements.end();}
		inline const_iterator begin() const		{return m_elements.begin();}
		inline const_iterator end() const		{return m_elements.end();}

		inline Element& get_element(iterator iter)						{return *iter;}
		inline const Element& get_element(const_iterator iter) const	{return *iter;}

	///	returns the number of elements that are stored in the interface.
		inline size_t size() const					{return m_size;}
		inline bool empty() const					{return size() == 0;}

		int get_target_proc() const					{return m_targetProc;}

	///	swaps the content of two interfaces.
	/** m_elements, m_size and m_targetProc are swapped.*/
		void swap(BasicInterface& interface)
		{
			using std::swap;
			m_elements.swap(interface.m_elements);
			swap(m_size, interface.m_size);
			swap(m_targetProc, interface.m_targetProc);
		}

	///	sort the entries in this interface.
		template <class TCompare>
		void sort(TCompare cmp)
		{
			using std::sort;
			m_elements.sort(cmp);
		}

	protected:
		ElemContainer	m_elements;
		size_t			m_size;
		int				m_targetProc;
};

////////////////////////////////////////////////////////////////////////
//	OrderedInterface
///	You may add elements to this interface and iterate over them
/**	You may retrieve a localID for each element in the interface*/
template <class TType,
		  template <class T, class Alloc> class TContainer = std::vector,
		  template <class T> class TAlloc = std::allocator>
class OrderedInterface
{
	protected:
		typedef typename type_traits<TType>::Elem	TElem;

		struct InterfaceEntry
		{
			InterfaceEntry(TElem e, size_t locID) : elem(e), localID(locID)	{}

			TElem	elem;
			size_t	localID;
		};

		template <class TElemCmp>
		struct InterfaceEntryCmp{
			InterfaceEntryCmp(TElemCmp elemCmp) : m_elemCmp(elemCmp) {}
			bool operator()(InterfaceEntry const& e1, InterfaceEntry const& e2)
			{return m_elemCmp(e1.elem, e2.elem);}
			TElemCmp	m_elemCmp;
		};

		typedef TContainer<InterfaceEntry, TAlloc<InterfaceEntry> >
														ElemContainer;

	public:
		typedef TType									Type;
		typedef typename type_traits<TType>::Elem		Element;

		typedef interface_tags::ordered_interface_tag	category_tag;

		typedef typename ElemContainer::iterator		iterator;
		typedef typename ElemContainer::const_iterator	const_iterator;

	public:
		OrderedInterface(int targetProc = -1) :
			m_size(0),
			m_targetProc(targetProc),
			m_idCounter(1)	{}

		inline iterator push_back(const Element& elem)
		{
			++m_size;
			return m_elements.insert(m_elements.end(),
									(InterfaceEntry(elem, get_free_id())));
		}

		inline iterator erase(iterator iter)
		{
			--m_size;
			return m_elements.erase(iter);
		}

		inline iterator begin()		{return m_elements.begin();}
		inline iterator end()		{return m_elements.end();}

		inline const_iterator begin() const {return m_elements.begin();}
		inline const_iterator end() const {return m_elements.end();}

		inline Element& get_element(iterator iter)	{return (*iter).elem;}
		inline size_t get_local_id(iterator iter)	{return (*iter).localID;}

		inline const Element& get_element(const_iterator iter) const	{return (*iter).elem;}
		inline size_t get_local_id(const_iterator iter)	const 			{return (*iter).localID;}


	///	returns the number of elements that are stored in the interface.
		inline size_t size() const					{return m_size;}
		inline bool empty()	const					{return size() == 0;}

		int get_target_proc() const					{return m_targetProc;}

	///	returns true if iter1 < iter2.
		static inline bool cmp(iterator iter1, iterator iter2,
								const std::input_iterator_tag&)
		{
			return (*iter1).localID < (*iter2).localID;
		}

	///	swaps the content of two interfaces.
	/** m_elements, m_size, m_targetProc and m_idCounter are swapped.*/
		void swap(OrderedInterface& interface)
		{
			using std::swap;
			m_elements.swap(interface.m_elements);
			swap(m_size, interface.m_size);
			swap(m_targetProc, interface.m_targetProc);
			swap(m_idCounter, interface.m_idCounter);
		}
		
	///	sort the entries in this interface.
		template <class TCompare>
		void sort(TCompare cmp)
		{
			using std::sort;
			InterfaceEntryCmp<TCompare> ieCmp(cmp);
			m_elements.sort(ieCmp);

		//	we have to reset all local indices
			m_idCounter = 1;
			for(iterator iter = m_elements.begin(); iter != m_elements.end(); ++iter)
				(*iter).localID = m_idCounter++;
		}


	protected:
	///	returns a free id in each call.
	/**	Those ids are not necessarily aligned.*/
		size_t get_free_id()
		{
			if(m_idCounter == 0)
			{
				m_idCounter = 1;
			//	we have to reset all entries
				for(iterator iter = begin(); iter != end(); ++iter)
					(*iter).localID = m_idCounter++;
			}

			return m_idCounter++;
		}

	protected:
		ElemContainer	m_elements;
		size_t m_size;
		int m_targetProc;
		size_t m_idCounter;
};


////////////////////////////////////////////////////////////////////////
//	Layout Categories
///	Layout tags allow to differentiate between layouts with different features
namespace layout_tags
{
///	marks a layout as a single-level layout
/**
 * Typedefs that have to be supported by a layout with this tag:
 *		- category_tag
 *		- iterator
 * 		- Interface
 *		- Type
 *		- Element
 *
 * Methods that have to be supported by a layout with this tag:
 *		- iterator begin()
 *		- iterator end()
 *		- Interface& interface(iterator iter)
 *		- int proc_id(iterator iter)
 *		- iterator erase(iterator iter)
 *		- int get_local_src_id()
 */
class single_level_layout_tag	{};

///	marks a layout as a multi-level layout
/**
 * Typedefs that have to be supported by a layout with this tag:
 *		- category_tag
 *		- iterator
 * 		- Interface
 *		- Element
 *
 * Methods that have to be supported by a layout with this tag:
 *		- size_t num_levels()
 *		- iterator begin(size_t level)
 *		- iterator end(size_t level)
 *		- Interface& interface(iterator iter)
 *		- int proc_id(iterator iter)
 *		- iterator erase(iterator iter, size_t level)
 *		- int get_local_src_id()
 */
class multi_level_layout_tag	{};

}//	end of namespace layout_tags


////////////////////////////////////////////////////////////////////////
//	SingleLevelLayout
///	the standard single-level-layout implementation
/**
 * A Layout is a collection of interfaces.
 * Each interface is associated with a process-id.
 *
 * This layout type supports the requirements of the
 * pcl::layout_tags::single_level_layout_tag category.
 *
 * The layout passes its local srcID on to created interfaces.
 *
 * Additionally it features methods that allow to add new interfaces
 *
 * In order to allow one method to operate both on a SingleLevelLayout
 * and a MultiLevelLayout, the (size_t level = 0) convenience parameter
 * has been added to some methods. Those parameters are ignored throughout
 * the whole implementation.
 */
template <class TInterface>
class SingleLevelLayout
{
	public:
		SingleLevelLayout()		{}

	////////////////////////////////////////////////
	//	typedefs required by implementation
	///	an interface-map is a list of interfaces, each associated with a process id.
		typedef std::map<int, TInterface>	InterfaceMap;

	////////////////////////////////////////////////
	//	typedefs required by layout-tag
	///	Layout category
		typedef layout_tags::single_level_layout_tag	category_tag;

	///	Interface type
		typedef TInterface						Interface;

	///	Type
		typedef typename Interface::Type		Type;

	///	Element type
		typedef typename Interface::Element		Element;

	///	An iterator that allows to iterate over the interfaces stored in the layout.
		typedef typename InterfaceMap::iterator			iterator;
		typedef typename InterfaceMap::const_iterator	const_iterator;

	public:
	////////////////////////////////////////////////
	//	methods required by the layout-tag

	///	returns the iterator to the first interface of the layout.
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator begin(size_t level = 0)					{return m_interfaceMap.begin();}
		inline const_iterator begin(size_t level = 0) const		{return m_interfaceMap.begin();}

	///	returns the iterator to the last interface of the layout.
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator end(size_t level = 0)					{return m_interfaceMap.end();}
		inline const_iterator end(size_t level = 0)	const		{return m_interfaceMap.end();}

	///	returns true if the layout has no interfaces.
	/**	Note that this method only tells whether there are interfaces or not.
	 * To check whether there are any elements use has_interface_elements.
	 */
		inline bool empty(size_t level = 0)	const 				{return begin() == end();}
		
	///	returns 1
		inline size_t num_levels() const						{return 1;}
		
	///	returns the interface to the given iterator.
		inline Interface& interface(iterator iter)						{return iter->second;}
		inline const Interface& interface(const_iterator iter) const	{return iter->second;}

	///	returns the target process of the interface given in iterator
		inline int proc_id(iterator iter) const 				{return iter->first;}
		inline int proc_id(const_iterator iter) const			{return iter->first;}

	///	erases the interface at the given iterator.
	/**	returns an iterator to the next interface.*/
		inline iterator erase(iterator iter, size_t level = 0)
			{
				iterator tIter = iter++;
				m_interfaceMap.erase(tIter);
				return iter;
			}
		
	///	clears the layout
		void clear()											{m_interfaceMap.clear();}

	///	returns the interface to the given process.
	/**	if the queried interface exist, it will be returned.
	 *	If not it will be created.
	 *	The new interfaces localSrcID will be set to the localSrcID of this layout.*/
		inline Interface& interface(int procID, size_t level = 0)
		{
			iterator iter = m_interfaceMap.find(procID);
			if(iter != m_interfaceMap.end())
				return iter->second;
			return m_interfaceMap.insert(make_pair(procID, Interface(procID))).first->second;
		}

		inline const Interface& interface(int procID, size_t level = 0) const
		{
			const_iterator iter = m_interfaceMap.find(procID);
			UG_ASSERT(iter != m_interfaceMap.end(), "trying to access an non-existing interface in a constant layout");
			return iter->second;
		}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID, size_t level = 0) const
			{return m_interfaceMap.find(procID) != m_interfaceMap.end();}

	///	returns the sum of the interface sizes
		inline size_t num_interface_elements() const
		{
			size_t sum = 0;
			for(const_iterator iter = begin(); iter != end(); ++iter)
				sum += interface(iter).size();
			return sum;
		}

	///	returns true if the layout contains interface elements
		inline bool has_interface_elements() const
		{
			for(iterator iter = begin(); iter != end(); ++iter){
				if(!interface(iter).empty())
					return true;
			}
			return false;
		}

	/// returns the number of interfaces in the layout
		inline size_t num_interfaces(size_t level = 0) const
		{
			return m_interfaceMap.size();
		}

	///	sort the entries in all interfaces of this layout
		template <class TCompare>
		void sort_interface_entries(TCompare cmp)
		{
			for(iterator iter = begin(); iter != end(); ++iter)
				interface(iter).sort(cmp);
		}

	private:
/*
	///	copy-constructor is not yet implemented
		SingleLevelLayout(const SingleLevelLayout& sll);

	///	assignement-operator is not yet implemented
		SingleLevelLayout& operator = (const SingleLevelLayout& sll);
*/
	protected:
	///	holds the interfaces in a map.
		InterfaceMap	m_interfaceMap;
};

////////////////////////////////////////////////////////////////////////
//	MultiLevelLayout
///	the standard multi-level-layout implementation
/**
 * A MultiLevelLayout is a collection of interfaces, which are
 * grouped in different levels.
 *
 * Each interface is associated with a process-id.
 * If the localSrcID of this layout == -1 then, when a new level is
 * created, it's local srcID is automatically
 * initialized with the level index. localSrcIDs of interfaces created on
 * those levels are then initialised with the level index, too.
 *
 * If the localSrcID >= 0 it is simply passed on to the level-layouts.
 *
 * This layout type supports the requirements of the
 * pcl::layout_tags::multi_level_layout_tag category.
 *
 * Additionally it features methods that allow to add new interfaces.
 */
template <class TInterface>
class MultiLevelLayout
{
	public:
	////////////////////////////////////////////////
	//	typedefs required by implementation
	///	on each level a single-level-layout is used
		typedef SingleLevelLayout<TInterface>		LevelLayout;

	////////////////////////////////////////////////
	//	typedefs required by layout-tag
	///	Layout category
		typedef layout_tags::multi_level_layout_tag	category_tag;

	///	Interface type
		typedef TInterface							Interface;

	///	Type
		typedef typename Interface::Type			Type;

	///	Element type
		typedef typename Interface::Element			Element;

	///	An iterator that allows to iterate over the interfaces stored in the layout.
		typedef typename LevelLayout::iterator			iterator;
		typedef typename LevelLayout::const_iterator	const_iterator;

	public:
		MultiLevelLayout()								{}
		MultiLevelLayout(const MultiLevelLayout& mll)	{assign_layout(mll);}

		~MultiLevelLayout()								{clear();}

		MultiLevelLayout& operator = (const MultiLevelLayout& mll)
		{
			assign_layout(mll);
			return *this;
		}

	////////////////////////////////////////////////
	//	methods required by the layout-tag

	///	returns the iterator to the first interface of the layout in the given level.
	/**	You should access the values of this iterator using the methods
	 *	Layout::interface and Layout::proc_id.
	 *	Make sure that the level matches the level in the associated end() call.*/
		inline iterator begin(size_t level)				{require_level(level); return m_vLayouts[level]->begin();}
		inline const_iterator begin(size_t level) const	{require_level(level); return m_vLayouts[level]->begin();}

	///	returns the iterator to the last interface of the layout in the given level.
	/**	You should access the values of this iterator using the methods
	 *	Layout::interface and Layout::proc_id.
	 *	Make sure that the level matches the level in the associated begin() call.*/
		inline iterator end(size_t level)				{require_level(level); return m_vLayouts[level]->end();}
		inline const_iterator end(size_t level) const	{require_level(level); return m_vLayouts[level]->end();}

	///	returns true if the layout has no interfaces on the given level.
	/**	Note that this method only tells whether there are interfaces or not.
	 * To check whether there are any elements use has_interface_elements.
	 */
		inline bool empty(size_t level)					{return begin(level) == end(level);}
		inline bool empty(size_t level) const			{return begin(level) == end(level);}

	///	returns true if the layout has no interfaces.
	/**	Note that this method only tells whether there are interfaces or not.
	 * To check whether there are any elements use has_interface_elements.
	 */
		inline bool empty()	const						{for(size_t l = 0; l < num_levels(); ++l){if(!empty(l)) return false;} return true;}

	///	returns the interface to the given iterator.
		inline Interface& interface(iterator iter)						{return iter->second;}
		inline const Interface& interface(const_iterator iter) const	{return iter->second;}

	///	returns the interface to the given iterator.
		inline int proc_id(const_iterator iter) const					{return iter->first;}

	///	erases the interface at the given iterator on the given level.
	/**	returns an iterator to the next interface.*/
		inline iterator erase(iterator iter, size_t level)	{return m_vLayouts[level]->erase(iter);}
		
	///	returns the number of levels.
		inline size_t num_levels()	const					{return m_vLayouts.size();}

	////////////////////////////////////////////////
	//	methods that enhance the layout-tag
	///	deletes all levels.
		void clear()
		{
			for(size_t i = 0; i < m_vLayouts.size(); ++i)
				delete m_vLayouts[i];
			m_vLayouts.clear();
		}

	///	returns the interface to the given process.
	/**	if the interface doesn't exist yet, it will be created.*/
		inline Interface& interface(int procID, size_t level)				{require_level(level); return m_vLayouts[level]->interface(procID);}
		inline const Interface& interface(int procID, size_t level) const	{require_level(level); return m_vLayouts[level]->interface(procID);}

	///	returns true if an interface to the given procID on the given level already exists.
		inline bool interface_exists(int procID, size_t level)			{require_level(level); return m_vLayouts[level]->interface_exists(procID);}
		inline bool interface_exists(int procID, size_t level) const	{require_level(level); return m_vLayouts[level]->interface_exists(procID);}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID)				{for(size_t i = 0; i < num_levels(); ++i){if(m_vLayouts[i]->interface_exists(procID)) return true;} return false;}
		inline bool interface_exists(int procID) const			{for(size_t i = 0; i < num_levels(); ++i){if(m_vLayouts[i]->interface_exists(procID)) return true;} return false;}

	///	returns the layout at the given level.
	/**	If level >= num_levels() then the layouts in between
		will be automatically created.*/
		inline LevelLayout& layout_on_level(int level)				{require_level(level); return *m_vLayouts[level];}
		inline const LevelLayout& layout_on_level(int level) const	{require_level(level); return *m_vLayouts[level];}

	///	returns the sum of the interface sizes
		inline size_t num_interface_elements() const
		{
			size_t sum = 0;
			for(size_t lvl = 0; lvl < num_levels(); ++lvl){
				for(iterator iter = begin(lvl); iter != end(lvl); ++iter)
					sum += interface(iter).size();
			}
			return sum;
		}

	///	returns true if the layout contains any interface entries
		inline bool has_interface_elements() const
		{
			for(size_t lvl = 0; lvl < num_levels(); ++lvl){
				for(iterator iter = begin(lvl); iter != end(lvl); ++iter){
					if(!interface(iter).empty())
						return true;
				}
			}
			return false;
		}

	/// returns the number of interfaces in the layout
		inline size_t num_interfaces(size_t level) const
		{
			return m_vLayouts[level].size();
		}

	///	sort the entries in all interfaces of this layout
		template <class TCompare>
		void sort_interface_entries(TCompare cmp)
		{
			for(size_t lvl = 0; lvl < num_levels(); ++lvl)
				layout_on_level(lvl).sort_interface_entries(cmp);
		}

	protected:
	///	adds num new levels.
		inline void new_levels(size_t num)		{for(size_t i = 0; i < num; ++i) m_vLayouts.push_back(new LevelLayout());}

	///	if the required level doesn't exist yet, it will created.
		inline void require_level(size_t level)			{if(level >= num_levels()) new_levels(level - num_levels() + 1);}
		inline void require_level(size_t level) const	{if(level >= num_levels()){UG_THROW("Level too high: " << level << ", while num_levels == " << num_levels());}}

	///	clears this layout and then copies all levels from the given layout
		void assign_layout(const MultiLevelLayout& mll)
		{
			clear();
			for(size_t i = 0; i < mll.m_vLayouts.size(); ++i)
				m_vLayouts.push_back(new LevelLayout(*mll.m_vLayouts[i]));
		}

	protected:
		std::vector<LevelLayout*>	m_vLayouts;
};

////////////////////////////////////////////////////////////////////////
//	ICommunicationPolicy
///	specializations are responsible to pack and unpack interface data during communication.
/**	Make sure that you use the same communication-policy for send and receive operations.
 *	Otherwise problems regarding buffer-sizes may occur.
 */
template <class TLayout>
class ICommunicationPolicy
{
	public:
		typedef TLayout						Layout;
		typedef typename Layout::Interface 	Interface;

		virtual ~ICommunicationPolicy()	{}

	////////////////////////////////
	//	COLLECT AND EXTRACT
	///	returns the size of the buffer in bytes, that will be required for interface-communication.
	/**	Determines the size of the buffer on which the extract and receive methods
	 *	for the given interface will operate.
	 *	If the buffer-size can't be calculated on both sides (sender and receiver)
	 *	this method should return -1. This will lead to an additional communication
	 *	step in which buffer-sizes will be exchanged.
	 *	If the buffer-size can be calculated on both sides, it makes sense to do so,
	 *	since this leads to less communication and overall improved performance.
	 *	The buffer-size has to exactly match the size of required memory. Make sure that you
	 *	completely fill the buffer during collect(...) and that you read all data during
	 *	extract(...).
	 *	The default implementation returns -1.
	 */	
		virtual int
		get_required_buffer_size(const Interface& interface)	{return -1;}
		
	////////////////////////////////
	//	COLLECT
	///	signals the beginning of a layout collection.
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		begin_layout_collection(const Layout* pLayout)			{return true;}

	///	signals the end of a layout collection
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		end_layout_collection()									{return true;}

	///	should write data which is associated with the interface elements to the buffer.
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface) = 0;

	////////////////////////////////
	//	EXTRACT
	///	signals the beginning of a layout extraction.
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		begin_layout_extraction(const Layout* pLayout)			{return true;}

	///	signals the end of a layout extraction
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		end_layout_extraction()									{return true;}

	///	signals that a new layout-level will now be processed.
	/**	This is primarily interesting for layout-extraction of multi-level-layouts.
	 *	Before extract is called for the interfaces of one level of a layout,
	 *	begin_level_extraction(level) is called.
	 *	If single-level-layouts are processed, this method is called
	 *	once with level = 0.
	 *	This method is called after begin_layout_extraction and before
	 *	the associated extract calls.*/
	 	virtual void begin_level_extraction(int level)			{}

	///	extract data from the buffer and assigns it to the interface-elements.
	/**	If this method is called between calls to begin_layout_extraction and
		end_layout_extraction, the interface that is passed to this method
		belongs to the layout.*/
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface) = 0;
};

}//	end of namespace pcl

#endif
