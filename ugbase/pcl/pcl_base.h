//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_BASE__
#define __H__PCL__PCL_BASE__

#include <iterator>
#include <vector>
#include <list>
#include <map>

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
 */
class basic_interface_tag									{};

/**	The ordered_interface_tag derives from the basic_interface_tag.
 *	Thus all classes and methods that may operate on interfaces
 *	with the basic_interface_tag will operate on interfaces with
 *	this tag, too.
 *
 *	Interfaces with this category have to feature a cmp compare
 *	method that compares two iterators. If iter1 < iter2 then
 *	the method should return true
 *
 *	Methods that have to be featured in such interfaces:
 *		- static bool cmp(iterator iter1, iterator iter2);
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
		BasicInterface()	: m_size(0)		{};
		inline iterator push_back(const Element& elem)	{++m_size; return m_elements.insert(m_elements.end(), elem);}
		inline void erase(iterator iter)				{--m_size; m_elements.erase(iter);}
		
		inline iterator begin()		{return m_elements.begin();}
		inline iterator end()		{return m_elements.end();}
		
		inline Element& get_element(iterator iter)	{return *iter;}
		
	///	returns the number of elements that are stored in the interface.
		inline size_t size()						{return m_size;}
		
	protected:
		ElemContainer	m_elements;
		size_t			m_size;
};


///	Internally used by pcl::OrderedInterface
/**	this struct holds an elem and a local id, as they are
 *	required if we use non-random-access-iterators*/
template <typename TIterTag, class TElem>
struct TOrderedInterfaceEntryByTag
{
	TOrderedInterfaceEntryByTag(const TElem& e, uint ID) :
		elem(e), localID(ID)	{}
		
	TElem	elem;
	uint	localID;
};

///	Internally used by pcl::OrderedInterface
/**	specialization for random-access-iterators.
 *	A local id is not required here.*/
template <class TElem> struct
TOrderedInterfaceEntryByTag<std::random_access_iterator_tag, TElem>
{
	TOrderedInterfaceEntryByTag(const TElem& e) : elem(e)	{}
	
	TElem elem;
};


////////////////////////////////////////////////////////////////////////
//	OrderedInterface
///	You may add elements to this interface and iterate over them
/**	You may compare iterators of elements in this interface using the
 *	cmp method.
 *
 *	This interface stores elements in a vector internally.
 *	This allows for easy compare but makes erase expensive.*/
template <class TType,
		  template <class T, class Alloc> class TContainer = std::vector,
		  template <class T> class TAlloc = std::allocator>
class OrderedInterface
{
	protected:
		typedef typename type_traits<TType>::Elem	TElem;
		
	//	this is required to pick the right InterfaceEntryType
		typedef typename
			std::iterator_traits<
				typename TContainer<TElem, TAlloc<TElem> >::iterator>
				::iterator_category						iterator_category;
			
		typedef TOrderedInterfaceEntryByTag<iterator_category, TElem>
														InterfaceEntry;
														
		typedef TContainer<InterfaceEntry, TAlloc<InterfaceEntry> >
														ElemContainer;
		
	public:
		typedef TType									Type;
		typedef typename type_traits<TType>::Elem		Element;
		
		typedef interface_tags::ordered_interface_tag	category_tag;
		
		typedef typename ElemContainer::iterator		iterator;
		typedef typename ElemContainer::const_iterator	const_iterator;
		
	public:
		OrderedInterface() : m_size(0), m_idCounter(1)	{}
		
		inline iterator push_back(const Element& elem)
		{
			return push_back(elem, iterator_category());
		}
		
		inline void erase(iterator iter)
		{
			--m_size;
			m_elements.erase(iter);
		}
		
		inline iterator begin()		{return m_elements.begin();}
		inline iterator end()		{return m_elements.end();}
		
		inline Element& get_element(iterator iter)	{return (*iter).elem;}
		
	///	returns the number of elements that are stored in the interface.
		inline size_t size()						{return m_size;}
		
	///	returns true if iter1 < iter2.
		static inline bool cmp(iterator iter1, iterator iter2)
		{
			return cmp(iter1, iter2, iterator_category());
		}
		
	protected:
		inline iterator push_back(const Element& elem, 
								const std::input_iterator_tag&)
		{
			++m_size;
			return m_elements.insert(m_elements.end(),
									(InterfaceEntry(elem, get_free_id())));
		}
		
		inline iterator push_back(const Element& elem, 
								const std::random_access_iterator_tag&)
		{
			++m_size;
			return m_elements.insert(m_elements.end(), elem);
		}
		
		static inline bool cmp(iterator iter1, iterator iter2,
								const std::input_iterator_tag&)
		{
			return (*iter1).localID < (*iter2).localID;
		}

		static inline bool cmp(iterator iter1, iterator iter2,
								const std::random_access_iterator_tag&)
		{
			return iter2 - iter1 > 0;
		}
				
	///	returns a free id in each call.
	/**	This method may only be called if
	 *	iterator_category != std::random_access_iterator_tag.*/
		uint get_free_id()
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
		uint m_size;
		uint m_idCounter;
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
 *		- Interface& interface(iterator& iter)
 *		- int proc_id(iterator& iter)
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
 *		- Interface& interface(iterator& iter)
 *		- int proc_id(iterator& iter)
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
 * Additionally it features methods that allow to add new interfaces
 */
template <class TInterface>
class SingleLevelLayout
{
	public:
		SingleLevelLayout()	{}
		
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
		typedef typename InterfaceMap::iterator		iterator;

	public:
	////////////////////////////////////////////////
	//	methods required by the layout-tag

	///	returns the iterator to the first interface of the layout.
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator begin()							{return m_interfaceMap.begin();}

	///	returns the iterator to the last interface of the layout.	
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator end()							{return m_interfaceMap.end();}

	///	returns the interface to the given iterator.
		inline Interface& interface(iterator& iter)		{return iter->second;}
		
	///	returns the interface to the given iterator.
		inline int proc_id(iterator& iter)				{return iter->first;}
	
	
	////////////////////////////////////////////////
	//	methods that enhance the layout-tag
	
	///	returns the interface to the given process.
	/**	if the interface doesn't exist yet, it will be created.*/
		inline Interface& interface(int procID)			{return m_interfaceMap[procID];}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID)		{return m_interfaceMap.find(procID) != m_interfaceMap.end();}

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
		typedef typename LevelLayout::iterator		iterator;

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

	///	returns the iterator to the last interface of the layout in the given level.	
	/**	You should access the values of this iterator using the methods
	 *	Layout::interface and Layout::proc_id.
	 *	Make sure that the level matches the level in the associated begin() call.*/
		inline iterator end(size_t level)				{require_level(level); return m_vLayouts[level]->end();}

	///	returns the interface to the given iterator.
		inline Interface& interface(iterator& iter)		{return iter->second;}
		
	///	returns the interface to the given iterator.
		inline int proc_id(iterator& iter)				{return iter->first;}
	
	///	returns the number of levels.
		inline int num_levels()							{return m_vLayouts.size();}

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
		inline Interface& interface(int procID, size_t level)	{require_level(level); return m_vLayouts[level].interface(procID);}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID, size_t level)	{require_level(level); return m_vLayouts[level].interface_exists(procID);}
	
	///	returns the layout at the given level.
	/**	If level >= num_levels() then the layouts in between
		will be automatically created.*/
		inline LevelLayout& layout_on_level(int level)			{require_level(level); return *m_vLayouts[level];}

	protected:
	///	adds num new levels.
		inline void new_levels(int num)			{for(int i = 0; i < num; ++i) m_vLayouts.push_back(new LevelLayout);}
		
	///	if the required level doesn't exist yet, it will created.
		inline void require_level(int level)	{if(level >= num_levels()) new_levels(level - num_levels() + 1);}
		
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
//	LayoutMap
///	lets you access layouts by type and key
/**
 * The LayoutMap helps you to organize your layouts
 * (e.g. master- and slave-layouts).
 *
 * You may query layouts for any type. That makes this class very
 * flexible, however it requires consistent use of types throughout the
 * whole program.
 *
 * You may use a LayoutMap as follows (VertexBase is an arbitrary type):
 *
 * \code
 * LayoutMap<SinlgeLevelLayout, BasicInterface, int> layoutMap;
 * assert(!layoutMap.has_layout<VertexBase>(0));
 * SingleLevelLayout<BasicInterface<VertexBase> > > l = layoutMap.get_layout<VertexBase>(0);
 * assert(layoutMap.has_layout<VertexBase>(0));
 * \endcode
 *
 * It may be good to use some typedefs in order to make the code a little
 * easier to read:
 * \code
 * typedef LayoutMap<SingleLevelLayout, BasicInterface, int> SomeLayoutMap;
 * SomeLayoutMap::Types<VertexBase>::Layout l = layoutMap.get_layout<VertexBase>(0);
 * \endcode
 *
 * The Types struct is very useful when it comes to using a LayoutMap in
 * template code, too.
 */
 /*
 <class TType,
		  template<class T, class Alloc = std::allocator<T> >
			class TContainer = std::vector>
*/
template <template <class TInterface> class TLayout,
			template <class TType,
				template <class T, class Alloc> class TContainer = std::vector,
				template <class T> class TAlloc = std::allocator>
				class TInterface,
			class TKey,
			template<class T, class Alloc>
				class TInterfaceElemContainer = std::vector,
			template<class T> 
				class TAlloc = std::allocator>
			
class LayoutMap
{
	public:
		typedef TKey	Key;
		
	///	defines the types that are used by a LayoutMap for a given TType.
		template <class TType>
		struct Types
		{
			typedef TInterface<TType,
								TInterfaceElemContainer,
								TAlloc>	Interface;
			typedef TLayout<Interface>				Layout;
			typedef typename Interface::Element		Element;
			typedef std::map<TKey, Layout>			Map;
		};
		
	public:
	///	checks whether the layout associated with the given key exists for the given type.
		template <class TType>
		bool has_layout(TKey key)
		{
			typename Types<TType>::Map& m = get_layout_map<TType>();
			return m.find(key) != m.end();											
		}
		
	///	creates the required layout if it doesn't exist already.
		template <class TType>
		typename Types<TType>::Layout& get_layout(TKey key)
		{
			typename Types<TType>::Map& m = get_layout_map<TType>();
			return m[key];
		}
		
	protected:
		template <class TType>
		inline typename Types<TType>::Map&
		get_layout_map()
		{
			static typename Types<TType>::Map layoutMap;
			return layoutMap;
		}
};

}//	end of namespace pcl

#endif
