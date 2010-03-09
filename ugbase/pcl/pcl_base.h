//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_BASE__
#define __H__PCL__PCL_BASE__

#include <vector>
#include <list>
#include <map>

namespace pcl
{
////////////////////////////////////////////////////////////////////////
//	Interface Categories
/**	Every interface has to feature this tag.
 *	This tag simply says that one may iterate over the elements in the
 *	interface.
 *
 *	Methods that have to be featured in the interface:
 *		- iterator begin();
 *		- iterator end();
 *		- const_iterator begin() const;
 *		- const_iterator end() const;
 *		- TElem& get_element(iterator iter);
 *		- const TElem& get_element(iterator iter) const;
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
 *	Methods that have to be featured for the interface:
 *		- static bool cmp(iterator iter1, iterator iter2);
 */
class ordered_interface_tag : public basic_interface_tag	{};



////////////////////////////////////////////////////////////////////////
//	BasicInterface
///	You may add elements to this interface and iterate over them
template <class TElem,
		  template<class T, class Alloc = std::allocator<T> >
			class TContainer = std::vector>
class BasicInterface
{
	protected:
		typedef TContainer<TElem>				ElemContainer;
		
	public:
		typedef basic_interface_tag				category_tag;
		
		typedef TElem									Element;
		typedef typename ElemContainer::iterator		iterator;
		typedef typename ElemContainer::const_iterator	const_iterator;
		
	public:
		BasicInterface()	: m_size(0)		{};
		inline iterator push_back(const TElem& elem)	{++m_size; return m_elements.insert(m_elements.end(), elem);}
		inline void erase(iterator iter)				{--m_size; m_elements.erase(iter);}
		
		inline iterator begin()		{return m_elements.begin();}
		inline iterator end()		{return m_elements.end();}
		
		inline TElem& get_element(iterator iter)	{return *iter;}
		
		inline size_t size()	{return m_size;}
		
	protected:
		ElemContainer	m_elements;
		size_t			m_size;
};


////////////////////////////////////////////////////////////////////////
//	OrderedInterface
///	You may add elements to this interface and iterate over them
/**	You may compare iterators of elements in this interface using the
 *	cmp method.
 *
 *	This interface stores elements in a vector internally.
 *	This allows for easy compare but makes erase expensive.*/
/*
template <class TElem, template<class> class TContainer = std::vector>
class OrderedInterface
{
	protected:
	//	this struct holds an elem and a local id, as they are
	//	required if we use non-random-access-iterators
		template <typename TIterTag> struct TInterfaceEntryByTag
		{
			TElem	elem;
			uint	localID;
		};
		
	//	specialization for random-access-iterators. A local id is not required here.
		template <> struct
		TInterfaceEntryByTag<std::random_access_iterator_tag>
		{
			TElem elem;
		};
		
	//	this is required to pick the right InterfaceEntryType
		typedef typename
			iterator_traits<typename TContainer<TElem>::iterator>::iterator_category
			iterator_category;
			
		typedef TInterfaceEntryByTag<iterator_category>	InterfaceEntry;
		typedef TContainer<InterfaceEntry>				ElemContainer;
		
	public:
		typedef ordered_interface_tag			category_tag;
		
		typedef ElemContainer::iterator			iterator;
		typedef ElemContainer::const_iterator	const_iterator;
		
	public:
		iterator add_element(const TElem& elem);
		
		iterator begin();
		iterator end();
		
		TElem& get_element(iterator iter);
		
	///	returns true if iter1 < iter2.
		static bool cmp(iterator iter1, iterator iter2);
		
	protected:
		ElemContainer	m_elements;
};
*/

////////////////////////////////////////////////////////////////////////
//	Layout
///	the standard layout implementation
/**
 * A Layout is nothing more than a collection of interfaces.
 * Each interface is associated with a process-id.
 */
template <class TInterface>
class Layout
{
	public:
	//	typedefs
	///	the interface type
		typedef TInterface					Interface;
	///	Element type
		typedef typename Interface::Element	Element;
		
	///	an interface-map is a list of interfaces, each associated with a process id.
		typedef std::map<int, Interface>	InterfaceMap;

	///	An iterator that allows to iterate over the interfaces stored in the layout.
		typedef typename InterfaceMap::iterator		iterator;

	public:
	//	methods
	///	returns the interface to the given process.
		inline Interface& interface(int procID)			{return m_interfaceMap[procID];}

	///	returns the interface to the given iterator.
		inline Interface& interface(iterator& iter)		{return iter->second;}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID)		{return m_interfaceMap.find(procID) != m_interfaceMap.end();}

	///	returns the interface to the given iterator.
		inline int proc_id(iterator& iter)				{return iter->first;}

	///	returns the iterator to the first interface of the layout.
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator begin()							{return m_interfaceMap.begin();}

	///	returns the iterator to the last interface of the layout.	
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator end()							{return m_interfaceMap.end();}

	protected:
	///	holds the interfaces in a map.
		InterfaceMap	m_interfaceMap;
};

}//	end of namespace pcl

#endif
