//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

#ifndef __H__PCL__PCL_BASE__
#define __H__PCL__PCL_BASE__

#include <vector>
#include <map>

namespace pcl
{

////////////////////////////////////////////////////////////////////////
//	group_traits
///	collection of types for a ElementGroup
/**
 * The group_traits are used throughout pcl in order to define
 * the correct parameter- and return-value-types.
 *
 * The following types have to be defined by a specialization of
 * the group_traits:
 * - Element:		THe type of the elements which the ElementGroup contains.
 * - LocalID:		The type of a local id.
 * - Interface:		The type of an Interface.
 * - Layout:		The type of a Layout
 *
 * if your element-group already defines thos types, you may use the
 * default implementation of the group_traits. If not you can create
 * a specialization for your group through template specialization.
 */
template <class TElementGroup>
class group_traits
{
	public:
		typedef typename TElementGroup::Element		Element;
		typedef typename TElementGroup::LocalID		LocalID;
		typedef typename TElementGroup::Interface	Interface;
		typedef typename TElementGroup::Layout		Layout;
};



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
	///	the interface type - just for conveniance
		typedef TInterface	Interface;

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
