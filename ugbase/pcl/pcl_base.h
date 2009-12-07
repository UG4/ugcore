//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d05

/**/////////////////////////////////////////////////////////////////////
//	Base of the 'parallel communication layer'.
//	...
*///////////////////////////////////////////////////////////////////////

#ifndef __H__PCL__PCL_BASE__
#define __H__PCL__PCL_BASE__

#include <vector>
#include <map>

namespace pcl
{

////////////////////////////////////////////////////////////////////////
//	pcl_traits
///	collection of types for a ElementGroup
/**
 * The pcl_traits are used throughout pcl in order to define
 * the correct parameter- and return-value-types.
 *
 * The following types have to be defined by a specialization of
 * the pcl_traits:
 * - Element:		The elements of which the ElementGroup consists.
 * - ElementRef:	A reference-type to an Element.
 * - LocalID:		The type of a local id.
 * - Interface:		The type of an Interface.
 * - Layout:		The type of a Layout
 */
template <class TElementGroup>
class group_traits;



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
		typedef InterfaceMap::iterator		iterator;

	public:
	//	methods
	///	returns the interface to the given process.
		inline Interface& interface(int procID)			{return m_interfaceMap[procID];}

	///	returns the interface to the given iterator.
		inline Interface& interface(iterator& iter)		{return iter->second;}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID)		{return return m_interfaceMap.find(procID) != m_streamMap.end();}

	///	returns the interface to the given iterator.
		inline int& proc_id(iterator& iter)				{return iter->first;}

	///	returns true if an interface to the given procID already exists.
		inline bool interface_exists(int procID)		{return return m_interfaceMap.find(procID) != m_streamMap.end();}
		
	///	returns the iterator to the first interface of the layout.
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator begin()							{return m_interfaceMap.begin();}

	///	returns the iterator to the last interface of the layout.	
	/**	You should access the values of this iterator using the methods
		Layout::interface and Layout::proc_id.*/
		inline iterator end()							{return m_interfaceMap.begin();}

	protected:
	///	holds the interfaces in a map.
		InterfaceMap	m_interfaceMap;
};

}//	end of namespace pcl

#endif
