// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m12 d09

#include "lib_grid/lg_base.h"
#include "pcl/pcl_communication.h"

#ifndef __H__LIB_GRID__COLLECTOR_COPY__
#define __H__LIB_GRID__COLLECTOR_COPY__

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CollectorCopy
///	copies values from a specified attachment to a stream.
/**
 * A specialization of pcl::group_traits for TElementGroup has to exists.
 * pcl::group_traits::LocalID has to match either VertexBase*, EdgeBase*,
 * Face* or Volume*.
 */
template <class TElementGroup, class TAttachment>
class CollectorCopy : public pcl::ICollector<TElementGroup>
{
	public:
		typedef typename pcl::group_traits<TElementGroup>::Element		Element;
		typedef typename pcl::group_traits<TElementGroup>::LocalID		LocalID;
		typedef typename pcl::group_traits<TElementGroup>::Interface	Interface;
		typedef typename TAttachment::ValueType Value;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_source before passing it to a communicator.*/
		CollectorCopy();
		
	///	Calls set_target internally.
		CollectorCopy(Grid& grid, TAttachment& attachment);
		
	///	The grid and the attachment from where the data shall be copied.
	/**	Make sure that attachment is attached to the correct
	 *	element-type of the grid.*/
		void set_source(Grid& grid, TAttachment& attachment);
		
	///	writes the data for the given interface to the buffer.
	/**	Derived from ICollector
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		virtual bool
		collect(std::ostream& buff, Interface& interface);
		
	protected:
		Grid::AttachmentAccessor<Element, TAttachment>	m_aaVal;
};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of the methods of CollectorCopy
////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
CollectorCopy<TElementGroup, TAttachment>::
CollectorCopy()
{
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
CollectorCopy<TElementGroup, TAttachment>::
CollectorCopy(Grid& grid, TAttachment& attachment)
{
	set_source(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
void
CollectorCopy<TElementGroup, TAttachment>::
set_source(Grid& grid, TAttachment& attachment)
{
	m_aaVal.access(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
bool
CollectorCopy<TElementGroup, TAttachment>::
collect(std::ostream& buff, Interface& interface)
{
	for(size_t i = 0; i < interface.size(); ++i)
		buff.write((char*)&m_aaVal[interface[i]], sizeof(Value));
}

};

#endif
