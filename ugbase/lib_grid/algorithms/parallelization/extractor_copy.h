// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m12 d09

#include "lib_grid/lg_base.h"
#include "pcl/pcl_communication.h"

#ifndef __H__LIB_GRID__EXTRACTOR_COPY__
#define __H__LIB_GRID__EXTRACTOR_COPY__

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	ExtractorCopy
///	copies values from a stream to a specified attachment.
/**
 * This extractor is suited for MASTER->SLAVE communication, since each
 * value is simply copied.
 * If you're using it for SLAVE->MASTER communication, please note that
 * a master can receive values from multiple slaves. In this case received 
 * values are simply overwritten.
 *
 * A specialization of pcl::group_traits for TElementGroup has to exists.
 * pcl::group_traits::LocalID has to match either VertexBase*, EdgeBase*,
 * Face* or Volume*.
 */
template <class TElementGroup, class TAttachment>
class ExtractorCopy : public pcl::IExtractor<TElementGroup>
{
	public:
		typedef typename pcl::group_traits<TElementGroup>::Element		Element;
		typedef typename pcl::group_traits<TElementGroup>::LocalID		LocalID;
		typedef typename pcl::group_traits<TElementGroup>::Interface	Interface;
		typedef typename TAttachment::ValueType Value;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_target before passing it to a communicator.*/
		ExtractorCopy();
		
	///	Calls set_target internally.
		ExtractorCopy(Grid& grid, TAttachment& attachment);
		
	///	The grid and the attachment to where the data shall be copied.
	/**	If the attachment is not already attached to the correct
	 *	element type, it will be attached here.*/
		void set_target(Grid& grid, TAttachment& attachment);
		
	///	reads the data from the buffer to the given interface .
	/**	Derived from IExtractor
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		virtual bool
		extract(std::istream& buff, Interface& interface);
				
	protected:
		Grid::AttachmentAccessor<Element, TAttachment>	m_aaVal;
};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of the methods of ExtractorCopy
////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
ExtractorCopy<TElementGroup, TAttachment>::
ExtractorCopy()
{
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
ExtractorCopy<TElementGroup, TAttachment>::
ExtractorCopy(Grid& grid, TAttachment& attachment)
{
	set_target(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
void
ExtractorCopy<TElementGroup, TAttachment>::
set_target(Grid& grid, TAttachment& attachment)
{
	if(!grid.has_attachment<Element>(attachment))
		grid.attach_to<Element>(attachment);
	m_aaVal.access(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TElementGroup, class TAttachment>
bool
ExtractorCopy<TElementGroup, TAttachment>::
extract(std::istream& buff, Interface& interface)
{
	for(size_t i = 0; i < interface.size(); ++i)
		buff.read((char*)&m_aaVal[interface[i]], sizeof(Value));
}

};

#endif
