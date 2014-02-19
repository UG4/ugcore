// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// March 2013

#include "pcl/pcl_base.h"
#include "common/serialization.h"

#ifndef __H__PLG__COMPOL_BINARY_OR_ATTACHMENT__
#define __H__PLG__COMPOL_BINARY_OR_ATTACHMENT__

namespace ug
{

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
///	Combines received and existing values using a binary or operation
/**
 * TLayout::Type has to be either Vertex, Edge, Face or Volume.
 */
template <class TLayout, class TAttachment>
class ComPol_BinaryOrAttachment : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename TAttachment::ValueType Value;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_source before passing it to a communicator.*/
		ComPol_BinaryOrAttachment();
		
	///	Calls set_target internally.
		ComPol_BinaryOrAttachment(Grid& grid, TAttachment& attachment);

		virtual ~ComPol_BinaryOrAttachment()	{}
		
	///	The grid and the attachment from where the data shall be copied.
	/**	Make sure that attachment is attached to the correct
	 *	element-type of the grid.*/
		void set_attachment(Grid& grid, TAttachment& attachment);
		
	///	writes the data for the given interface to the buffer.
	/**	Derived from ICollector
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface);
		
	///	reads the data from the buffer to the given interface .
	/**	Derived from IExtractor
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface);
		
	protected:
		Grid::AttachmentAccessor<GeomObj, TAttachment>	m_aaVal;
};


/// @}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
ComPol_BinaryOrAttachment<TNodeLayout, TAttachment>::
ComPol_BinaryOrAttachment()
{
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
ComPol_BinaryOrAttachment<TNodeLayout, TAttachment>::
ComPol_BinaryOrAttachment(Grid& grid, TAttachment& attachment)
{
	set_attachment(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
void ComPol_BinaryOrAttachment<TNodeLayout, TAttachment>::
set_attachment(Grid& grid, TAttachment& attachment)
{
	m_aaVal.access(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
bool ComPol_BinaryOrAttachment<TNodeLayout, TAttachment>::
collect(ug::BinaryBuffer& buff, const Interface& interface)
{
	for(typename Interface::const_iterator iter = interface.begin();
		iter != interface.end(); ++iter)
	{
		Serialize(buff, m_aaVal[interface.get_element(iter)]);
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
bool ComPol_BinaryOrAttachment<TNodeLayout, TAttachment>::
extract(ug::BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(typename Interface::const_iterator iter = interface.begin();
		iter != interface.end(); ++iter)
	{
		Deserialize(buff, v);
		m_aaVal[interface.get_element(iter)] |= v;
	}
	return true;
}

};

#endif
