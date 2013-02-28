// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m12 d09

#include "pcl/pcl_base.h"

#ifndef __H__PLG__COMPOL_COPY_ATTACHMENT__
#define __H__PLG__COMPOL_COPY_ATTACHMENT__

namespace ug
{

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
///	copies values from a specified attachment to a stream and back.
/**
 * TLayout::Type has to be either VertexBase, EdgeBase, Face or Volume.
 */
template <class TLayout, class TAttachment>
class ComPol_CopyAttachment : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename TAttachment::ValueType Value;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_source before passing it to a communicator.*/
		ComPol_CopyAttachment();
		
	///	Calls set_target internally.
		ComPol_CopyAttachment(Grid& grid, TAttachment& attachment);

		virtual ~ComPol_CopyAttachment()	{}
		
	///	The grid and the attachment from where the data shall be copied.
	/**	Make sure that attachment is attached to the correct
	 *	element-type of the grid.*/
		void set_attachment(Grid& grid, TAttachment& attachment);
		
	///	writes the data for the given interface to the buffer.
	/**	Derived from ICollector
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		virtual bool
		collect(ug::BinaryBuffer& buff, Interface& interface);
		
	///	reads the data from the buffer to the given interface .
	/**	Derived from IExtractor
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		virtual bool
		extract(ug::BinaryBuffer& buff, Interface& interface);
		
	protected:
		Grid::AttachmentAccessor<GeomObj, TAttachment>	m_aaVal;
};


/// @}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of the methods of CollectorCopy
////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
ComPol_CopyAttachment<TNodeLayout, TAttachment>::
ComPol_CopyAttachment()
{
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
ComPol_CopyAttachment<TNodeLayout, TAttachment>::
ComPol_CopyAttachment(Grid& grid, TAttachment& attachment)
{
	set_attachment(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
void ComPol_CopyAttachment<TNodeLayout, TAttachment>::
set_attachment(Grid& grid, TAttachment& attachment)
{
	m_aaVal.access(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
bool ComPol_CopyAttachment<TNodeLayout, TAttachment>::
collect(ug::BinaryBuffer& buff, Interface& interface)
{
	for(typename Interface::iterator iter = interface.begin();
		iter != interface.end(); ++iter)
		buff.write((char*)&m_aaVal[interface.get_element(iter)], sizeof(Value));
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
bool ComPol_CopyAttachment<TNodeLayout, TAttachment>::
extract(ug::BinaryBuffer& buff, Interface& interface)
{
	for(typename Interface::iterator iter = interface.begin();
		iter != interface.end(); ++iter)
		buff.read((char*)&m_aaVal[interface.get_element(iter)], sizeof(Value));
	return true;
}

};

#endif
