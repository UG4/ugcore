// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m12 d09

#include "pcl/pcl_base.h"

#ifndef __H__PLG__COMPOL_GATHER_VEC_ATTACHMENT__
#define __H__PLG__COMPOL_GATHER_VEC_ATTACHMENT__

namespace ug
{

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
///	Gathers the values stored in vector-attachments
/**
 * TLayout::Type has to be either VertexBase, EdgeBase, Face or Volume.
 * TAttachment has to be of the type std::vector<SomeType> or compatible.
 */
template <class TLayout, class TAttachment>
class ComPol_GatherVecAttachment : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename TAttachment::ValueType Vector;
		typedef typename Vector::value_type		Value;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_source before passing it to a communicator.*/
		ComPol_GatherVecAttachment();
		
	///	Calls set_target internally.
		ComPol_GatherVecAttachment(Grid& grid, TAttachment& attachment);

		virtual ~ComPol_GatherVecAttachment()	{}
		
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
		Grid::AttachmentAccessor<GeomObj, TAttachment>	m_aaVec;
};


/// @}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of the methods of CollectorCopy
////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
ComPol_GatherVecAttachment<TNodeLayout, TAttachment>::
ComPol_GatherVecAttachment()
{
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
ComPol_GatherVecAttachment<TNodeLayout, TAttachment>::
ComPol_GatherVecAttachment(Grid& grid, TAttachment& attachment)
{
	set_attachment(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
void ComPol_GatherVecAttachment<TNodeLayout, TAttachment>::
set_attachment(Grid& grid, TAttachment& attachment)
{
	m_aaVec.access(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
bool ComPol_GatherVecAttachment<TNodeLayout, TAttachment>::
collect(ug::BinaryBuffer& buff, const Interface& interface)
{
	for(typename Interface::const_iterator iter = interface.begin();
		iter != interface.end(); ++iter)
	{
		Serialize(buff, m_aaVec[interface.get_element(iter)]);
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TNodeLayout, class TAttachment>
bool ComPol_GatherVecAttachment<TNodeLayout, TAttachment>::
extract(ug::BinaryBuffer& buff, const Interface& interface)
{
	Vector tvec;
	for(typename Interface::const_iterator iter = interface.begin();
		iter != interface.end(); ++iter)
	{
		Element e = interface.get_element(iter);
		tvec.clear();
		Deserialize(buff, tvec);
		m_aaVec[e].reserve(m_aaVec[e].size() + tvec.size());
		for(size_t i = 0; i < tvec.size(); ++i)
			m_aaVec[e].push_back(tvec[i]);
	}
	return true;
}

};

#endif
