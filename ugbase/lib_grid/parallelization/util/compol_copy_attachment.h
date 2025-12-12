/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__PLG__COMPOL_COPY_ATTACHMENT__
#define __H__PLG__COMPOL_COPY_ATTACHMENT__

#include "lib_grid/grid/grid.h"  // for Grid::AttachmentAccessor
// #include "pcl/pcl_base.h"
#include "pcl/pcl_communication_structs.h"  // for ICommunicationPolicy
#include "common/serialization.h"


namespace ug {

/// \addtogroup lib_grid_parallelization
/// @{

////////////////////////////////////////////////////////////////////////
///	copies values from a specified attachment to a stream and back.
/**
 * TLayout::Type has to be either Vertex, Edge, Face or Volume.
 */
template <typename TLayout, typename TAttachment>
class ComPol_CopyAttachment : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		using Layout = TLayout;
		using GeomObj = typename Layout::Type;
		using Element = typename Layout::Element;
		using Interface = typename Layout::Interface;
		using Value = typename TAttachment::ValueType;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_source before passing it to a communicator.*/
		ComPol_CopyAttachment();
		
	///	Calls set_target internally.
		ComPol_CopyAttachment(Grid& grid, TAttachment attachment);

		~ComPol_CopyAttachment() override = default;
		
	///	The grid and the attachment from where the data shall be copied.
	/**	Make sure that attachment is attached to the correct
	 *	element-type of the grid.*/
		void set_attachment(Grid& grid, TAttachment& attachment);
		
	///	writes the data for the given interface to the buffer.
	/**	Derived from ICollector
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		bool
		collect(BinaryBuffer& buff, const Interface& interface) override;
		
	///	reads the data from the buffer to the given interface .
	/**	Derived from IExtractor
	 *	Make sure that all members of the interface are members of the
	 *	grid too.*/
		bool
		extract(BinaryBuffer& buff, const Interface& interface) override;
		
		void extract_on_constrained_elems_only(bool enable);

	protected:
		Grid::AttachmentAccessor<GeomObj, TAttachment>	m_aaVal;
		bool m_extractOnConstrainedElemsOnly;
};


/// @}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of the methods of CollectorCopy
////////////////////////////////////////////////////////////////////////
template <typename TNodeLayout, typename TAttachment>
ComPol_CopyAttachment<TNodeLayout, TAttachment>::
ComPol_CopyAttachment() :
	m_extractOnConstrainedElemsOnly(false)
{
}

////////////////////////////////////////////////////////////////////////
template <typename TNodeLayout, typename TAttachment>
ComPol_CopyAttachment<TNodeLayout, TAttachment>::
ComPol_CopyAttachment(Grid& grid, TAttachment attachment) :
	m_extractOnConstrainedElemsOnly(false)
{
	set_attachment(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <typename TNodeLayout, typename TAttachment>
void ComPol_CopyAttachment<TNodeLayout, TAttachment>::
set_attachment(Grid& grid, TAttachment& attachment)
{
	m_aaVal.access(grid, attachment);
}

////////////////////////////////////////////////////////////////////////
template <typename TNodeLayout, typename TAttachment>
bool ComPol_CopyAttachment<TNodeLayout, TAttachment>::
collect(BinaryBuffer& buff, const Interface& interface)
{
	for(typename Interface::const_iterator iter = interface.begin();
		iter != interface.end(); ++iter)
	{
		Serialize(buff, m_aaVal[interface.get_element(iter)]);
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
template <typename TNodeLayout, typename TAttachment>
bool ComPol_CopyAttachment<TNodeLayout, TAttachment>::
extract(BinaryBuffer& buff, const Interface& interface)
{
	if(m_extractOnConstrainedElemsOnly){
		for(typename Interface::const_iterator iter = interface.begin();
			iter != interface.end(); ++iter)
		{
			typename TAttachment::ValueType val;
			Deserialize(buff, val);
			if(interface.get_element(iter)->is_constrained())
				m_aaVal[interface.get_element(iter)] = val;
		}
	}
	else{
		for(typename Interface::const_iterator iter = interface.begin();
			iter != interface.end(); ++iter)
		{
			Deserialize(buff, m_aaVal[interface.get_element(iter)]);
		}
	}
	return true;
}

template <typename TNodeLayout, typename TAttachment>
void ComPol_CopyAttachment<TNodeLayout, TAttachment>::
extract_on_constrained_elems_only(bool enable)
{
	m_extractOnConstrainedElemsOnly = enable;
}

};

#endif
