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
 * TLayout::Type has to be either Vertex, Edge, Face or Volume.
 * TAttachment has to be of the type std::vector<SomeType> or compatible.
 */
template <class TLayout, class TAttachment>
class ComPol_GatherVecAttachment : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		using Layout = TLayout;
		using GeomObj = typename Layout::Type;
		using Element = typename Layout::Element;
		using Interface = typename Layout::Interface;
		using Vector = typename TAttachment::ValueType;
		using Value = typename Vector::value_type;
		
	///	Initialises the collector with an invalid grid.
	/**	be sure to call set_source before passing it to a communicator.*/
		ComPol_GatherVecAttachment();
		
	///	Calls set_target internally.
		ComPol_GatherVecAttachment(Grid& grid, TAttachment& attachment);

		~ComPol_GatherVecAttachment() override = default;
		
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
collect(BinaryBuffer& buff, const Interface& interface)
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
extract(BinaryBuffer& buff, const Interface& interface)
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
