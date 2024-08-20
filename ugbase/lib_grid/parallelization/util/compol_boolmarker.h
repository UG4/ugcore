/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__compol_boolmarker__
#define __H__UG__compol_boolmarker__

#include "pcl/pcl_communication_structs.h"
#include "lib_grid/tools/bool_marker.h"

namespace ug
{

///	adds marking at extracting side
template <class TLayout>
class ComPol_BoolMarker_AddMarks : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	///	Construct the communication policy with a ug::BoolMarker.
		ComPol_BoolMarker_AddMarks(BoolMarker& marker)
			 :	m_rMarker(marker)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)		{return interface.size() * sizeof(byte_t);}

	///	writes 1 for marked and 0 for unmarked interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				bool isMarked = m_rMarker.is_marked(elem);
				byte_t refMark = (isMarked ? 1 : 0);
				buff.write((char*)&refMark, sizeof(byte_t));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			byte_t val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte_t));
				if(val)	m_rMarker.mark(elem);
			}
			return true;
		}

	protected:
		BoolMarker& m_rMarker;
};

///	removes marks at extracting side, if no mark was received
template <class TLayout>
class ComPol_BoolMarker_RemoveMarks : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	///	Construct the communication policy with a ug::BoolMarker.
		ComPol_BoolMarker_RemoveMarks(BoolMarker& marker)
			 :	m_rMarker(marker)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)		{return interface.size() * sizeof(byte_t);}

	///	writes 1 for marked and 0 for unmarked interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				bool isMarked = m_rMarker.is_marked(elem);
				byte_t refMark = (isMarked ? 1 : 0);
				buff.write((char*)&refMark, sizeof(byte_t));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			byte_t val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte_t));
				if(!val)	m_rMarker.unmark(elem);
			}
			return true;
		}

	protected:
		BoolMarker& m_rMarker;
};

}//	end of namespace

#endif
