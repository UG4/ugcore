/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__compol_selection__
#define __H__UG__compol_selection__
#include "pcl/pcl_communication_structs.h"
#include "common/types.h"

namespace ug
{


template <class TLayout>
class ComPol_Selection : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	///	Construct the communication policy with a ug::Selector.
	/**	Through the parameters select and deselect one may specify whether
	 * a process selects and/or deselects elements based on the received
	 * selection status.*/
		ComPol_Selection(ISelector& sel, bool select = true, bool deselect = true)
			 :	m_sel(sel), m_bSelectAllowed(select), m_bDeselectAllowed(deselect)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(byte_t);}

	///	writes 1 for selected and 0 for unselected interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				byte_t refMark = m_sel.get_selection_status(elem);
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
				if(val && select_allowed())
					m_sel.select(elem, val);
				else if(!val && deselect_allowed())
					m_sel.deselect(elem);
			}
			return true;
		}

	protected:
		inline bool select_allowed()	{return m_bSelectAllowed;}
		inline bool deselect_allowed()	{return m_bDeselectAllowed;}

		ISelector& m_sel;
		bool m_bSelectAllowed;
		bool m_bDeselectAllowed;
};



////////////////////////////////////////////////////////////////////////////////
/**	Enables only the part of the selection-state for which the specified stateBits
 * hold 1.
 */
template <class TLayout>
class ComPol_EnableSelectionStateBits : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		ComPol_EnableSelectionStateBits(ISelector& sel, byte_t stateBits)
			 :	m_sel(sel), m_stateBits(stateBits)
		{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(byte_t);}

	///	writes writes the selection states of the interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				byte_t refMark = m_sel.get_selection_status(elem);
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
				m_sel.select(elem, m_sel.get_selection_status(elem)
								   | (val & m_stateBits));
			}
			return true;
		}

		ISelector& m_sel;
		byte_t m_stateBits;
};

}//	end of namespace

#endif
