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

#ifndef __H__UG__compol_subset__
#define __H__UG__compol_subset__
#include "pcl/pcl_communication_structs.h"

namespace ug
{


template <typename TLayout>
class ComPol_Subset : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		using Layout = TLayout;
		using GeomObj = typename Layout::Type;
		using Element = typename Layout::Element;
		using Interface = typename Layout::Interface;
		using InterfaceIter = typename Interface::const_iterator;

	///	Construct the communication policy with a ug::SubsetHandler.
	/**	If overwrite is specified, the subset index in the target handler is
	 * simply overwritten (default is false).*/
		ComPol_Subset(ISubsetHandler& sel, bool overwrite = false)
			 :	m_sh(sel), m_overwriteEnabled(overwrite)
		{}

		~ComPol_Subset() override = default;

		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			return interface.size() * sizeof(int);
		}

	///	writes 1 for selected and 0 for unassigned interface entries
		virtual bool
		collect(BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				int si = m_sh.get_subset_index(elem);
				buff.write((char*)&si, sizeof(int));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(BinaryBuffer& buff, const Interface& interface)
		{
			int nsi = -1;
			bool retVal = true;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&nsi, sizeof(int));
				if(m_overwriteEnabled){
					m_sh.assign_subset(elem, nsi);
				}
				else{
					if(m_sh.get_subset_index(elem) == -1){
						m_sh.assign_subset(elem, nsi);
					}
					else if(m_sh.get_subset_index(elem) != nsi){
					//	if the subset indices do not match, we have a problem here.
						retVal = false;
					}
				}
			}
			return retVal;
		}

	protected:
		ISubsetHandler&	m_sh;
		bool m_overwriteEnabled;
};

}//	end of namespace

#endif
