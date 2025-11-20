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

#ifndef __H__UG__compol_interface_status__
#define __H__UG__compol_interface_status__

#include <map>
#include <vector>
#include "common/assert.h"
#include "../distributed_grid.h"
#include "pcl/pcl_communication_structs.h"

namespace ug
{

///	Exchanges information on the interface-status of connected elements.
/**	The policy fills bool-vectors for each interface, telling whether
 * the associated elements on the connected processes contain the given status.
 * (at least one of both has to contain it).
 *
 * The vectors are organized in a map. The rank of the process, to which the
 * interface connects, is used as the key.
 */
template <typename TLayout>
class ComPol_InterfaceStatus : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		using Layout = TLayout;
		using GeomObj = typename Layout::Type;
		using Element = typename Layout::Element;
		using Interface = typename Layout::Interface;
		using InterfaceIter = typename Interface::const_iterator;

	public:
	/**
	 * \param status	(optional - default: INT_NONE). An or-combination of
	 * 					the constants enumerated in InterfaceNodeTypes.
	 * \{*/
		ComPol_InterfaceStatus(uint status = INT_NONE) :
			m_pLayout(nullptr),
			m_distGridMgr(nullptr),
			m_status(status),
			m_curLevel(0)	{}

		ComPol_InterfaceStatus(DistributedGridManager* distGridMgr,
								uint status = INT_NONE) :
			m_pLayout(nullptr),
			m_distGridMgr(distGridMgr),
			m_status(status),
			m_curLevel(0)	{}

		~ComPol_InterfaceStatus() override = default;

		void set_status(uint status)	{m_status = status;}
	/**	\} */

	///	returns the current status
		uint get_status()				{return m_status;}

	///	set the distributed grid manager
		void set_distributed_grid_manager(DistributedGridManager* distGridMgr)
		{m_distGridMgr = distGridMgr;}


	///	returns the vector which holds the results for the interface to the given process
	/**	This vector should not be accessed until after communication.
	 * After communication the vector has the size of the specified interface in the
	 * layout, for which communication was performed.*/
		const std::vector<bool>& get_result_vec(int procRank, int level = 0)
		{
			UG_ASSERT(m_pLayout, "communication has to be performed"
								" before results are checked.");

			UG_ASSERT(m_pLayout->interface_exists(procRank, level),
					"interface to proc " << procRank << " on level "
					<< level << " does not exist.");

			if(level >= m_vecMaps.size())
				m_vecMaps.resize(level + 1);

			std::vector<bool>& result = m_vecMaps[level][procRank];

			UG_ASSERT(result.size() == m_pLayout->interface(procRank, level).size(),
						"result size and interface size do not match."
						" perform communication directly before querying the results.");

			return result;
		}

		const std::map<int, std::vector<bool> >& get_result_map(int level = 0)
		{
			UG_ASSERT(m_pLayout, "communication has to be performed"
								" before results are checked.");

			if(level >= (int)m_vecMaps.size())
				m_vecMaps.resize(level + 1);

			return m_vecMaps[level];
		}

	////////////////////////
	//	COMMUNICATION STUFF
		int
		get_required_buffer_size(const Interface& interface) override {
			return interface.size() * sizeof(char);
		}

		bool
		begin_layout_extraction(const Layout* pLayout) override {
		//	clear and prepare the maps and set the layout
			m_pLayout = pLayout;

			m_vecMaps.clear();
			m_vecMaps.resize(m_pLayout->num_levels());

			return true;
		}

		bool
		collect(BinaryBuffer& buff, const Interface& interface) override {
		//	iterate over all elements in the interface and write
		//	whether they contain the status or not.
			UG_ASSERT(m_distGridMgr, "Please set the distributed grid manager"
									" before performing communication.");

			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				char val = (char)m_distGridMgr->contains_status(elem, m_status);
				buff.write(&val, sizeof(char));
			}

			return true;
		}

		void begin_level_extraction(int level) override {
			UG_ASSERT(level < (int)m_vecMaps.size(), "Internal pcl error.");
			m_curLevel = level;
		}

		bool
		extract(BinaryBuffer& buff, const Interface& interface) override {
		//	read the info from the buff and push bools to the
		//	vector associated with the interface.

			std::vector<bool>&
				vec = m_vecMaps[m_curLevel][interface.get_target_proc()];
			vec.clear();
			vec.reserve(interface.size());

			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				char val;
				buff.read(&val, sizeof(char));

				vec.push_back(((bool)val) ||
							  m_distGridMgr->contains_status(elem, m_status));
			}

			return true;
		}

	private:
		using VecMap = std::map<int, std::vector<bool> >;

		std::vector<VecMap>		m_vecMaps;
		const TLayout*			m_pLayout;
		DistributedGridManager*	m_distGridMgr;
		uint					m_status;
		int						m_curLevel;
};

}//	end of namespace

#endif
