// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 19.05.2011 (m,d,y)

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
template <class TLayout>
class ComPol_InterfaceStatus : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

	public:
	/**
	 * \param status	(optional - default: INT_NONE). An or-combination of
	 * 					the constants enumerated in InterfaceNodeTypes.
	 * \{*/
		ComPol_InterfaceStatus(uint status = INT_NONE) :
			m_pLayout(NULL),
			m_distGridMgr(NULL),
			m_status(status),
			m_curLevel(0)	{}

		ComPol_InterfaceStatus(DistributedGridManager* distGridMgr,
								uint status = INT_NONE) :
			m_pLayout(NULL),
			m_distGridMgr(distGridMgr),
			m_status(status),
			m_curLevel(0)	{}

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
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			return interface.size() * sizeof(char);
		}

		virtual bool
		begin_layout_extraction(const Layout* pLayout)
		{
		//	clear and prepare the maps and set the layout
			m_pLayout = pLayout;

			m_vecMaps.clear();
			m_vecMaps.resize(m_pLayout->num_levels());

			return true;
		}

		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
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

		virtual void begin_level_extraction(int level)
		{
			UG_ASSERT(level < (int)m_vecMaps.size(), "Internal pcl error.");
			m_curLevel = level;
		}

		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
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
		typedef std::map<int, std::vector<bool> > VecMap;

		std::vector<VecMap>		m_vecMaps;
		const TLayout*			m_pLayout;
		DistributedGridManager*	m_distGridMgr;
		uint					m_status;
		int						m_curLevel;
};

}//	end of namespace

#endif
