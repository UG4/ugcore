/**
 * \file rsamg_debug.h
 * \author Martin Rupp
 *
 *
 * \date 16.06.10
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_NODEINFO_H__
#define __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_NODEINFO_H__

//#include "maxheap.h"
#include "../boxsort.h"
#include "../amg_debug_helper.h"

namespace ug {

//  structs
#define AMG_UNKNOWN							(-1)
#define AMG_ASSIGNED_RATING					(-2)
#define AMG_COARSE_RATING					(-3)
#define AMG_UNASSIGNED_FINE_INDIRECT_RATING	(-4)
#define AMG_FINE_RATING						(-6)
class AMGNode
{
private:

	//int newIndex;		
	inline void set_assigned(){rating = AMG_ASSIGNED_RATING;}		
	inline void set_coarse()	{rating = AMG_COARSE_RATING;}
	inline void set_fine_direct(){rating = AMG_FINE_RATING;	}
	
	inline void set_unassigned_fine_indirect(){rating = AMG_UNASSIGNED_FINE_INDIRECT_RATING;}
	inline void set_fine_indirect_level(int level) { rating = AMG_FINE_RATING+1-level;}
	
public:
	AMGNode() { rating = AMG_UNKNOWN; }
	int rating;
	inline bool is_coarse() const {	return rating == AMG_COARSE_RATING;}
	inline bool is_fine_direct() const {return (rating == AMG_FINE_RATING);}
	inline bool operator ==(int comp)
	{
		return comp == rating;
	}
	
	inline bool is_unassigned_fine_indirect() const {return rating == AMG_UNASSIGNED_FINE_INDIRECT_RATING;}
	
	
	inline bool is_fine_indirect_level(int level) const { return rating == AMG_FINE_RATING+1-level;}
	

	inline bool is_assigned() const {return (rating <= AMG_ASSIGNED_RATING);}

	friend std::ostream &operator << (std::ostream &out, const AMGNode &n)
	{
		out << "Rating: " << n.rating;
		if(n.rating < 0)
		{
			if(n.is_coarse()) out << " (coarse)";
			else if(n.is_unassigned_fine_indirect()) out << "(unknown indirect fine)";
			else if(n.is_fine_direct())
				out << " (fine direct)";
			else
				out << " (indirect level " << AMG_FINE_RATING+1-n.rating << ")";
		}
		out << " ";
		return out;
	}
	void print() const
	{
		std::cout << *this << std::endl;
	} // << " newindex: " << newIndex << endl;
	
	inline bool operator > (const AMGNode &other) const
	{
		if(rating == other.rating)
			return this < &other; // we somehow want a STABLE sort, for that coarsening is in the direction of the numbering of the elements
		else
			return rating > other.rating;
	}

	inline size_t get_val() const
	{
		UG_ASSERT(rating >= 0 && rating < 1000, "rating is " << rating << ", out of bounds [0, 1000]");
		return (int)rating;
	}

	friend class AMGNodes;
};


class AMGNodes
{
public:
	typedef AMGNode value_type;
	AMGNodes()
	{
		m_unassigned = 0;
		m_iNrOfCoarse = 0;
	}

	AMGNodes(size_t N)
	{
		m_unassigned = 0;
		m_iNrOfCoarse = 0;
		resize(N);
	}


	void resize(size_t N)
	{
		//size_t oldsize = size();
		m_nodes.resize(N);

	}

	AMGNode &operator[] (size_t i)
	{
		return m_nodes[i];
	}

	const AMGNode &operator[] (size_t i) const
	{
		return m_nodes[i];
	}


	size_t size() const
	{
		return m_nodes.size();
	}

	inline void set_assigned(size_t i)
	{
		m_nodes[i].set_assigned();
		UG_ASSERT(m_unassigned != 0, i);
		m_unassigned--;
	}

	inline void external_set_coarse(size_t i)
	{
		//UG_ASSERT(!is_slave(i), i);
		if(m_nodes[i].is_assigned() == false)
		{
			UG_ASSERT(m_unassigned != 0, i);
			m_unassigned--;
		}
		else if(m_nodes[i].is_coarse()) return;
		m_nodes[i].set_coarse();
		m_iNrOfCoarse++;
	}

	inline void set_coarse(size_t i)
	{
		external_set_coarse(i);
	}
	inline void set_fine_direct(size_t i)
	{
		m_nodes[i].set_fine_direct();
		UG_ASSERT(m_unassigned != 0, i);
		m_unassigned--;
	}

	inline void set_unassigned_fine_indirect(size_t i)
	{
		UG_ASSERT(m_nodes[i].is_fine_direct(), "");
		m_nodes[i].set_unassigned_fine_indirect();
		m_unassigned++;
	}

	inline void set_fine_indirect_level(size_t i, int level)
	{
		UG_ASSERT(m_nodes[i].is_unassigned_fine_indirect(), "");
		m_nodes[i].set_fine_indirect_level(level);
		m_unassigned--;
	}

	inline void set_isolated(size_t i)
	{
		if(m_nodes[i].rating >=0)
			m_unassigned--;
		m_nodes[i].set_fine_direct();
	}


	size_t get_unassigned() const
	{
		return m_unassigned;
	}

	size_t get_nr_of_coarse() const
	{
		return m_iNrOfCoarse;
	}

	void set_rating(size_t i, int rating)
	{
		UG_ASSERT(rating >=0, rating << " i = " << i);
		//UG_ASSERT(!is_slave(i), i);
		m_unassigned++;
		if(m_nodes[i].is_coarse()) m_iNrOfCoarse--;
		m_nodes[i].rating = rating;
	}

	void print_ratings(const cAMG_helper &amghelper, size_t level) const
	{
		UG_LOG("Coarsen ratings:\n")
		for(size_t i=0; i<size(); i++)
			UG_LOG(i << " [" << amghelper.GetOriginalIndex(level, i) << "] " << m_nodes[i] << std::endl);
	}




#ifdef UG_PARALLEL
public:
	bool is_inner(size_t i) const
	{
		return !bMaster[i] && !bSlave[i];
	}
	bool is_master(size_t i) const
	{
		return bMaster[i];
	}

	bool is_slave(size_t i) const
	{
		return bSlave[i];
	}

	stdvector<bool> bMaster;
	stdvector<bool> bSlave;
#else
public:
	bool is_inner(size_t i) const
	{
		return true;
	}
	bool is_master(size_t i) const
	{
		return true;
	}

	bool is_slave(size_t i) const
	{
		return false;
	}
#endif

private:
	stdvector<AMGNode> m_nodes;
	int m_unassigned;
	int m_iNrOfCoarse;			//< nr of nodes to assign
};

//	typedef maxheap<AMGNode> nodeinfo_pq_type;
typedef BoxPriorityQueue<AMGNode> nodeinfo_pq_type;

}


#endif // __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_NODEINFO_H__
