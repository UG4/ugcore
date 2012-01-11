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
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_nodes.h"
#endif

namespace ug {

//  structs
#define AMG_UNKNOWN							(-1)
#define AMG_ASSIGNED_RATING					(-2)
#define AMG_PARALLEL_DONT_CARE				(-3)
#define AMG_COARSE_RATING					(-4)
#define AMG_UNASSIGNED_FINE_INDIRECT_RATING	(-5)
#define AMG_DIRICHLET_RATING				(-6)
#define AMG_FINE_RATING						(-10)

class AMGNode
{
private:

	//int newIndex;		
	inline void set_assigned(){rating = AMG_ASSIGNED_RATING;}		
	inline void set_coarse()	{rating = AMG_COARSE_RATING;}
	inline void set_fine_direct(){rating = AMG_FINE_RATING;	}
	inline bool is_uncalculated_fine() { return false; }


	inline void set_unassigned_fine_indirect(){rating = AMG_UNASSIGNED_FINE_INDIRECT_RATING;}
	// indirect level 1 == fine direct
	inline void set_fine_indirect_level(int level) { rating = AMG_FINE_RATING+1-level;}
	inline void set_dirichlet() { rating = AMG_DIRICHLET_RATING; }

public:
	inline void set_parallel_dont_care() { rating = AMG_PARALLEL_DONT_CARE; }
	AMGNode() { rating = AMG_ASSIGNED_RATING; }
	int rating;
	inline bool is_parallel_dont_care() const { return rating == AMG_PARALLEL_DONT_CARE; }
	inline bool is_coarse() const {	return rating == AMG_COARSE_RATING;}
	inline bool is_dirichlet() const { return rating == AMG_DIRICHLET_RATING; }
	inline bool is_fine_direct() const {return (rating == AMG_FINE_RATING);}
	inline bool is_fine() const {return rating==AMG_FINE_RATING;}
	inline bool operator ==(int comp)
	{
		return comp == rating;
	}
	
	inline bool is_unassigned_fine_indirect() const {return rating == AMG_UNASSIGNED_FINE_INDIRECT_RATING;}
	inline bool is_aggressive_fine() const { return rating < AMG_FINE_RATING; }
	
	inline bool is_fine_indirect_level(int level) const { return rating == AMG_FINE_RATING+1-level;}
	

	inline bool is_assigned() const {return (rating <= AMG_ASSIGNED_RATING);}

	friend std::ostream &operator << (std::ostream &out, const AMGNode &n)
	{
		out << "Rating: " << n.rating;
		if(n.rating < 0)
		{
			if(n.is_coarse()) out << " (coarse)";
			else if (n.is_dirichlet()) out << " (dirichlet=fine) ";
			else if(n.is_unassigned_fine_indirect()) out << "(unknown indirect fine)";
			else if(n.is_parallel_dont_care()) out << "(parallel dont care)";
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
#ifdef UG_PARALLEL
	AMGNodes(ParallelNodes &_PN) : PN(_PN)
#else
	AMGNodes()
#endif
	{
		m_unassigned = 0;
		m_unassignedIndirectFine = 0;
		m_iNrOfCoarse = 0;
		m_iNrOfIndirectFine = 0;
	}

	#ifdef UG_PARALLEL
	AMGNodes(size_t N, ParallelNodes &_PN) : PN(_PN)
#else
	AMGNodes(size_t N)
#endif
	{
		m_unassigned = 0;
		m_unassignedIndirectFine = 0;
		m_iNrOfCoarse = 0;
		m_iNrOfIndirectFine = 0;
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

	void assign(size_t i)
	{
		if(!m_nodes[i].is_assigned())
		{
			UG_ASSERT(m_unassigned != 0, i);
			m_unassigned--;
		}
	}

	inline void external_set_coarse(size_t i)
	{
		assign(i);
		if(m_nodes[i].is_coarse()) return;
		m_nodes[i].set_coarse();
		m_iNrOfCoarse++;
	}

	inline void set_coarse(size_t i)
	{
		external_set_coarse(i);
	}
	inline void set_fine_direct(size_t i)
	{
		assign(i);
		if(m_nodes[i].is_coarse()) m_iNrOfCoarse--;
		m_nodes[i].set_fine_direct();
	}

	inline void set_fine(size_t i)
	{
		set_fine_direct(i);
	}

	inline void set_unassigned_fine_indirect(size_t i)
	{
		assign(i);
		if(m_nodes[i].is_coarse()) m_iNrOfCoarse--;
		m_unassignedIndirectFine++;
		m_nodes[i].set_unassigned_fine_indirect();
	}

	inline void set_fine_indirect_level(size_t i, int level)
	{
		UG_ASSERT(m_nodes[i].is_unassigned_fine_indirect(), i << " is not marked as indirect fine, but as " << m_nodes[i]);
		m_nodes[i].set_fine_indirect_level(level);
		m_unassignedIndirectFine--;
		m_iNrOfIndirectFine++;
	}

	inline void set_dirichlet(size_t i)
	{
		assign(i);
		if(m_nodes[i].is_coarse()) m_iNrOfCoarse--;
		m_nodes[i].set_dirichlet();
	}

	inline void set_parallel_dont_care(size_t i)
	{
		assign(i);
		m_nodes[i].set_parallel_dont_care();
	}


	size_t get_unassigned() const
	{
		return m_unassigned;
	}

	size_t get_unassigned_indirect_fine() const
	{
		return m_unassignedIndirectFine;
	}

	size_t get_nr_of_coarse() const
	{
		return m_iNrOfCoarse;
	}

	size_t get_nr_of_indirect_fine() const
	{
		return m_iNrOfIndirectFine;
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

	void print()
	{
		UG_LOG(std::endl);
		for(size_t i=0; i<size(); i++)
			UG_LOG(i << ": " << m_nodes[i] << std::endl);
		UG_LOG(std::endl);
	}




#ifdef UG_PARALLEL
public:
	bool needs_assignment(size_t i) const
	{
		return is_master_or_inner(i);
	}
	bool is_master_or_inner(size_t i) const
	{
		return PN.is_master_or_inner(i);
	}
	bool is_inner(size_t i) const
	{
		return PN.is_inner(i);
	}
	bool is_master(size_t i) const
	{
		return PN.is_master(i);
	}

	bool is_slave(size_t i) const
	{
		return PN.is_slave(i);
	}

#else
public:
	bool needs_assignment(size_t i) const { return true; }
	bool is_master_or_inner(size_t i) const { return true; }
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
#ifdef UG_PARALLEL
	ParallelNodes &PN;
#endif
	stdvector<AMGNode> m_nodes;
	size_t m_unassigned;
	size_t m_unassignedIndirectFine;
	size_t m_iNrOfCoarse;			//< nr of nodes to assign
	size_t m_iNrOfIndirectFine;
};

//	typedef maxheap<AMGNode> nodeinfo_pq_type;
typedef BoxPriorityQueue<AMGNode> nodeinfo_pq_type;

}


#endif // __H__UG__LIB_ALGEBRA__RSAMG_SOLVER__RSAMG_NODEINFO_H__
