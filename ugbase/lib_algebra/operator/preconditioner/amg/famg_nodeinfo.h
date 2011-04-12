/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__LIB_DISCRETIZATION__FAMG_SOLVER__FAMG_NODEINFO_H__
#define __H__LIB_DISCRETIZATION__FAMG_SOLVER__FAMG_NODEINFO_H__

namespace ug {

//  structs
#define FAMG_UNINTERPOLATEABLE	(-1001)
#define FAMG_FINE_RATING		(-1002)
#define FAMG_COARSE_RATING		(-1003)
#define FAMG_DIRICHLET_RATING	(-1003)


class FAMGNode
{
	friend class FAMGNodes;

private:
	inline void set_fine(){rating = FAMG_FINE_RATING;}
	inline void set_coarse(){rating = FAMG_COARSE_RATING;}
	inline void set_uninterpolateable()
	{	rating = FAMG_UNINTERPOLATEABLE;	}
	inline void set_dirichlet() { rating = FAMG_DIRICHLET_RATING; }


public:
	FAMGNode() { rating = 0.0; }
	double rating;
	//int newIndex;		

	
	inline bool is_fine() const { return rating == FAMG_FINE_RATING; }
	inline bool is_coarse() const { return rating == FAMG_COARSE_RATING; }
	inline bool is_uninterpolateable() const { return rating == FAMG_UNINTERPOLATEABLE; }
	inline bool is_dirichlet() const { return rating == FAMG_DIRICHLET_RATING; }
	
	inline bool could_be_coarse() const
	{
		return is_valid_rating() || is_coarse() || is_uninterpolateable(); // is_dirichlet ??
	}
	inline bool is_valid_rating() const
	{
		return rating >= 0;
	}
	inline bool operator > (const FAMGNode &other) const
	{
		if(rating == other.rating)
			return this < &other;
		else
			return rating < other.rating;
	}

	inline double get_val() const
	{
		UG_ASSERT(rating >= 0 && rating < 1000, "rating is " << rating << ", out of bounds [0, 1000]");
		return rating;
	}

	friend std::ostream &operator << (std::ostream &out, const FAMGNode &n)
	{
		out << "Rating: " << n.rating;
		if(n.is_fine()) out << " (fine)";
		else if(n.is_coarse()) out << " (coarse)";
		if(n.is_uninterpolateable()) out << " (uninterpolateable)";
		return out;
	}
};

//! this structure holds coarse/fine information as well as master/slave
//! for use in amg. it tries to replace the missing "node" object here
class FAMGNodes
{
public:
	FAMGNodes(SparseMatrix<double> &_P, size_t level, ug::cAMG_helper &amghelper)
	: P(_P), m_level(level), m_amghelper(amghelper)
	{
		m_iNrOfCoarse = 0;
	}

	void create(size_t size)
	{
		m_iNrOfCoarse = 0;
		nodes.clear(); 		nodes.resize(size);
		newIndex.clear(); 	newIndex.resize(size, -1);
#ifdef UG_PARALLEL
		OLtype.clear();		OLtype.resize(size, 0);
#endif
	}

	FAMGNode &operator [] (size_t i) { return nodes[i]; }
	const FAMGNode &operator [] (size_t i) const { return nodes[i]; }

	template<typename vec_type>
	int get_rating(const vec_type &M)
	{
		AMG_PROFILE_FUNC();
		size_t rating = 0;
		for(size_t i=0; i<M.size(); ++i)
		{
			const FAMGNode &ninfo = nodes[M[i].from];
			if(ninfo.is_fine())
			{
				UG_DLOG(LIB_ALG_AMG, 2, " pair " << get_original_index(M[0].from) << ", " << get_original_index(M[1].from) << " is invalid, since " << get_original_index(M[i].from) << " is fine. ");
				return -1;
			}
			/*else if(ninfo.is_uninterpolateable())
				rating += 1000;
			else*/
			if(!ninfo.is_coarse()) rating++;
		}
		return rating;
	}

	// updates the rating of node <node>
	// and writes it to rating
	// possible neighbors are always stored as follows:
	//	- possible_neighbors[i][0] is the best currently known pair for interpolating i
	//	- this can change
	//	a) if a parent node gets fine -> removal of this pair
	//	b) if a parent node gets coarse
	//	returns true if value has been updated, otherwise false
	template<typename neighborstruct>
	bool update_rating(size_t node, stdvector<neighborstruct> &PN)
	{
		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 2, " update rating of node " << node << "... ");
		if(nodes[node].is_valid_rating() == false)
		{
			UG_DLOG(LIB_ALG_AMG, 2, nodes[node].rating << " is not a valid rating."); return false;
		}

		int mini = -1;
		double minrating = 10000;
		//double minF = 1e12;
		UG_DLOG(LIB_ALG_AMG, 2, "has " << PN.size() << " possible parent nodes. ");
		for(size_t i=0; i<PN.size(); )
		{
			double irating = get_rating(PN[i].parents) + PN[i].F/1000;

			if(irating < 0)
			{
				std::swap(PN[i], PN.back()); // remove this pair
				PN.resize(PN.size()-1);
				UG_DLOG(LIB_ALG_AMG, 2, " removed pair " << i << " ");
				continue;
			}

			// choose the pairs with best minimization F.
			if(irating < minrating) // || (irating == minrating && PN[i].F < minF))
			{
				minrating = irating;
				//minF = PN[i].F;
				mini = i;
			}
			i++;
		}
		UG_DLOG(LIB_ALG_AMG, 2, "now has " << PN.size() << " possible parent nodes. ");
		for(size_t i=0; i<PN.size(); i++)
			UG_DLOG(LIB_ALG_AMG, 2, PN[i].parents[0].from << "-" << PN[i].parents[1].from << " ");

		if(mini != -1)
		{
			UG_DLOG(LIB_ALG_AMG, 2, " has rating " << minrating << "\n");
			if(mini != 0) std::swap(PN[0], PN[mini]);
			if(nodes[node].rating != minrating)
			{
				UG_DLOG(LIB_ALG_AMG, 2, " new rating! ");
				nodes[node].rating = minrating;
				return true;
			}
			else return false;
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 2, " is uninterpolateable! ");
			if(nodes[node].is_uninterpolateable())
				return false;
			else
			{
				set_uninterpolateable(node);
				return true;
			}
		}
	}


	size_t size() const { return nodes.size(); }



#ifdef UG_PARALLEL
public:
	bool is_inner_node(size_t i) const
	{
		return OLtype[i] == 0;
	}

	bool is_master(size_t i) const
	{
		return OLtype[i] & 1;
	}

	bool is_slave(size_t i, int OLlevel=0) const
	{
		return OLtype[i] & (1 << (OLlevel+1));
	}

	bool is_master_on(size_t i, int pid)
	{
		if(is_slave(i, 0))
			return m_masterOn[i] == pid;
		else
			return false;
	}

	bool i_must_assign(size_t i) const
	{
		return is_inner_node(i) || is_master(i);
	}

	bool i_can_set_coarse(size_t i) const
	{
		return nodes[i].is_coarse() ||	// node already coarse
			((i_must_assign(i) || is_slave(i, 0)) && nodes[i].could_be_coarse());  // or i can set coarse
	}

	void print_OL_types() const
	{
		for(size_t i=0; i<size(); i++)
		{
			UG_LOG("Index " << i << ": ");
			if(is_inner_node(i))
				UG_LOG("inner")
			else
			{
				UG_LOG("non-inner ");
				if(is_master(i))
					UG_LOG("master ")
				for(size_t j=0; j<5; j++)
					if(is_slave(i, j)) UG_LOG("slave" << j << " ");
			}
			UG_LOG("\n");
		}
	}

	void calculate_unassigned()
	{
		m_iUnassigned = 0;
		for(size_t i=0; i<size(); i++)
		{
			if(i_must_assign(i))
				m_iUnassigned++;

		}
	}
	stdvector<int> OLtype;
	std::map<size_t, int> m_masterOn;
#else
public:
	bool is_inner_node(size_t i)
	{
		return true;
	}

	bool is_master(size_t i)
	{
		return true;
	}

	bool is_master_on(size_t i, int pid)
	{
		return pid == 0;
	}

	bool is_slave(size_t i, int OLlevel=0)
	{
		return false;
	}

	bool i_must_assign(size_t i)
	{
		return true;
	}

	bool i_can_set_coarse(size_t i)
	{
		return nodes[i].is_valid_rating() || nodes[i].is_coarse();
	}
	void calculate_unassigned()
	{
		m_iUnassigned = size();
	}
#endif

	void set_fine(size_t index)
	{
		UG_ASSERT(nodes[index].is_coarse() == false, "try to set node " << index << " fine, but is coarse???");
		if(nodes[index].is_fine() == false && i_must_assign(index))
			m_iUnassigned--;
		nodes[index].set_fine();
	}

	void external_set_coarse(size_t index)
	{
		UG_ASSERT(i_can_set_coarse(index), "Index " << index << " can not be coarse.\n");
		if(nodes[index].is_coarse() == false)
		{
			nodes[index].set_coarse();

			m_iUnassigned--;
			// todo: perhaps we should not do this here
			newIndex[index] = m_iNrOfCoarse++;
			P(index, newIndex[index]) = 1.0;
		}
	}

	void set_coarse(size_t index)
	{
		UG_ASSERT(i_can_set_coarse(index), "Index " << index << " can not be coarse.\n");
		external_set_coarse(index);
	}

	void set_uninterpolateable(size_t index)
	{
		UG_ASSERT(nodes[index].is_valid_rating() || nodes[index].is_uninterpolateable(),
				"node " << index << " is already set to " << nodes[index].rating);
		nodes[index].set_uninterpolateable();
	}

	size_t get_nr_of_coarse()
	{
		return m_iNrOfCoarse;
	}

	size_t get_unassigned()
	{
		return m_iUnassigned;
	}

	size_t get_original_index(size_t i)
	{
		UG_ASSERT(m_level != (size_t)-1, "");
		return m_amghelper.GetOriginalIndex(m_level, i);
	}

//private:
public:
	stdvector<int> newIndex;
	stdvector<FAMGNode> nodes; // !!! this HAS to be a consecutive array, because it is used in the heap

private:
	size_t m_iNrOfCoarse;			// number of coarse nodes so far
	size_t m_iUnassigned;			// number of still unassigned nodes

	SparseMatrix<double> &P;
	size_t m_level;
	ug::cAMG_helper &m_amghelper;
};


template<typename neighborstruct>
void UpdateNeighbors(const cgraph &SymmNeighGraph, size_t node, stdvector<stdvector<neighborstruct> > &possible_neighbors,
		FAMGNodes &nodes, maxheap<FAMGNode> &heap)
{
	for(cgraph::const_row_iterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
	{
		size_t neigh = (*conn);
		UG_ASSERT(node != neigh, "");
		if(!nodes[neigh].is_valid_rating())
			continue;

		if(nodes.update_rating(neigh, possible_neighbors[neigh]))
		{
			if(nodes[neigh].is_uninterpolateable())
			{
				UG_DLOG(LIB_ALG_AMG, 2, " remove from heap! ");
				heap.remove(neigh);
			}
			else
			{
				UG_DLOG(LIB_ALG_AMG, 2, " update in heap! ");
				heap.update(neigh);
			}
		}
	}
}

template<typename neighborstruct, typename heap_type>
void GetRatings(stdvector<stdvector<neighborstruct> > &possible_neighbors,
		FAMGNodes &nodes, heap_type &heap)
{
	UG_DLOG(LIB_ALG_AMG, 2, "\nGetRatings...\n\n");
	for(size_t i=0; i<nodes.size(); i++)
	{
		UG_DLOG(LIB_ALG_AMG, 2, "node " << nodes.get_original_index(i) << ": ");
		if(nodes[i].rating == 0)
		{
			nodes.update_rating(i, possible_neighbors[i]);
			if(nodes[i].is_valid_rating())
				heap.insert_item(i);
		}
		UG_DLOG(LIB_ALG_AMG, 2, "rating = " << nodes[i].rating << "\n");
	}
}


}



#endif // __H__LIB_DISCRETIZATION__FAMG_SOLVER__FAMG_NODEINFO_H__
