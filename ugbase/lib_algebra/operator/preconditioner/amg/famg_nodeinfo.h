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
#define FAMG_DIRICHLET_RATING		(-1003)


class famg_nodeinfo
{
public:
	famg_nodeinfo() { rating = 0.0; }
	double rating;
	//int newIndex;		
	inline void set_fine(){rating = FAMG_FINE_RATING;}
	inline void set_coarse(){rating = FAMG_COARSE_RATING;}
	inline void set_uninterpolateable(){rating = FAMG_UNINTERPOLATEABLE;	}
	inline void set_dirichlet() { rating = FAMG_DIRICHLET_RATING; }
	
	inline bool is_fine() const { return rating == FAMG_FINE_RATING; }
	inline bool is_coarse() const { return rating == FAMG_COARSE_RATING; }
	inline bool is_uninterpolateable() const { return rating == FAMG_UNINTERPOLATEABLE; }
	inline bool is_dirichlet() const { return rating == FAMG_DIRICHLET_RATING; }
	
	inline bool is_valid_rating()
	{
		return rating >= 0;
	}
	inline bool operator > (const famg_nodeinfo &other) const
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

	friend ostream &operator << (ostream &out, const famg_nodeinfo &n)
	{
		out << "Rating: " << n.rating;
		if(n.is_fine()) out << " (fine)";
		else if(n.is_coarse()) out << " (coarse)";
		if(n.is_uninterpolateable()) out << " (uninterpolateable)";
		return out;
	}
};

class famg_nodes
{
public:
	famg_nodes(size_t size)
	{
		nodes.resize(size);
	}

	famg_nodeinfo &operator [] (size_t i) { return nodes[i]; }
	const famg_nodeinfo &operator [] (size_t i) const { return nodes[i]; }

	template<typename vec_type>
	int get_rating(const vec_type &M)
	{
		size_t rating = 0;
		for(size_t i=0; i<M.size(); ++i)
		{
			const famg_nodeinfo &ninfo = nodes[M[i].from];
			if(ninfo.is_fine())
			{
				UG_DLOG(LIB_ALG_AMG, 2, " pair " << GetOriginalIndex(M[0].from) << ", " << GetOriginalIndex(M[1].from) << " is invalid, since " << GetOriginalIndex(M[i].from) << " is fine. ");
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
				swap(PN[i], PN.back()); // remove this pair
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
			if(mini != 0) swap(PN[0], PN[mini]);
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
				nodes[node].set_uninterpolateable();
				return true;
			}
		}
	}

	size_t size() const { return nodes.size(); }
private:
	stdvector<famg_nodeinfo> nodes; // !!! this HAS to be a consecutive array
};


template<typename neighborstruct>
void UpdateNeighbors(const cgraph &SymmNeighGraph, size_t node, stdvector<stdvector<neighborstruct> > &possible_neighbors,
		famg_nodes &nodes, maxheap<famg_nodeinfo> &heap)
{
	for(cgraph::cRowIterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
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
		famg_nodes &nodes, heap_type &heap)
{
	UG_DLOG(LIB_ALG_AMG, 2, "\nGetRatings...\n\n");
	for(size_t i=0; i<nodes.size(); i++)
	{
		UG_DLOG(LIB_ALG_AMG, 2, "node " << GetOriginalIndex(i) << ": ");
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
