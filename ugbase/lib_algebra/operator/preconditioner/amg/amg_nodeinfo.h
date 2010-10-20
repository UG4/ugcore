/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 16.06.10
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_NODEINFO_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_NODEINFO_H__

//#include "maxheap.h"
#include "boxsort.h"

namespace ug {

//  structs
#define ASSIGNED_RATING					(-1000000000)
#define COARSE_RATING					(-1000000001)
#define FINE_RATING						(-2000000000)
#define FINE_RATING_INDIRECT_UNASSIGNED	(-1000000002)
struct amg_nodeinfo
{
	int rating;
	//int newIndex;		
	inline void setAssigned(){rating = ASSIGNED_RATING;}		
	inline void setCoarse()	{rating = COARSE_RATING;}
	inline void setFineDirect(){rating = FINE_RATING;	}
	
	inline void setFineIndirect(){rating = FINE_RATING_INDIRECT_UNASSIGNED;}
	
	inline bool isCoarse() const {	return rating == COARSE_RATING;}
	inline bool isFineDirect() const {return (rating == FINE_RATING);}
	
	inline bool isUnassignedFineIndirect() const {return rating == FINE_RATING_INDIRECT_UNASSIGNED;}
	
	
	inline bool isFineIndirectLevel(int level) const { return rating == FINE_RATING+1-level;}
	inline void setFineIndirectLevel(int level) { rating = FINE_RATING+1-level;}
	
	inline bool isAssigned() const {return (rating <= ASSIGNED_RATING);}
	
	friend ostream &operator << (ostream &out, const amg_nodeinfo &n)
	{
		out << "Rating: " << n.rating;
		if(n.rating < 0)
		{
			if(n.isCoarse()) out << " (coarse)";
			else if(n.isUnassignedFineIndirect()) out << "(unknown indirect fine)";
			else if(FINE_RATING-n.rating == 0)
				out << " (fine direct)";
			else
				out << " (indirect level " << FINE_RATING+1-n.rating << ")";
		}
		out << " ";
		return out;
	}
	void print() const
	{
		cout << *this << endl;
	} // << " newindex: " << newIndex << endl;
	
	inline bool operator > (const amg_nodeinfo &other) const
	{
		if(rating == other.rating)
			return this < &other; // we somehow want a STABLE sort, for that coarsening is in the direction of the numbering of the elements
		else
			return rating > other.rating;
	}

	inline size_t get_val() const
	{
		UG_ASSERT(rating > 0 && rating < 1000, "rating is " << rating << ", out of bounds [0, 1000]");
		return (int)rating;
	}
};

//	typedef maxheap<amg_nodeinfo> nodeinfo_pq_type;
typedef BoxSort<amg_nodeinfo> nodeinfo_pq_type;

}


#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_NODEINFO_H__
