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


struct famg_nodeinfo
{
	int rating;
	//int newIndex;		
	inline void set_fine(){rating = FAMG_FINE_RATING;}
	inline void set_coarse(){rating = FAMG_COARSE_RATING;}
	inline void set_uninterpolateable(){rating = FAMG_UNINTERPOLATEABLE;	}
	
	inline void is_fine() { return rating == FAMG_FINE_RATING; }
	inline void is_coarse() { return rating == FAMG_COARSE_RATING; }
	inline void is_uninterpolateable() { return rating == FAMG_UNINTERPOLATEABLE; }
	
	inline void is_valid_rating()
	{
		return rating >= 0;
	}
	inline bool operator > (const amg_nodeinfo &other) const
	{
		if(rating == other.rating)
			return this < &other;
		else
			return rating > other.rating;
	}

	inline size_t get_val() const
	{
		UG_ASSERT(rating >= 0 && rating < 1000, "rating is " << rating << ", out of bounds [0, 1000]");
		return (int)rating;
	}
};

}


#endif // __H__LIB_DISCRETIZATION__FAMG_SOLVER__FAMG_NODEINFO_H__
