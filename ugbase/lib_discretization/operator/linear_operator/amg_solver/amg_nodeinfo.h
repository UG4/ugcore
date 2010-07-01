/*
 *  amg_nodeinfo.h
 *  flexamg
 *
 *  Created by Martin Rupp on 16.06.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#pragma once

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
	
	inline bool isCoarse(){	return rating == COARSE_RATING;}
	inline bool isFineDirect()	{return (rating == FINE_RATING);}
	
	inline bool isUnassignedFineIndirect(){return rating == FINE_RATING_INDIRECT_UNASSIGNED;}
	
	
	inline bool isFineIndirectLevel(int level){ return rating == FINE_RATING+1-level;}
	inline void setFineIndirectLevel(int level){ rating = FINE_RATING+1-level;}
	
	inline bool isAssigned(){return (rating <= ASSIGNED_RATING);}
	
	friend ostream &operator << (ostream &out, amg_nodeinfo &n)
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
	void print()
	{
		cout << *this << endl;
	} // << " newindex: " << newIndex << endl;
	
	inline bool operator > (const amg_nodeinfo &other)
	{
		if(rating == other.rating)
			return this < &other; // we somehow want a STABLE sort, for that coarsening is in the direction of the numbering of the elements
		else
			return rating > other.rating;
	}
};


}