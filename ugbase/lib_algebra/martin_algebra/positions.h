/*
 *  positions.h
 *  flexamg
 *
 *  Created by Martin Rupp on 03.03.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#pragma once
namespace ug{
////////////////////////////////////////////////////////////////////////////////
// stuff for position output
struct postype
{
	double x, y, z;
	friend std::ostream &operator << (std::ostream &out, const postype &p)
	{
		if(flexamg_dimensions == 2)
			out << p.x << " " << p.y;
		else
			out << p.x << " " << p.y << " " << p.z;
		return out;
	}
};


postype GetPosForIndex(int i);

void writePosToStream(ostream &out);
void writeToPosFile(const char *filename);

extern int *parentIndex[32];
static int GetOriginalIndex(int level, int i)
{
	while(level > 0)
		i = parentIndex[level--][i];
	return i;
}

static postype GetPosForIndexAtLevel(int i, int level)
{
	return GetPosForIndex(GetOriginalIndex(level, i));
}


} // namespace ug
