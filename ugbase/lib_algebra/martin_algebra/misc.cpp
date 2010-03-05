/*
 *  misc.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */

#include "misc.h"
#include "positions.h"
#include <iostream>
#include <fstream>
#include <vector>

namespace ug{
string nrstring(double d)
{
	char s[255];
	sprintf(s, "%g", d);
	return string(s);
}

string nrstring(int i)
{
	char s[255];
	sprintf(s, "%d", i);
	return string(s);
}


int *parentIndex[32];

std::vector<postype> positions;
postype GetPosForIndex(int i)
{
	return positions[i];	
}

void writeToPosFile(const char *filename)
{
	fstream file(filename, ios::out);
	writePosToStream(file);
}

void writePosToStream(ostream &out)
{
	out << positions.size() << endl;
	for(int i=0; i< positions.size() ; i++)
		out << positions[i] << endl;
}

} // namespace ug