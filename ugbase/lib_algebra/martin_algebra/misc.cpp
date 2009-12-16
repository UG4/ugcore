/*
 *  misc.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */

#include "misc.h"
#include <iostream>
#include <fstream>

const char *boldredcolor = "\x1b[1;31m";
const char *boldgreencolor = "\x1b[1;32m";
const char *boldbluecolor = "\x1b[1;34m";

const char *redcolor = "\x1b[0;31m";
const char *greencolor = "\x1b[0;32m";
const char *bluecolor = "\x1b[0;34m";

const char *normalcolor = "\x1b[0;0m";

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

void spaceout(int n)
{
	for(int i=0;i<n; i++)
		cout << " ";
}

int *parentIndex[32];

pos2d *positions;
int iNrOfPositions;
pos2d GetPosForIndex(int i)
{
	return positions[i];	
}
void writePosToStream(ostream &out)
{
	out << iNrOfPositions << endl;
	for(int i=0; i<iNrOfPositions; i++)
		out << GetPosForIndex(i).x << " " << GetPosForIndex(i).y << endl;
}
void writeToPosFile(const char *filename)
{
	fstream file(filename, ios::out);
	writePosToStream(file);
}
