/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "common/math/ugmath.h"
#include "subset_color_util.h"

namespace ug {

vector3 GetColorFromDefaultPalette(int index)
{
//	values taken from http://en.wikipedia.org/wiki/Web_colors
	float stdColors[][3] = {{150, 150, 255},//My blue
							{255, 0, 0},	//Red
							{0, 255, 0},	//Lime
							{0, 0, 255},	//Blue
							{255, 0, 255},	//Magenta
							{255, 255, 0},	//Yellow
							{0, 255, 255},	//Aqua
							{255, 192, 203},//Pink
							{152, 251, 152},//PaleGreen
							{176, 224, 230},//PowderBlue
							{240, 230, 140},//Khaki
							{255, 99, 71},	//Tomato
							{0, 191, 255},	//DeepSkyBlue
							{255, 160, 122}	//LightSalmon
							};

	const int numCols = 14;

	if(index >= 0 && index < numCols)
		return vector3(stdColors[index][0] / 255.f, stdColors[index][1] / 255.f, stdColors[index][2] / 255.f);

	index -= numCols;

	float val = 2.f* 3.14159265 * (float)index / 3.148 + (float)index / 15.f;
	vector3 vCol(1.f + cos(val), 1.f + sin(0.6* val), 1.f - cos(0.373*val));

	VecNormalize(vCol, vCol);
	return vCol;
}

////////////////////////////////////////////////////////////////////////
//	AssignSubsetColors
void AssignDefaultSubsetColors(ISubsetHandler& sh)
{
	for(int i = 0; i < sh.num_subsets(); ++i)
	{
		SubsetInfo& si = sh.subset_info(i);
		vector3 col = GetColorFromDefaultPalette(i);
		si.color.x() = col.x();
		si.color.y() = col.y();
		si.color.z() = col.z();
		si.color.w() = 1.f;
	}
}


void AssignSubsetColorsRedToGreen(ISubsetHandler& sh, int firstSi, int numSi)
{
	const vector3 red(1.0, 0, 0);
	const vector3 green(0, 1.0, 0);

	if(firstSi < 0)
		firstSi = 0;
	if((numSi < 0) || (firstSi + numSi > sh.num_subsets()))
	   numSi = sh.num_subsets() - firstSi;

	for(int i = 0; i < numSi; ++i)
	{
		number ia = 1;
		if(numSi > 1)
			ia = (number)i / (number)(numSi - 1);

		vector3 c;
		VecScaleAdd(c, (1. - ia), red, ia, green);
		
		SubsetInfo& si = sh.subset_info(i);
		si.color.x() = c.x();
		si.color.y() = c.y();
		si.color.z() = c.z();
		si.color.w() = 1.f;
	}
}

void AssignSubsetColorsBlueToGreen(ISubsetHandler& sh, int firstSi, int numSi)
{
	const vector3 blue(0, 0, 1.0);
	const vector3 green(0, 1.0, 0);

	if(firstSi < 0)
		firstSi = 0;
	if((numSi < 0) || (firstSi + numSi > sh.num_subsets()))
	   numSi = sh.num_subsets() - firstSi;

	for(int i = 0; i < numSi; ++i)
	{
		number ia = 1;
		if(numSi > 1)
			ia = (number)i / (number)(numSi - 1);

		vector3 c;
		VecScaleAdd(c, (1. - ia), blue, ia, green);
		
		SubsetInfo& si = sh.subset_info(i);
		si.color.x() = c.x();
		si.color.y() = c.y();
		si.color.z() = c.z();
		si.color.w() = 1.f;
	}
}

}//	end of namespace
