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

#include "orientation_util.h"

namespace ug{

bool EdgeOrientationMatches(EdgeVertices* ev, Face* f)
{
	return OrientationMatches(ev, f);
}

bool OrientationMatches(EdgeVertices* ev, Face* f)
{
//	find the first vertex of ed in f
	size_t i;
	for(i = 0; i < f->num_vertices(); ++i)
	{
		if(f->vertex(i) == ev->vertex(0))
			break;
	}

	if(i < f->num_vertices())
	{
	//	the first one has been found.
	//	check whether the second vertex of ed is the
	//	same as the next vertex of f
		if(ev->vertex(1) == f->vertex((i+1)%f->num_vertices()))
			return true;//	the orientation is the same
	}

//	the orientation is not the same.
	return false;
}


UG_API 
bool OrientationMatches(FaceVertices* fv, Volume* v)
{
//	find the matching face desc and compare
	FaceDescriptor fd;
	for(size_t iface = 0; iface < v->num_faces(); ++iface){
		v->face_desc(iface, fd);
		if(CompareVertices(fv, &fd)){
		//	check if their orientation matches
		//	find the first vertex of fv in f
			size_t i;
			for(i = 0; i < fd.num_vertices(); ++i)
			{
				if(fd.vertex(i) == fv->vertex(0))
					break;
			}

			if(i < fd.num_vertices())
			{
			//	the first one has been found.
			//	check whether the second vertex of ed is the
			//	same as the next vertex of f
				if(fv->vertex(1) == fd.vertex((i+1) % fd.num_vertices()))
					return true;//	the orientation is the same
			}

		//	the orientation is not the same.
			return false;
		}
	}

//	the orientation is not the same.
	return false;
}

}//	end of namespace
