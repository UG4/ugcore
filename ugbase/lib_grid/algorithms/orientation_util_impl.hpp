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

#ifndef __H__UG_orientation_util_impl
#define __H__UG_orientation_util_impl


#include <stack>
#include <vector>

#include "common/math/ugmath.h"
#include "geom_obj_util/face_util.h"
#include "geom_obj_util/volume_util.h"

namespace ug {


////////////////////////////////////////////////////////////////////////
//	InvertOrientation
template <typename iter_t>
void InvertOrientation(Grid& grid, iter_t elemsBegin,
					   iter_t elemsEnd)
{
	for(iter_t iter = elemsBegin; iter != elemsEnd; ++iter)
		grid.flip_orientation(*iter);
}


////////////////////////////////////////////////////////////////////////
//	FixFaceOrientation
template <typename TFaceIterator>
void FixFaceOrientation(Grid& grid, TFaceIterator facesBegin,
						TFaceIterator facesEnd)
{
	using namespace std;

	{
	//	make sure that boundary faces are oriented outwards
		Grid::volume_traits::secure_container vols;
		for(TFaceIterator iter = facesBegin; iter != facesEnd; ++iter){
			grid.associated_elements(vols, *iter);
			if(vols.size() == 1 && !OrientationMatches(*iter, vols[0])){
				grid.flip_orientation(*iter);
			}
		}
	}

//	the rest only considers manifold faces
//	we use marks to differentiate between processed and unprocessed faces
	grid.begin_marking();

//	we have to mark all faces between facesBegin and facesEnd initially,
//	to differentiate them from the other faces in the grid.
	for(TFaceIterator iter = facesBegin; iter != facesEnd; ++iter){
		if(NumAssociatedVolumes(grid, *iter) == 0)
			grid.mark(*iter);
	}

//	this edge descriptor will be used multiple times
	EdgeDescriptor ed;
	
//	containers to store neighbour elements.
	vector<Face*> vNeighbours;
	
//	stack that stores candidates
	stack<Face*> stkFaces;
	
//	we'll iterate through all faces
	while(facesBegin != facesEnd)
	{
	//	candidates are empty at this point.
	//	if the face is unprocessed it is a new candidate
		if(grid.is_marked(*facesBegin))
		{
		//	mark it as candidate (by removing the mark)
			grid.unmark(*facesBegin);
			stkFaces.push(*facesBegin);
			
		//	while the stack is not empty
			while(!stkFaces.empty())
			{
			//	get the candidate
				Face* f = stkFaces.top();
				stkFaces.pop();

			//	get the neighbours for each side
				for(size_t i = 0; i < f->num_edges(); ++i)
				{
					f->edge_desc(i, ed);
					GetNeighbours(vNeighbours, grid, f, i);

				//	fix orientation of unprocessed neighbours.
					for(size_t j = 0; j < vNeighbours.size(); ++j)
					{
						Face* fn = vNeighbours[j];
						if(grid.is_marked(fn))
						{
						//	check whether the orientation of f and fn differs.
							if(EdgeOrientationMatches(&ed, fn))
							{
							//	the orientation of ed is the same as the orientation
							//	of an edge in fn.
							//	the faces thus have different orientation.
								grid.flip_orientation(fn);
							}

						//	mark the face as processed and add it to the stack
							grid.unmark(fn);
							stkFaces.push(fn);
						}
					}
				}
			}
		}
		
	//	check the next face
		++facesBegin;
	}

	grid.end_marking();
}


template <typename TAAPosVRT>
bool
CheckOrientation(Volume* vol, TAAPosVRT& aaPosVRT)
{
//	some type definitions
	using vector_t = typename TAAPosVRT::ValueType;
	
//	First calculate the center of the volume
	vector_t volCenter = CalculateCenter(vol, aaPosVRT);
	
//	now check for each side whether it points away from the center.
	size_t numFaces = vol->num_faces();
	FaceDescriptor fd;
	vector_t normal;
	for(size_t i = 0; i < numFaces; ++i){
		vol->face_desc(i, fd);
		CalculateNormal(normal, &fd, aaPosVRT);
		
	//	in order to best approximate quadrilateral faces, we'll calculate the
	//	center of the face and compare that to the volCenter.
		vector_t faceCenter = CalculateCenter(&fd, aaPosVRT);
		
	//	now compare normal and center
		vector_t dir;
		VecSubtract(dir, faceCenter, volCenter);
		if(VecDot(dir, normal) < 0)
			return false;
	}
	
//	all center / normal checks succeeded. Orientation is fine.
	return true;
}

template <typename TVolIterator, typename TAAPosVRT>
int
FixOrientation(Grid& grid, TVolIterator volsBegin, TVolIterator volsEnd,
			   TAAPosVRT& aaPosVRT)
{
	int numFlips = 0;
//	iterate through all volumes
	for(TVolIterator iter = volsBegin; iter != volsEnd; ++iter){
	//	check whether the orientation is fine
		if(!CheckOrientation(*iter, aaPosVRT)){
			grid.flip_orientation(*iter);
			++numFlips;
		}
	}
	
	return numFlips;
}


}//	end of namespace

#endif