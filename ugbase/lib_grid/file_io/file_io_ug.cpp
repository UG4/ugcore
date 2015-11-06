/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Martin Stepniewski, Sebastian Reiter
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

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <set>
#include "file_io_ug.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/algorithms/polychain_util.h"
#include "lib_grid/callbacks/callbacks.h"


using namespace std;

namespace ug
{

static	bool CollectLines(Grid& grid, const SubsetHandler& shFace, EdgeSelector& LineSel);

static	bool CollectInnerVertices(Grid& grid, VertexSelector& InnVrtSel, VertexSelector& SurfVrtSel);

static	bool CollectSurfaceVertices(Grid& grid, const SubsetHandler& shFace,
									VertexSelector& SurfVrtSel);

static	bool CollectAllVerticesForNG(Grid& grid, VertexSelector& NgVrtSel,
							VertexSelector& SurfVrtSel, VertexSelector& InnVrtSel);

static bool WriteLGM(Grid& grid,
					  const char* lgmFilename,
					  const char* problemName,
					  const char* lgmName,
					  int convex,
					  const SubsetHandler& shFaces,
					  const SubsetHandler& shVolumes,
					  EdgeSelector& LineSel,
					  VertexSelector& SurfVrtSel,
					  Grid::EdgeAttachmentAccessor<AInt>& aaLineIndex,
					  Grid::VertexAttachmentAccessor<AInt>& aaSurfVrtIndex,
					  Grid::VertexAttachmentAccessor<AVector3>& aaPos);

static	bool GetRightLeftUnitIndex(int& rightIndex, int& leftIndex, Grid& grid, Face* face,
								   const SubsetHandler& shVolume);

static bool WriteNG(Grid& grid,
					 const SubsetHandler& shFaces,
					 const SubsetHandler& shVolumes,
					 const char* ngFilename,
					 VertexSelector& SurfVrtSel,
					 VertexSelector& InnVrtSel,
					 VertexSelector& NgVrtSel,
					 EdgeSelector& LineSel,
					 Grid::EdgeAttachmentAccessor<AInt>& aaLineIndex,
					 Grid::VertexAttachmentAccessor<AInt>& aaInnVrtIndex,
					 Grid::VertexAttachmentAccessor<AInt>& aaNgVrtIndex,
					 Grid::FaceAttachmentAccessor<AInt>& aaFaceIndex,
					 Grid::VertexAttachmentAccessor<AVector3>& aaPos);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	ConvertTETGENToUG
/// converts tetgen files (*.node, *.face and *.ele) to UG files (*.lgm, *.ng)
bool ExportGridToUG(const Grid& g, const SubsetHandler& shFace, const SubsetHandler& shVolume,
					const char* fileNamePrefix, const char* lgmName,
					const char* problemName, int convex)
{
//	the original grid may not be altered
	Grid grid = g;

//	we need subset-handlers that operate on the local grid
	SubsetHandler shFaces(grid, SHE_FACE);
	SubsetHandler shVolumes(grid, SHE_VOLUME);
	shFaces = shFace;
	shVolumes = shVolume;

//	fix orientation of faces
	for(int i = 0; i < shFaces.num_subsets(); ++i)
		FixFaceOrientation(grid, shFaces.begin<Face>(i), shFaces.end<Face>(i));

//	make sure that all face-subsets are in consecutive order
	{
		bool foundEmpty = false;
		for(int i = 0; i < shFaces.num_subsets(); ++i){
			if(shFaces.num<Face>(i) == 0)
				foundEmpty = true;
			else{
				if(foundEmpty){
				//	there must be no empty face-subsets between filled ones.
					UG_LOG("WARNING in ExportGridToUG: Empty face subset found between filled ones. Aborting...\n");
					return false;
				}
			}
		}
	}

//	make sure that all volume-subsets are in consecutive order
	{
		bool foundEmpty = false;
		for(int i = 0; i < shVolumes.num_subsets(); ++i){
			if(shVolumes.num<Volume>(i) == 0)
				foundEmpty = true;
			else{
				if(foundEmpty){
				//	there must be no empty volume-subsets between filled ones.
					UG_LOG("WARNING in ExportGridToUG: Empty volume subset found between filled ones. Aborting...\n");
					return false;
				}
			}
		}
	}

//	initialization
	EdgeSelector	LineSel(grid);
	VertexSelector 	NgVrtSel(grid);
	VertexSelector	SurfVrtSel(grid);
	VertexSelector 	InnVrtSel(grid);

	grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	grid.enable_options(EDGEOPT_STORE_ASSOCIATED_FACES);
//	for write ng
	grid.enable_options(FACEOPT_STORE_ASSOCIATED_EDGES);
	grid.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);

	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPosition);

//	selection of lines, surface-vertices, inner vertices
	CollectLines(grid, shFaces, LineSel);
	CollectSurfaceVertices(grid, shFaces, SurfVrtSel);
	CollectInnerVertices(grid, InnVrtSel, SurfVrtSel);
	CollectAllVerticesForNG(grid, NgVrtSel, SurfVrtSel, InnVrtSel);


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//	INDEX ASSIGNMENT SECTION >>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//	assign index to every line
		AInt aLineIndex;
		grid.attach_to_edges(aLineIndex);
		Grid::EdgeAttachmentAccessor<AInt> aaLineIndex(grid, aLineIndex);
		int counter = 0;
		for(EdgeIterator EIter = LineSel.begin(); EIter != LineSel.end(); ++EIter, ++counter)
			aaLineIndex[*EIter] = counter;

	//	assign index to every face regarding EACH surface
		AInt aFaceIndex;
		grid.attach_to_faces(aFaceIndex);
		Grid::FaceAttachmentAccessor<AInt> aaFaceIndex(grid, aFaceIndex);
		for(int i = 0; i < shFaces.num_subsets(); ++i)
		{
			counter = 0;
		//	assign triangle indices
			for(TriangleIterator iter = shFaces.begin<Triangle>(i);
				iter != shFaces.end<Triangle>(i); ++iter, ++counter)
				aaFaceIndex[*iter] = counter;

		//	assign quadrilaterals. Increase index by 2, since 2 triangles are written for each quad
			for(QuadrilateralIterator iter = shFaces.begin<Quadrilateral>(i);
				iter != shFaces.end<Quadrilateral>(i); ++iter, counter += 2)
				aaFaceIndex[*iter] = counter;
		}

	//	assign index to every grid-vertex in *.ng file order
		AInt aNgVrtIndex;
		grid.attach_to_vertices(aNgVrtIndex);
		Grid::VertexAttachmentAccessor<AInt> aaNgVrtIndex(grid, aNgVrtIndex);
		counter = 0;
		for(VertexIterator VIter = NgVrtSel.begin(); VIter != NgVrtSel.end(); ++VIter, ++counter)
			aaNgVrtIndex[*VIter] = counter;

	// 	assign index to every inner-vertex
		AInt aInnVrtIndex;
		grid.attach_to_vertices(aInnVrtIndex);
		Grid::VertexAttachmentAccessor<AInt> aaInnVrtIndex(grid, aInnVrtIndex);
		counter = 0;
		for(VertexIterator VIter = InnVrtSel.begin(); VIter != InnVrtSel.end(); ++VIter, ++counter)
			aaInnVrtIndex[*VIter] = counter;

	//	assign index to every surface-vertex
		AInt aSurfVrtIndex;
		grid.attach_to_vertices(aSurfVrtIndex);
		Grid::VertexAttachmentAccessor<AInt> aaSurfVrtIndex(grid, aSurfVrtIndex);
		counter = 0;
		for(VertexIterator SVIter = SurfVrtSel.begin(); SVIter != SurfVrtSel.end(); ++SVIter, ++counter)
			aaSurfVrtIndex[*SVIter] = counter;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	string lgmFilename(fileNamePrefix);
	lgmFilename.append(".lgm");
//	write *.lgm file
	WriteLGM(grid, lgmFilename.c_str(), problemName, lgmName,
			 convex, shFaces, shVolumes, LineSel, SurfVrtSel, aaLineIndex, aaSurfVrtIndex, aaPos);

	string ngFilename(fileNamePrefix);
	ngFilename.append(".ng");

//	write *.ng file
	WriteNG(grid, shFaces, shVolumes, ngFilename.c_str(), SurfVrtSel, InnVrtSel, NgVrtSel, LineSel,
			aaLineIndex, aaInnVrtIndex, aaNgVrtIndex, aaFaceIndex, aaPos);

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	CollectLines
/// collects lines in a selector
static bool CollectLines(Grid& grid, const SubsetHandler& shFace, EdgeSelector& LineSel)
{
//	store associated faces in this vector
	vector<Face*> vFaces;

//	iterate through all edges in the grid and identify lines by comparing the subset-membership
//	of the associated faces
	for(EdgeIterator EIter = grid.edges_begin(); EIter != grid.edges_end(); ++EIter)
	{
		CollectFaces(vFaces, grid, *EIter);

		//if(vFaces.size() > 1)
		{
			size_t i = 0;
			Face* f1 = NULL;
		//	find the first face that is assigned to a subset
			for(; i < vFaces.size(); ++i)
			{
				if(shFace.get_subset_index(vFaces[i]) != -1){
					f1 = vFaces[i];
				//	if the iteration leaves with a break, we have to increase i manually
					++i;
					break;
				}
			}
						
		//	compare with others. only check the ones that are assigned to a subset.
		//	if only one face is in a subset, we'll need a line anyway
			bool gotTwoSubsetFaces = false;
			for(; i < vFaces.size(); ++i)
			{
				if(shFace.get_subset_index(vFaces[i]) == -1)
					continue;

				gotTwoSubsetFaces = true;
				if(shFace.get_subset_index(f1) != shFace.get_subset_index(vFaces[i]))
				{
					LineSel.select(*EIter);
					break;
				}
			}
			
			if((!gotTwoSubsetFaces) && f1){
			//	the edge is adjacent to only one subset-face.
			//	it thus has to be a line
				LineSel.select(*EIter);
			}
		}
	}

	return true;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	CollectInnerVertices
/// collects inner-vertices in a selector
static bool CollectInnerVertices(Grid& grid, VertexSelector& InnVrtSel, VertexSelector& SurfVrtSel)
{
//	iterate through all grid-vertices and select them with VertexSelector InnVrtSel
	for(VertexIterator VIter = grid.vertices_begin(); VIter != grid.vertices_end(); ++VIter)
	{
		InnVrtSel.select(*VIter);
	}

//	iterate through all surface-vertices and deselect them in VertexSelector InnVrtSel, so that there only remain
//	the inner vertices
	for(VertexIterator VIter = SurfVrtSel.begin(); VIter != SurfVrtSel.end(); ++VIter)
	{
		InnVrtSel.deselect(*VIter);
	}

	return true;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	CollectSurfaceVertices
/// collects surface-vertices in a selector
static bool CollectSurfaceVertices(Grid& grid, const SubsetHandler& shFace,
									VertexSelector& SurfVrtSel)
{
//	iterate through all faces that are assigned to a subset and select them
	for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter)
	{
		Face* f = *iter;
		if(shFace.get_subset_index(f) != -1)
		{
			for(uint i = 0; i < f->num_vertices(); ++i)
				SurfVrtSel.select(f->vertex(i));
		}
	}

	return true;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	CollectAllVerticesForNG
/// collects all vertices in *.ng file order in a selector
static bool CollectAllVerticesForNG(Grid& grid, VertexSelector& NgVrtSel,
							VertexSelector& SurfVrtSel, VertexSelector& InnVrtSel)
{
//	collecting all vertices for an *.ng file requires first to select the surface-vertices and then the inner ones

//	iterate through all surface-vertices and select them
	for(VertexIterator VIter = SurfVrtSel.begin(); VIter != SurfVrtSel.end(); ++VIter)
	{
		NgVrtSel.select(*VIter);
	}

//	iterate through all inner vertices and select them
	for(VertexIterator VIter = InnVrtSel.begin(); VIter != InnVrtSel.end(); ++VIter)
	{
		NgVrtSel.select(*VIter);
	}

	return true;
}

/**	THIS METHOD USES Grid::mark
 * This method is required since double-triangles have to be avoided
 * during triangulation of quadrilaterals which share two edges whith
 * another quadrilateral (this e.g. occurs in fracture geometries with neumann
 * boundaries).
 * Note that only faces in the same subset as q are regarded.
 */
static int GetFirstRegularVertex(Grid& grid, const SubsetHandler& sh, Face* q)
{
	grid.begin_marking();
	vector<Edge*> edges;
	vector<Face*> faces;
	int si = sh.get_subset_index(q);

//	check whether a connected face exists, which shares two sides with q.
//	iterate over all edges, collect associated faces and mark them.
//	if a face was already marked, then it is shared by at least
//	two sides. We will then immediately mark all associated vertices of
//	that face.
	CollectAssociated(edges, grid, q);
	for(size_t i_edge = 0; i_edge != edges.size(); ++i_edge){
		CollectAssociated(faces, grid, edges[i_edge]);
	//	add all faces to nbrs. Check whether it is already contained.
		for(size_t i_face = 0; i_face < faces.size(); ++i_face){
			Face* f = faces[i_face];
			if(sh.get_subset_index(f) != si)
				continue;
			if(f == q)
				continue;

		//	check whether the face was already encountered
			if(grid.is_marked(f)){
			//	it was already there!
			//	mark all of its vertices
				for(size_t i_vrt = 0; i_vrt < f->num_vertices(); ++i_vrt)
					grid.mark(f->vertex(i_vrt));
			}
			else{
			//	its there for the first time. mark it.
				grid.mark(f);
			}
		}
	}

//	find the first unmarked vertex in q
	int firstUnmarked = -1;
	for(int i_vrt = 0; i_vrt < (int)q->num_vertices(); ++i_vrt){
		if(!grid.is_marked(q->vertex(i_vrt))){
			firstUnmarked = i_vrt;
			break;
		}
	}
	grid.end_marking();
	return firstUnmarked;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	WriteLGM
/// writes an *lgm file
static bool WriteLGM(Grid& grid,
			  const char* lgmFilename,
			  const char* problemName,
			  const char* lgmName,
			  int convex,
			  const SubsetHandler& shFaces,
			  const SubsetHandler& shVolumes,
			  EdgeSelector& LineSel,
			  VertexSelector& SurfVrtSel,
			  Grid::EdgeAttachmentAccessor<AInt>& aaLineIndex,
			  Grid::VertexAttachmentAccessor<AInt>& aaSurfVrtIndex,
			  Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
//	initialization
	ofstream out(lgmFilename);
	if(!out)
		return false;
	out.setf(ios::scientific);
	out.precision(24);

//	write the header
	out << "#Domain-Info" << endl;
	out << "name = " << lgmName << endl;
	out << "problemname = " << problemName << endl;
	out << "convex = " << convex << endl << endl;

//	write the units
	out << "#Unit-Info" << endl;
	for(int i = 0; i < shVolumes.num_subsets(); ++i)
	{
	//	if the name is empty then write a standard name
	//	be cautious with the unit indices:	id = 0 is reserved for outer space ONLY by UG
		if(shVolumes.subset_info(i).name.size() > 0)
			out << "unit " << i + 1 << " " << shVolumes.subset_info(i).name << endl;
		else {
			out << "unit " << i + 1 << " " << "unit_" << i + 1 << endl;
		}

	}
	out << endl;

//	write the lines
	out << "#Line-Info" << endl;
	for(EdgeIterator iter = LineSel.begin(); iter != LineSel.end(); ++iter)
	{
		Edge* e = *iter;
		out << "line " << aaLineIndex[e] << ": points: ";
		out << aaSurfVrtIndex[e->vertex(0)] << " " << aaSurfVrtIndex[e->vertex(1)] << ";" << endl;
	}
	out << endl;

//	write the surfaces
	out << "#Surface-Info" << endl;
	{
	//	used to find vertices that belong to the i-th subset.
		VertexSelector tmpSurfVrtSel(grid);
		for(int i = 0; i < shFaces.num_subsets(); ++i)
		{
			tmpSurfVrtSel.clear();

			//	identify left and right unit index of each surface
			int tmpLeft, tmpRight;
			if(!GetRightLeftUnitIndex(tmpRight, tmpLeft, grid, *shFaces.begin<Face>(i), shVolumes))
			{
				LOG("- GetRightLeftUnitIndex failed during lgm-write.\n");
				LOG("- In surface " << i << " face 0\n");
				LOG("- This can happen due to volume-elements with bad orinentation.\n");
				LOG("- IMPLEMENT a geometrical method for fallback!\n");
			}

			out << "surface " << i << ": left=" << tmpLeft << "; right=" << tmpRight << "; points:";
			for(ConstFaceIterator FIter = shFaces.begin<Face>(i);
				FIter != shFaces.end<Face>(i); ++FIter)
			{
				Face* f = *FIter;
				for(uint j = 0; j < f->num_vertices(); ++j)
				{
					Vertex* v = f->vertex(j);
					if(!tmpSurfVrtSel.is_selected(v))
					{
						tmpSurfVrtSel.select(v);
						out << " " << aaSurfVrtIndex[v];
					}
				}
			}

		//	write lines
		//	we have to avoid to write lines twice which lie inside a subset.
			grid.begin_marking();
			out << "; lines:";
			{
				vector<Edge*> vEdges;
				for(ConstFaceIterator FIter = shFaces.begin<Face>(i);
					FIter != shFaces.end<Face>(i); ++FIter)
				{
					CollectEdges(vEdges, grid, *FIter);
					for(vector<Edge*>::iterator EIter = vEdges.begin();
						EIter != vEdges.end(); ++EIter)
					{
						Edge* e = *EIter;
						if(LineSel.is_selected(e) && (!grid.is_marked(e)))
						{
							grid.mark(e);
							out << " " << aaLineIndex[e];
						}
					}
				}
			}
			grid.end_marking();
			
		// 	write triangles
			out << "; triangles:";
			{
				vector<Vertex*> vVertices;
				for(ConstTriangleIterator TIter = shFaces.begin<Triangle>(i);
					TIter != shFaces.end<Triangle>(i); ++TIter)
				{
					Triangle* t = *TIter;
					for(uint j = 0; j < t->num_vertices(); ++j)
					{
						out << " " << aaSurfVrtIndex[t->vertex(j)];
					}
					out << ";";
				}
			}

		//	write quadrilaterals as triangles
		//	special care has to be taken since quadrliaterals can possibly
		//	be connected at two sides. If one isn't careful one could thus produce
		//	two identical triangles during triangulation.
			{
			//	we need Grid::mark to check which point is not contained in a neighbour
				for(ConstQuadrilateralIterator iter = shFaces.begin<Quadrilateral>(i);
					iter != shFaces.end<Quadrilateral>(i); ++iter)
				{
					Quadrilateral* q = *iter;

					int firstRegular = GetFirstRegularVertex(grid, shFaces, q);

				//	if we didn't find one, we can't split the face. Schedule a warning.
					if(firstRegular == -1){
						UG_LOG("Can't split quadrilateral due to too many associated degenerated faces.\n");
					}
					else{
						for(int j = firstRegular; j < firstRegular + 3; ++j)
							out << " " << aaSurfVrtIndex[q->vertex(j%4)];

						out << ";";

						for(int j = firstRegular + 2; j < firstRegular + 5; ++j)
							out << " " << aaSurfVrtIndex[q->vertex(j%4)];

						out << ";";

					/*
						for(uint j = 0; j < 3; ++j)
							out << " " << aaSurfVrtIndex[q->vertex(j)];

						out << ";";

						for(uint j = 2; j < 5; ++j)
							out << " " << aaSurfVrtIndex[q->vertex(j%4)];

						out << ";";
					*/
					}
				}
			}
		//	surface written...
			out << endl;
		}
	}

//	write the vertices position data
	out << endl << "#Point-Info" << endl;
	for(VertexIterator iter = SurfVrtSel.begin(); iter != SurfVrtSel.end(); ++iter)
	{
		out << aaPos[*iter].x() << " " << aaPos[*iter].y() << " " << aaPos[*iter].z() << ";" << endl;
	}

//	close out file
	out.close();

	return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	GetRightLeftUnitIndices
/// identifies the right and left unit index to a given surface
static bool GetRightLeftUnitIndex(int& rightIndex, int& leftIndex, Grid& grid, Face* face,
									const SubsetHandler& shVolume)
{
//	initialization
	vector<Volume*> vVolumes;
	FaceDescriptor fd;

	rightIndex = 0;
	leftIndex = 0;

//	collect all volumes which are adjacent to face into vector<Volume*> vVolumes
	CollectVolumes(vVolumes, grid, face);
//	iterate through all volumes adjacent to the face and identify left and right unit index
	for(uint j = 0; j < vVolumes.size(); ++j)
	{
		Volume* v = vVolumes[j];

	//	find the face in the volume that matches the face
		for(uint i = 0; i < v->num_faces(); ++i)
		{
			v->face_desc(i, fd);
			if(CompareVertices(face, &fd))
			{
			//	we found a matching face.
			//	check whether the orientation is the same as in the face or not
				int i0 = GetVertexIndex(face, fd.vertex(0));
				int i1 = GetVertexIndex(face, fd.vertex(1));
				if(i1 == (i0 + 1) % (int)fd.num_vertices())
				{
					if(rightIndex)
						return false;
				//	the orientation is the same. the volume is on the right.
					rightIndex = shVolume.get_subset_index(v) + 1;
				}
				else
				{
					if(leftIndex)
						return false;
				//	the orientation is not the same. the volume is on the left.
					leftIndex = shVolume.get_subset_index(v) + 1;
				}

				break;
			}
		}

	//	stop iterating through all adjacent volumes, if indices have already been assigned
		if(leftIndex && rightIndex)
			break;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	WriteNG
/// writes an *ng file
static bool WriteNG(Grid& grid,
					 const SubsetHandler& shFaces,
					 const SubsetHandler& shVolumes,
					 const char* ngFilename,
					 VertexSelector& SurfVrtSel,
					 VertexSelector& InnVrtSel,
					 VertexSelector& NgVrtSel,
					 EdgeSelector& LineSel,
					 Grid::EdgeAttachmentAccessor<AInt>& aaLineIndex,
					 Grid::VertexAttachmentAccessor<AInt>& aaInnVrtIndex,
					 Grid::VertexAttachmentAccessor<AInt>& aaNgVrtIndex,
					 Grid::FaceAttachmentAccessor<AInt>& aaFaceIndex,
					 Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
//	initialization
	vector<Edge*> vEdges;
	vector<Face*> vFaces;
	vector<bool> faceFlags(shFaces.num_subsets());

	ofstream out(ngFilename);
	if(!out)
		return false;
	out.setf(ios::scientific);
	out.precision(24);
//UG_LOG("INFO: WRITING DEBUG INFO TO .ng FILE!\n");//ONLY FOR DEBUG
//size_t nodeCount = 0;//ONLY FOR DEBUG

//	write the boundary nodes section
	out << "# boundary nodes" << endl;
	for(VertexIterator VIter = SurfVrtSel.begin(); VIter != SurfVrtSel.end(); ++VIter)
	{
		Vertex* v = *VIter;

		for(int i = 0; i < shFaces.num_subsets(); ++i)
		{
			faceFlags[i] = false;
		}

//out << "# node-id: " << nodeCount << endl;//ONLY FOR DEBUG
//++nodeCount;//ONLY FOR DEBUG

	//	write the boundary-vertex position data
		out << "B " << aaPos[v].x() << " " << aaPos[v].y() << " " << aaPos[v].z() << endl;
		out << "	";

	//	iterate through all lines
		for(EdgeIterator EIter = LineSel.begin(); EIter != LineSel.end(); ++EIter)
		{
			Edge* e = *EIter;

		//	check if the current boundary-vertex is also a line-vertex
			if(GetVertexIndex(e, v) != -1)
			{
			//	write the line-index and line-vertex-index (corresponds to local position)
				out << " L " << aaLineIndex[e] << " " << GetVertexIndex(e, v);
			}
		}

	//	surface data
	//	iterate through all associated boundary-faces of vertex v
		for(Grid::AssociatedFaceIterator FIter = grid.associated_faces_begin(v);
			FIter != grid.associated_faces_end(v); ++FIter)
		{
			Face* f = *FIter;

		//	write the surface section including its index and the triangle index of
		//	ONE triangle-representative (per face-subset) associated to vertex v
			if(shFaces.get_subset_index(f) != -1)
			{
				if(faceFlags[shFaces.get_subset_index(f)] == false)
				{
				//	write local coordinates of vertex v on the triangle-representative
					int faceInd = aaFaceIndex[f];
					int vrtInd = GetVertexIndex(f, v);
					vector2 vCoord(0, 0);
					if(f->num_vertices() == 4){
						int firstRegular = GetFirstRegularVertex(grid, shFaces, f);
						if(firstRegular == -1){
							UG_LOG("Can't split quadrilateral due to too many associated degenerated faces.\n");
							continue;
						}
						vrtInd = (vrtInd + 4 - firstRegular)%4;
					}

					switch(vrtInd)
					{
						case 0:	vCoord = vector2(0, 0); break;
						case 1:	vCoord = vector2(1.0, 0); break;
						case 2:	vCoord = vector2(0, 1.0); break;
						case 3:	vCoord = vector2(1.0, 0); faceInd++; break; //index into the second sub-triangle
					}

					out << " S " << shFaces.get_subset_index(f) << " " << faceInd;
					out << " " << vCoord.x() << " " << vCoord.y();

					faceFlags[shFaces.get_subset_index(f)] = true;
				}
			}
		}

	out << ";" << endl;
	}

//	write inner nodes
	out << endl;
	out << "# inner nodes" << endl;

	if(InnVrtSel.num() > 0)
	{
	//	if there are inner vertices, iterate through all of them
		for(VertexIterator VIter = InnVrtSel.begin(); VIter != InnVrtSel.end(); ++VIter)
		{
			Vertex* v = *VIter;
//out << "# node-id: " << nodeCount << endl;//ONLY FOR DEBUG
//++nodeCount;//ONLY FOR DEBUG
		//	write the position data of the inner vertices
			out << "I " << aaPos[v].x() << " " << aaPos[v].y() << " " << aaPos[v].z() << ";" << endl;
		}
	}

//	write elements
	out << endl;
	out << "# elements" << endl;

	for(int i = 0; i < shVolumes.num_subsets(); ++i)
	{
	// 	iterate through all volumes
		for(ConstVolumeIterator VIter = shVolumes.begin<Volume>(i);
			VIter != shVolumes.end<Volume>(i); ++VIter)
		{
			Volume* v = *VIter;

		//	be cautious with the unit indices:	id = 0 is reserved for outer space ONLY by UG
			out << " E " << i + 1;

		//	iterate through all element vertices
			for(uint j = 0; j < v->num_vertices(); ++j)
			{
			//	write volume-vertex-index (in *.ng file order)
				out << " " << aaNgVrtIndex[v->vertex(j)];
			}

		//	collect all associated boundary-faces of volume v and write their *.ng file indices
			vFaces.clear();
			CollectFaces(vFaces, grid, v);
			for(vector<Face*>::iterator iter = vFaces.begin(); iter != vFaces.end(); ++iter)
			{
				Face* f = *iter;
				if(shFaces.get_subset_index(f) != -1)
				{
					out << " F " << aaNgVrtIndex[f->vertex(0)];
					for(uint j = 1; j < f->num_vertices(); ++j)
						out << " " << aaNgVrtIndex[f->vertex(j)];
				}
			}

			out << ";" << endl;
		}
	}

	out.close();
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	2D EXPORT
////////////////////////////////////////////////////////////////////////////////////////////////
//	ExportMeshToUG_2D and helper functions
static bool FaceIsOnRightSide(Face* f, Edge* e)
{
//  If the vertices in e are in the same order as the vertices in f,
//  then f is on the left side of e (since vertices are specified counter-clockwise).
    size_t ind1 = GetVertexIndex(f, e->vertex(0));
    size_t ind2 = GetVertexIndex(f, e->vertex(1));

    if(ind2 == (ind1 + 1) % f->num_vertices())
        return false;

    return true;
}

bool ExportGridToUG_2D(Grid& grid, const char* fileName, const char* lgmName,
					   const char* problemName, int convex,
					   SubsetHandler* psh)
{
	string lgmFileName = fileName;
	lgmFileName.append(".lgm");

	string ngFileName = fileName;
	ngFileName.append(".ng");

//	open the file
	ofstream out(lgmFileName.c_str());
	if(!out)
	{
		LOG("Failure in ExportGridToUG_2D: couldn't open " << lgmFileName << " for write" << endl);
		return false;
	}

	out.precision(12);

//	vectors are used to collect associated elements
	vector<Face*> vFaces;
	
//	Position accessor
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	initialise all indices with -1
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
	SetAttachmentValues(aaInt, grid.vertices_begin(), grid.vertices_end(), -1);

//	write the header
	out << "#Domain-Info" << endl;
	out << "name = " << lgmName << endl;
	out << "problemname = " << problemName << endl;
	out << "convex = " << convex << endl << endl;

	bool bUnitsSupplied;
	out << "#Unit-Info" << endl;

	if(psh)
	{
	//	there are units
		bUnitsSupplied = true;

	//	write the units
		bool bGotEmpty = false;// helps us to print a warning if required.
		for(int i = 0; i < psh->num_subsets(); i++)
		{
		//	units that do not contain faces are ignored.
			if(psh->num<Face>(i) == 0){
				bGotEmpty = true;
				continue;
			}

			if(bGotEmpty){
			//	schedule a warning
				UG_LOG("WARNING in ExportGridToUG_2D: Empty units between filled ones.\n");
			}

			SubsetInfo& sh = psh->subset_info(i);
			string unitName = sh.name;
			size_t k = sh.name.find(".");
			if(k != string::npos)
				unitName.erase(k, unitName.size() - k);
			if(unitName.size() > 0)
				out << "unit " << i+1 << " " << unitName.c_str() << endl;
			else
				out << "unit " << i+1 << " unit_" << i+1 << endl;
		}
	}
	else
	{
	//	no subsets are attached. This means we'll have to create one unit for all the triangles
		bUnitsSupplied = false;
		out << "unit 1 unit_1" << endl;
	}
	out << endl;

//	determine the vertices that go into the lgm file and assign indices to them.
//	This are all vertices wich lie on a border edge or wich are connected to triangles or quads of different subsets.
	{
		int pointIndexCounter = 0;

		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); iter++)
		{
			if(IsBoundaryVertex2D(grid, *iter))
			{
             	aaInt[*iter] = pointIndexCounter++;
			}
			else if(bUnitsSupplied)
			{
			//	check if the point is connected to triangles and quads of different subsets
                int subsetIndex = -1;
				
				CollectFaces(vFaces, grid, *iter);
				
				if(!vFaces.empty())
					subsetIndex = psh->get_subset_index(vFaces[0]);
				
				for(size_t i = 1; i < vFaces.size(); ++i){
					if(psh->get_subset_index(vFaces[i]) != subsetIndex){
						aaInt[*iter] = pointIndexCounter++;
						break;
					}
				}
			}
		}
	}

//	write the lines
	out << "#Line-Info" << endl;
//	we attach an int value to the edges wich describes their line-index. -1 if the edge is not a line
//	initialise line indices with -1
	grid.attach_to_edges(aInt);
	Grid::EdgeAttachmentAccessor<AInt> aaLineInt(grid, aInt);
	SetAttachmentValues(aaLineInt, grid.edges_begin(), grid.edges_end(), -1);

	int numLines = 0;
	vector<vector<Vertex*> > lines;
	{
	//	each edge which is connected to two different subsets or which is
	//	a boundary edge has to be written as a line
		if(bUnitsSupplied)
		{

		//	write the edge-subsets as lines
			for(int i = 0; i < psh->num_subsets(); ++i)
			{
			//	at this point we assume that lines are oriented
				if(!psh->empty<Edge>(i)){
				//	we have to follow the polygonal chain in subset i.
				//	We assume that the chain is regular
				//	create a callback for this subset, since we use it several times
					IsInSubset cbIsInSubset(*psh, i);
					std::pair<Vertex*, Edge*> curSec =
						GetFirstSectionOfPolyChain(grid, psh->begin<Edge>(i),
													psh->end<Edge>(i),
													cbIsInSubset);

					if(curSec.first == NULL)
						continue;

				//	find left and right unit
				//	curSec.second contains the first edge of the polychain.
					CollectFaces(vFaces, grid, curSec.second);

					if(vFaces.size() < 1)
						continue;

				//	we now have to check which subset lies on which side of the line
					int subLeft = psh->get_subset_index(vFaces[0]) + 1;
					int subRight = 0;
					if(vFaces.size() > 1)
						subRight = psh->get_subset_index(vFaces[1]) + 1;
					
					if(FaceIsOnRightSide(vFaces[0], curSec.second))
						swap(subLeft, subRight);

				//	since the edge and the first line segment can be flipped,
				//	we eventually have to swap again
					if(curSec.first != curSec.second->vertex(0))
						swap(subLeft, subRight);

					out << "line " << numLines << ": left="<< subLeft
						<< "; right=" << subRight << "; points:";
							
					size_t numEdgesInLine = 0;
					lines.push_back(vector<Vertex*>());
					vector<Vertex*>& curLine = lines.back();

					while(curSec.first)
					{
						Vertex* curVrt = curSec.first;
						Edge* curEdge = curSec.second;
						
					//	write the vertex index
						if(aaInt[curVrt] == -1){
							throw(UGError("RegularEdge subsets have not been assigned correctly."));
						}
						
						if(aaLineInt[curEdge] != -1){
							throw(UGError("Implementation error: One edge is in different lines."));
						}
						
						aaLineInt[curEdge] = numLines;
						++numEdgesInLine;
						
						out << " " << aaInt[curVrt];
						curLine.push_back(curVrt);
						
						std::pair<Vertex*, Edge*> nextSec =
							GetNextSectionOfPolyChain(grid, curSec, cbIsInSubset);

					//	check whether the next section is still valid.
					//	if not, write the last vertex and exit.
						if(!nextSec.first){
							out << " " << aaInt[GetConnectedVertex(curEdge, curVrt)];
							curLine.push_back(GetConnectedVertex(curEdge, curVrt));
							break;
						}
						else if(aaLineInt[nextSec.second] != -1){
						//	the edge has already been considered
							out << " " << aaInt[nextSec.first];
							curLine.push_back(nextSec.first);
							break;
						}
						
						curSec = nextSec;
					}
					out << ";" << endl;
					++numLines;
					
				//	make sure that all edges in the subset have been written to the line
					if(numEdgesInLine != psh->num<Edge>(i)){
						stringstream ss;
						ss  << "Couldn't write all edges of subset " << i << " to a line. "
							<< "This can happen if the subset is not a regular polychain. "
							<< "It is either separated into different parts or irregular (contains junctions).";
									  
						throw(UGError(ss.str().c_str()));
					}
				}
			}

/*
			for(EdgeIterator iter = grid.edges_begin(); iter != grid.edges_end(); iter++)
			{
				Edge* e = *iter;
			//	check for adjacent faces of different subset indices
				CollectFaces(vFaces, grid, e);
				
				if(vFaces.size() == 2)
				{
					if(psh->get_subset_index(vFaces[0]) != psh->get_subset_index(vFaces[1]))
					{
					//	e is a line.
					//	assign the index
						aaLineInt[e] = numLines;

					//	write the line
						int subLeft = psh->get_subset_index(vFaces[0]) + 1;
						int subRight = psh->get_subset_index(vFaces[1]) + 1;
						if(!FaceIsOnRightSide(vFaces[0], e))
							swap(subLeft, subRight);

						out << "line " << numLines << ": left="<< subLeft
							<< "; right=" << subRight << "; points: "
							<< aaInt[e->vertex(0)] << " " << aaInt[e->vertex(1)] << ";" << endl;

						numLines++;

                    //  endvertices of lines have to be marked as boundary-vertices
						assert((aaInt[e->vertex(0)] != -1) && "This point should be marked as boundary vertex!");
						assert((aaInt[e->vertex(1)] != -1) && "This point should be marked as boundary vertex!");
					}
				}
			}
*/
		}
		else{
		//	check for boundary lines
			for(EdgeIterator iter = grid.edges_begin(); iter != grid.edges_end(); iter++)
			{
				Edge* e = *iter;
				CollectFaces(vFaces, grid, e);
			//	check if it is a boundary edge
				if(vFaces.size() == 1)
				{
				//	assign the index
					aaLineInt[e] = numLines;

					int unitIndex = 1;
					if(bUnitsSupplied)
						unitIndex = psh->get_subset_index(vFaces[0]) + 1;

					int subLeft = 0;
					int subRight = 0;
					if(FaceIsOnRightSide(vFaces[0], e))
						subRight = unitIndex;
					else
						subLeft = unitIndex;

					out << "line " << numLines << ": left="<< subLeft
						<< "; right=" << subRight << "; points: "
						<< aaInt[e->vertex(0)] << " " << aaInt[e->vertex(1)] << ";" << endl;

					numLines++;

				//  endvertices of lines have to be marked as boundary-vertices
					assert((aaInt[e->vertex(0)] != -1) && "This point should be marked as boundary vertex!");
					assert((aaInt[e->vertex(1)] != -1) && "This point should be marked as boundary vertex!");
				}
			}
		}
	}

//	write the vertices position data
	{
		out << endl << "#Point-Info" << endl;
		for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++)
		{
		//	only write the point if it lies on a line
			if(aaInt[*iter] != -1)
				out << aaPos[*iter].x() << " " << aaPos[*iter].y() << ";" << endl;
		}
	}

//	lgm-write done
	out.close();

//	now we have to write the ng file
//	open the file
	out.open(ngFileName.c_str());

//	enable scientific number format�
	//out.setf(ios::scientific);

	if(!out)
	{
		LOG("Failure in ExportGridToUG_2D: couldn't open " << ngFileName << " for write" << endl);
		grid.detach_from_vertices(aInt);
		grid.detach_from_edges(aInt);
		return false;
	}

//	write the nodes
	{
		int numNGVertexs = 0;
	//	first we'll write the boundary vertices
		{
			for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++)
			{
				Vertex* p = *iter;
			//	check if p is a boundary node
				if(aaInt[p] != -1)
				{
				//	it is. write it to the file
					out << "B ";
				//	position:
					out << aaPos[p].x() << " " << aaPos[p].y();
				//	connected lines
				//	make sure that each is only used once
					set<int> encounteredLines;
					for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(p);
						eIter != grid.associated_edges_end(p); ++eIter)
					{
						if(aaLineInt[(*eIter)] != -1){
						//	get the local coordinate of the vertex on its associated line
							int lineIndex = aaLineInt[(*eIter)];
							if(encounteredLines.count(lineIndex) > 0)
								continue;// the line was already encountered.

							vector<Vertex*>& curLine = lines.at(lineIndex);
							assert(curLine.size() > 1 && "A line has to contain at least 2 vertices.");
							int localCoord = -1;
							for(size_t i = 0; i < curLine.size(); ++i){
								if(curLine[i] == p)
									localCoord = i;
							}
							if(localCoord >= 0){
								out << " L " << lineIndex << " " << localCoord;
								encounteredLines.insert(lineIndex);
							}
							else{
								UG_LOG("WARNING in ExportGridToUG_2D: "
										<< "Couldn't find local coordinate of boundary vertex at"
										<< aaPos[p] << endl);
							}
						}
					}
				//	done
					out << ";" << endl;
					numNGVertexs++;
				}
			}
		}

	//	write the inner vertices
		{
			for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); iter++)
			{
				Vertex* p = *iter;
				if(aaInt[p] == -1)
				{
				//	write the point
					out << "I " << aaPos[p].x() << " " << aaPos[p].y() << ";" << endl;
					aaInt[p] = numNGVertexs++;
				}
			}
		}
	}

//	write the elements
	{
	//	loop through all the triangles
		for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); iter++)
		{
			Face* f = *iter;

			int unitIndex = 1;
			if(bUnitsSupplied)
				unitIndex = psh->get_subset_index(f) + 1;

			out << "E " << unitIndex;
			for(size_t i = 0; i < f->num_vertices(); ++i)
				out << " " << aaInt[f->vertex(i)];

		//	check if one of its edges is a line.
			for(size_t i = 0; i < f->num_edges(); ++i)
			{
				Edge* e = grid.get_edge(f, i);
				if(aaLineInt[e] != -1)
				{
					out << " S " << aaInt[e->vertex(0)] << " " << aaInt[e->vertex(1)];
				}
			}
		//	done
			out << ";" << endl;
		}
	}

//	write complete
	out.close();

//	clean up
	grid.detach_from_vertices(aInt);
	grid.detach_from_edges(aInt);

	return true;
}

}//	end of namespace
