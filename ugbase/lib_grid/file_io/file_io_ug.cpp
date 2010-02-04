// created by Martin Stepniewski, Sebastian Reiter
// s.b.reiter@googlemail.com
//	y09 m05 d27

#include "file_io_ug.h"
#include "lib_grid/lib_grid.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>

using namespace std;

namespace ug
{

static	bool CollectLines(Grid& grid, SubsetHandler& shFace, EdgeSelector& LineSel);

static	bool CollectInnerVertices(Grid& grid, VertexSelector& InnVrtSel, VertexSelector& SurfVrtSel);

static	bool CollectSurfaceVertices(Grid& grid, SubsetHandler& shFace, VertexSelector& SurfVrtSel);

static	bool CollectAllVerticesForNG(Grid& grid, VertexSelector& NgVrtSel,
							VertexSelector& SurfVrtSel, VertexSelector& InnVrtSel);

static	bool WriteLGM(Grid& grid,
					  const char* lgmFilename,
					  const char* problemName,
					  const char* lgmName,
					  int convex,
					  SubsetHandler& sh,
					  EdgeSelector& LineSel,
					  VertexSelector& SurfVrtSel,
					  Grid::EdgeAttachmentAccessor<AInt>& aaLineIndex,
					  Grid::VertexAttachmentAccessor<AInt>& aaSurfVrtIndex,
					  Grid::VertexAttachmentAccessor<AVector3>& aaPos);

static	bool GetRightLeftUnitIndex(int& rightIndex, int& leftIndex, Grid& grid, Face* face,
								   SubsetHandler& shVolume);

static bool WriteNG(Grid& grid,
					 SubsetHandler& sh,
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
bool ExportGridToUG(Grid& grid, SubsetHandler& sh, const char* fileNamePrefix,
					const char* lgmName, const char* problemName, int convex)
{
//	initialization
	EdgeSelector	LineSel(grid);
	VertexSelector 	NgVrtSel(grid);
	VertexSelector	SurfVrtSel(grid);
	VertexSelector 	InnVrtSel(grid);

//TODO: disable those options again, if they weren't enabled before.

	grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	grid.enable_options(EDGEOPT_STORE_ASSOCIATED_FACES);
//	for write ng
	grid.enable_options(FACEOPT_STORE_ASSOCIATED_EDGES);
	grid.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);


	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPosition);

//	selection of lines, surface-vertices, inner vertices
	CollectLines(grid, sh, LineSel);
	CollectSurfaceVertices(grid, sh, SurfVrtSel);
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
		for(EdgeBaseIterator EIter = LineSel.begin(); EIter != LineSel.end(); ++EIter, ++counter)
			aaLineIndex[*EIter] = counter;

	//	assign index to every face regarding EACH surface
		AInt aFaceIndex;
		grid.attach_to_faces(aFaceIndex);
		Grid::FaceAttachmentAccessor<AInt> aaFaceIndex(grid, aFaceIndex);
		for(uint i = 0; i < sh.num_subsets(); ++i)
		{
			counter = 0;
		//	assign triangle indices
			for(TriangleIterator iter = sh.begin<Triangle>(i); iter != sh.end<Triangle>(i); ++iter, ++counter)
				aaFaceIndex[*iter] = counter;

		//	assign quadrilaterals. Increase index by 2, since 2 triangles are written for each quad
			for(QuadrilateralIterator iter = sh.begin<Quadrilateral>(i); iter != sh.end<Quadrilateral>(i); ++iter, counter += 2)
				aaFaceIndex[*iter] = counter;
		}

	//	assign index to every grid-vertex in *.ng file order
		AInt aNgVrtIndex;
		grid.attach_to_vertices(aNgVrtIndex);
		Grid::VertexAttachmentAccessor<AInt> aaNgVrtIndex(grid, aNgVrtIndex);
		counter = 0;
		for(VertexBaseIterator VIter = NgVrtSel.begin(); VIter != NgVrtSel.end(); ++VIter, ++counter)
			aaNgVrtIndex[*VIter] = counter;

	// 	assign index to every inner-vertex
		AInt aInnVrtIndex;
		grid.attach_to_vertices(aInnVrtIndex);
		Grid::VertexAttachmentAccessor<AInt> aaInnVrtIndex(grid, aInnVrtIndex);
		counter = 0;
		for(VertexBaseIterator VIter = InnVrtSel.begin(); VIter != InnVrtSel.end(); ++VIter, ++counter)
			aaInnVrtIndex[*VIter] = counter;

	//	assign index to every surface-vertex
		AInt aSurfVrtIndex;
		grid.attach_to_vertices(aSurfVrtIndex);
		Grid::VertexAttachmentAccessor<AInt> aaSurfVrtIndex(grid, aSurfVrtIndex);
		counter = 0;
		for(VertexBaseIterator SVIter = SurfVrtSel.begin(); SVIter != SurfVrtSel.end(); ++SVIter, ++counter)
			aaSurfVrtIndex[*SVIter] = counter;
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//	set up unit names
	{
		int unitCount = 1;
		for(uint i = 0; i < sh.num_subsets(); ++i)
		{
			if(sh.subset_info(i).name.size() == 0 && sh.num_elements<Volume>(i) > 0)
			{
				stringstream ss;

			//	be cautious with unit indices: id = 0 is reserved for outer space ONLY by UG
				ss << "unit_" << unitCount;
				sh.subset_info(i).name = ss.str();
				++unitCount;
			}
		}
	}

	string lgmFilename(fileNamePrefix);
	lgmFilename.append(".lgm");
//	write *.lgm file
	WriteLGM(grid, lgmFilename.c_str(), problemName, lgmName,
			 convex, sh, LineSel, SurfVrtSel, aaLineIndex, aaSurfVrtIndex, aaPos);

	string ngFilename(fileNamePrefix);
	ngFilename.append(".ng");

//	write *.ng file
	WriteNG(grid, sh, ngFilename.c_str(), SurfVrtSel, InnVrtSel, NgVrtSel, LineSel,
			aaLineIndex, aaInnVrtIndex, aaNgVrtIndex, aaFaceIndex, aaPos);

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	CollectLines
/// collects lines in a selector
static bool CollectLines(Grid& grid, SubsetHandler& shFace, EdgeSelector& LineSel)
{
//	store associated faces in this vector
	vector<Face*> vFaces;

//	iterate through all edges in the grid and identify lines by comparing the subset-membership
//	of the associated faces
	for(EdgeBaseIterator EIter = grid.edges_begin(); EIter != grid.edges_end(); ++EIter)
	{
		CollectFaces(vFaces, grid, *EIter);

		if(vFaces.size() > 1)
		{
			uint i = 0;
			Face* f1 = NULL;
		//	find the first face that is assigned to a subset
			for(; i < vFaces.size(); ++i)
			{
				f1 = vFaces[i];
				if(shFace.get_subset_index(f1) != -1)
					break;
			}

		//	compare with others. only check the ones that are assigned to a subset.
			for(; i < vFaces.size(); ++i)
			{
				if(shFace.get_subset_index(vFaces[i]) == -1)
					continue;

				if(shFace.get_subset_index(f1) != shFace.get_subset_index(vFaces[i]))
				{
					LineSel.select(*EIter);
					break;
				}
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
	for(VertexBaseIterator VIter = grid.vertices_begin(); VIter != grid.vertices_end(); ++VIter)
	{
		InnVrtSel.select(*VIter);
	}

//	iterate through all surface-vertices and deselect them in VertexSelector InnVrtSel, so that there only remain
//	the inner vertices
	for(VertexBaseIterator VIter = SurfVrtSel.begin(); VIter != SurfVrtSel.end(); ++VIter)
	{
		InnVrtSel.deselect(*VIter);
	}

	return true;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	CollectSurfaceVertices
/// collects surface-vertices in a selector
static bool CollectSurfaceVertices(Grid& grid, SubsetHandler& shFace, VertexSelector& SurfVrtSel)
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
	for(VertexBaseIterator VIter = SurfVrtSel.begin(); VIter != SurfVrtSel.end(); ++VIter)
	{
		NgVrtSel.select(*VIter);
	}

//	iterate through all inner vertices and select them
	for(VertexBaseIterator VIter = InnVrtSel.begin(); VIter != InnVrtSel.end(); ++VIter)
	{
		NgVrtSel.select(*VIter);
	}

	return true;
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
			  SubsetHandler& sh,
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

//	write the header
	out << "#Domain-Info" << endl;
	out << "name = " << lgmName << endl;
	out << "problemname = " << problemName << endl;
	out << "convex = " << convex << endl << endl;

//	write the units
	out << "#Unit-Info" << endl;
	for(uint i = 0; i < sh.num_subsets(); ++i)
	{
		if(sh.num_elements<Volume>(i) > 0)
		{
		//	be cautious with the unit indices:	id = 0 is reserved for outer space ONLY by UG
			out << "unit " << i + 1 << " " << sh.subset_info(i).name << endl;
		}
	}
	out << endl;

//	write the lines
	out << "#Line-Info" << endl;
	for(EdgeBaseIterator iter = LineSel.begin(); iter != LineSel.end(); ++iter)
	{
		EdgeBase* e = *iter;
		out << "line " << aaLineIndex[e] << ": points: ";
		out << aaSurfVrtIndex[e->vertex(0)] << " " << aaSurfVrtIndex[e->vertex(1)] << ";" << endl;
	}
	out << endl;

//	write the surfaces
	out << "#Surface-Info" << endl;
	{
	//	used to find vertices that belong to the i-th subset.
		VertexSelector tmpSurfVrtSel(grid);
		for(uint i = 0; i < sh.num_subsets(); ++i)
		{
			tmpSurfVrtSel.clear_selection();

			//	identify left and right unit index of each surface
			int tmpLeft, tmpRight;
			if(!GetRightLeftUnitIndex(tmpRight, tmpLeft, grid, *sh.begin<Face>(i), sh))
			{
				LOG("GetRightLeftUnitIndex failed during lgm-write.\n");
				LOG("This can happen due to elements with bad orinentation.\n");
				LOG("IMPLEMENT a geometrical method for fallback!\n");
			}

			out << "surface " << i << ": left=" << tmpLeft << "; right=" << tmpRight << "; points:";
			for(FaceIterator FIter = sh.begin<Face>(i); FIter != sh.end<Face>(i); ++FIter)
			{
				Face* f = *FIter;
				for(uint j = 0; j < f->num_vertices(); ++j)
				{
					VertexBase* v = f->vertex(j);
					if(!tmpSurfVrtSel.is_selected(v))
					{
						tmpSurfVrtSel.select(v);
						out << " " << aaSurfVrtIndex[v];
					}
				}
			}

		//	write lines
			out << "; lines:";
			{
				vector<EdgeBase*> vEdges;
				for(FaceIterator FIter = sh.begin<Face>(i); FIter != sh.end<Face>(i); ++FIter)
				{
					CollectEdges(vEdges, grid, *FIter);
					for(vector<EdgeBase*>::iterator EIter = vEdges.begin(); EIter != vEdges.end(); ++EIter)
					{
						EdgeBase* e = *EIter;
						if(LineSel.is_selected(e))
						{
							out << " " << aaLineIndex[e];
						}
					}
				}
			}

		// 	write triangles
			out << "; triangles:";
			{
				vector<VertexBase*> vVertices;
				for(TriangleIterator TIter = sh.begin<Triangle>(i); TIter != sh.end<Triangle>(i); ++TIter)
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
			{
				for(QuadrilateralIterator iter = sh.begin<Quadrilateral>(i); iter != sh.end<Quadrilateral>(i); ++iter)
				{
					Quadrilateral* q = *iter;
					for(uint j = 0; j < 3; ++j)
						out << " " << aaSurfVrtIndex[q->vertex(j)];

					out << ";";

					for(uint j = 2; j < 5; ++j)
						out << " " << aaSurfVrtIndex[q->vertex(j%4)];

					out << ";";
				}
			}
		//	surface written...
			out << endl;
		}
	}

//	write the points position data
	out << endl << "#Point-Info" << endl;
	for(VertexBaseIterator iter = SurfVrtSel.begin(); iter != SurfVrtSel.end(); ++iter)
	{
		out << aaPos[*iter].x << " " << aaPos[*iter].y << " " << aaPos[*iter].z << ";" << endl;
	}

//	close out file
	out.close();

	return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	GetRightLeftUnitIndices
/// identifies the right and left unit index to a given surface
static bool GetRightLeftUnitIndex(int& rightIndex, int& leftIndex, Grid& grid, Face* face, SubsetHandler& shVolume)
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
			v->face(i, fd);
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
					 SubsetHandler& sh,
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
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<bool> faceFlags(sh.num_subsets());

	ofstream out(ngFilename);
	if(!out)
		return false;
	out.setf(ios::scientific);

//	write the boundary nodes section
	out << "# boundary nodes" << endl;
	for(VertexBaseIterator VIter = SurfVrtSel.begin(); VIter != SurfVrtSel.end(); ++VIter)
	{
		VertexBase* v = *VIter;

		for(uint i = 0; i < sh.num_subsets(); ++i)
		{
				faceFlags[i] = false;
		}

	//	write the boundary-vertex position data
		out << "B " << aaPos[v].x << " " << aaPos[v].y << " " << aaPos[v].z << endl;
		out << "	";

	//	iterate through all lines
		for(EdgeBaseIterator EIter = LineSel.begin(); EIter != LineSel.end(); ++EIter)
		{
			EdgeBase* e = *EIter;

		//	check if the current boundary-vertex is also a line-vertex
			if(GetVertexIndex(e, v) != -1)
			{
			//	write the line-index and line-vertex-index (corresponds to local position)
				out << " L " << aaLineIndex[e] << " " << GetVertexIndex(e, v);
			}
		}

	//	surface data
	//	iterate through all associated boundary-faces of vertex v
		for(FaceIterator FIter = grid.associated_faces_begin(v); FIter != grid.associated_faces_end(v); ++FIter)
		{
			Face* f = *FIter;

		//	write the surface section including its index and the triangle index of
		//	ONE triangle-representative (per face-subset) associated to vertex v
			if(sh.get_subset_index(f) != -1)
			{
				if(faceFlags[sh.get_subset_index(f)] == false)
				{
				//	write local coordinates of vertex v on the triangle-representative
					int faceInd = aaFaceIndex[f];
					int vrtInd = GetVertexIndex(f, v);
					vector2 vCoord(0, 0);
					switch(vrtInd)
					{
						case 0:	vCoord = vector2(0, 0); break;
						case 1:	vCoord = vector2(1.0, 0); break;
						case 2:	vCoord = vector2(0, 1.0); break;
						case 3:	vCoord = vector2(1.0, 0); faceInd++; break; //index into the second sub-triangle
					}

					out << " S " << sh.get_subset_index(f) << " " << faceInd;
					out << " " << vCoord.x << " " << vCoord.y;

					faceFlags[sh.get_subset_index(f)] = true;
				}
			}
		}

	out << ";" << endl;
	}

//	write inner nodes
	out << endl;
	out << "# inner nodes" << endl;

	if(InnVrtSel.num_selected() > 0)
	{
	//	if there are inner vertices, iterate through all of them
		for(VertexBaseIterator VIter = InnVrtSel.begin(); VIter != InnVrtSel.end(); ++VIter)
		{
			VertexBase* v = *VIter;

		//	write the position data of the inner vertices
			out << "I " << aaPos[v].x << " " << aaPos[v].y << " " << aaPos[v].z << ";" << endl;
		}
	}

//	write elements
	out << endl;
	out << "# elements" << endl;

	for(uint i = 0; i < sh.num_subsets(); ++i)
	{
	// 	iterate through all volumes
		for(VolumeIterator VIter = sh.begin<Volume>(i); VIter != sh.end<Volume>(i); ++VIter)
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
				if(sh.get_subset_index(f) != -1)
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

}//	end of namespace
