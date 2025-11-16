/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include <fstream>
#include <vector>
#include <cstring>
#include <string>
#include "file_io_obj.h"
#include "../lg_base.h"
#include "common/util/loader/loader_obj.h"
#include "common/error.h"

using namespace std;

namespace ug
{
//TODO: disable automatic triangulation. create edges. create quads.
////////////////////////////////////////////////////////////////////////
//	LoadGridFromOBJ
///	Loads a file from a wavefront '.obj' file. Fills optional subset-infos.
bool LoadGridFromOBJ(Grid& grid, const char* filename, AVector3& aPos,
		AVector2* paTexCoord,
		ISubsetHandler* pSubsetHandler,
		std::vector<OBJMaterial>* pvMaterials)
{
	LoaderObj loader;

	if(!loader.load_file(filename, false))
		return false;

	if(loader.num_points() == 0)
		return false;

//	attach data and attachment accessors
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<AVector3> aaPosVRT(grid, aPos);

//	create vertices and assign position data.
//	Store pointers in a vector so that they can be accessed by index.
	vector<RegularVertex*> vVertices;
	vVertices.reserve(loader.num_points());
	{
		for(uint i = 0; i < (uint)loader.num_points(); ++i)
		{
			RegularVertex* pVrt = *grid.create<RegularVertex>();
			vVertices.push_back(pVrt);
			aaPosVRT[pVrt] = *loader.point(i);
		}
	//	if a subset handler was specified, we'll push all vertices into subset 0
		if(pSubsetHandler){
			pSubsetHandler->assign_subset(grid.begin<RegularVertex>(), grid.end<RegularVertex>(), 0);
		}
	}

//	iterate through the objects in loader and add edges and faces to grid.
//	assign subsets and material indices at the same time.
	{
		int objCounter = 0;
		for(auto oIter = loader.objects_begin(); oIter != loader.objects_end(); ++oIter, ++objCounter)
		{
			LoaderObj::Object& obj = *(*oIter);
		//	set the subset-info of this object - if requested.
			if(pSubsetHandler != nullptr)
			{
				SubsetInfo si = pSubsetHandler->subset_info(objCounter);
				si.materialIndex = obj.m_iMaterialIndex;
				si.name = obj.m_strName;

				pSubsetHandler->set_subset_info(objCounter, si);
			}

		//	create edges and faces.
		//	creation differs if a subset handler is supplied.
			if(pSubsetHandler == nullptr)
			{
			//	create edges
				for(size_t i = 0; i < obj.m_vEdgeList.size(); i+=2){
					for(size_t j = i; j < i+2; ++j){
						if(obj.m_vEdgeList[j] < 0 || obj.m_vEdgeList[j] >= (int)vVertices.size()){
							UG_THROW("Bad vertex index in edge " << i / 2 << ": " << obj.m_vEdgeList[i] + 1);
						}
					}
					grid.create<RegularEdge>(EdgeDescriptor(vVertices[obj.m_vEdgeList[i]],
													vVertices[obj.m_vEdgeList[i+1]]));
				}
			//	create triangles
				for(size_t i = 0; i < obj.m_vTriangleList.size(); i+=3){
					for(size_t j = i; j < i+3; ++j){
						if(obj.m_vTriangleList[j] < 0 || obj.m_vTriangleList[j] >= (int)vVertices.size()){
							UG_THROW("Bad vertex index in triangle " << i / 3 << ": " << obj.m_vTriangleList[i] + 1);
						}
					}
					grid.create<Triangle>(TriangleDescriptor(vVertices[obj.m_vTriangleList[i]],
										vVertices[obj.m_vTriangleList[i+1]], vVertices[obj.m_vTriangleList[i+2]]));
				}

			//	create quads
				for(size_t i = 0; i < obj.m_vQuadList.size(); i+=4){
					for(size_t j = i; j < i+4; ++j){
						if(obj.m_vQuadList[j] < 0 || obj.m_vQuadList[j] >= (int)vVertices.size()){
							UG_THROW("Bad vertex index in quadrilateral " << i / 4 << ": " << obj.m_vQuadList[i] + 1);
						}
					}
					grid.create<Quadrilateral>(QuadrilateralDescriptor(vVertices[obj.m_vQuadList[i]],
										vVertices[obj.m_vQuadList[i+1]], vVertices[obj.m_vQuadList[i+2]],
										vVertices[obj.m_vQuadList[i+3]]));
				}
			}
			else
			{
			//	create edges
				for(size_t i = 0; i < obj.m_vEdgeList.size(); i+=2)
				{
					for(size_t j = i; j < i+2; ++j){
						if(obj.m_vEdgeList[j] < 0 || obj.m_vEdgeList[j] >= (int)vVertices.size()){
							UG_THROW("Bad vertex index in edge " << i / 2 << ": " << obj.m_vEdgeList[i] + 1);
						}
					}
					RegularEdge* e = *grid.create<RegularEdge>(EdgeDescriptor(vVertices[obj.m_vEdgeList[i]],
													vVertices[obj.m_vEdgeList[i+1]]));
					pSubsetHandler->assign_subset(e, objCounter);
				}

			//	create triangles
				for(size_t i = 0; i < obj.m_vTriangleList.size(); i+=3)
				{
					for(size_t j = i; j < i+3; ++j){
						if(obj.m_vTriangleList[j] < 0 || obj.m_vTriangleList[j] >= (int)vVertices.size()){
							UG_THROW("Bad vertex index in triangle " << i / 3 << ": " << obj.m_vTriangleList[i] + 1);
						}
					}
					Triangle* tri = *grid.create<Triangle>(TriangleDescriptor(vVertices[obj.m_vTriangleList[i]],
										vVertices[obj.m_vTriangleList[i+1]], vVertices[obj.m_vTriangleList[i+2]]));
					pSubsetHandler->assign_subset(tri, objCounter);
				}

			//	create quads
				for(size_t i = 0; i < obj.m_vQuadList.size(); i+=4)
				{
					for(size_t j = i; j < i+4; ++j){
						if(obj.m_vQuadList[j] < 0 || obj.m_vQuadList[j] >= (int)vVertices.size()){
							UG_THROW("Bad vertex index in quadrilateral " << i / 4 << ": " << obj.m_vQuadList[i] + 1);
						}
					}
					Quadrilateral* q = *grid.create<Quadrilateral>(QuadrilateralDescriptor(vVertices[obj.m_vQuadList[i]],
										vVertices[obj.m_vQuadList[i+1]], vVertices[obj.m_vQuadList[i+2]],
										vVertices[obj.m_vQuadList[i+3]]));
					pSubsetHandler->assign_subset(q, objCounter);
				}
			}

		//	assign texture coords if the corresponding attachment is supplied
			if(paTexCoord != nullptr)
			{
			//	since you can specify texture coords for each face from the .obj specification,
			//	and since libGrid only supports one texture coord per vertex, we reduce
			//	the supplied information.
				Grid::VertexAttachmentAccessor<AVector2> aaTexVRT(grid, *paTexCoord);
				for(int i= 0; i < (int)obj.m_vTriangleListTex.size(); ++i)
					aaTexVRT[vVertices[i]] = *loader.uv(obj.m_vTriangleListTex[i]);
			}
		}
	}

//	copy materials to pvMaterials - if supplied
	if(pvMaterials != nullptr)
	{
		pvMaterials->resize(loader.num_materials());
		for(int i = 0; i < loader.num_materials(); ++i)
		{
			const LoaderObj::Material& mat = loader.get_material(i);
			(*pvMaterials)[i].m_strName = mat.m_strName;
			(*pvMaterials)[i].m_strTextureDiffuse = mat.m_strTextureDiffuse;
			(*pvMaterials)[i].m_vDiffuse = mat.m_vDiffuse;
			(*pvMaterials)[i].m_fAlpha = mat.m_fAlpha;
		}
	}

	return true;
}


static void WriteEdges(ofstream& out, EdgeIterator edgesBegin, EdgeIterator edgesEnd,
						int indexDimension, Grid::VertexAttachmentAccessor<AInt>& aaInt)
{
	while(edgesBegin != edgesEnd)
	{
		Edge* e = *edgesBegin;
		out << "f";
		for(uint i = 0; i < 2; ++i)
		{
			out << " " << aaInt[e->vertex(i)];
			for(int j = 1; j < indexDimension; ++j)
				out << '/' << aaInt[e->vertex(i)];
		}
		++edgesBegin;
		out << endl;
	}
}

static void WriteFaces(ofstream& out, FaceIterator facesBegin, FaceIterator facesEnd,
						int indexDimension, Grid::VertexAttachmentAccessor<AInt>& aaInt)
{
	while(facesBegin != facesEnd)
	{
		Face* f = *facesBegin;
		out << "f";
		for(uint i = 0; i < f->num_vertices(); ++i)
		{
			out << " " << aaInt[f->vertex(i)];
			for(int j = 1; j < indexDimension; ++j)
				out << '/' << aaInt[f->vertex(i)];
		}
		++facesBegin;
		out << endl;
	}
}

bool SaveGridToOBJ(Grid& grid, const char* filename, AVector3& aPos,
		AVector2* paTexCoord,
		ISubsetHandler* pSubsetHandler,
		std::vector<OBJMaterial>* pvMaterials)
{
//	check whether the grid has vertex-position data. If not, return false
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("ERROR in SaveGridToOBJ: Position attachment missing.\n");
		return false;
	}

	string mtrlFullFilename;
	string mtrlFilename;

//	check if there are materials
	if(pvMaterials)
	{
	//	we have to construct the material-file name.
		mtrlFullFilename = filename;
		if(mtrlFullFilename.at(mtrlFullFilename.size() - 4) == '.')
			mtrlFullFilename.erase(mtrlFullFilename.size() - 4, 4);
		mtrlFullFilename.append(".mtl");

	//	extract the name of the file and store it in mtrlFilename
		string::size_type ind1, ind2;
		ind1 = mtrlFullFilename.find_last_of('\\');
		ind2 = mtrlFullFilename.find_last_of('/');
		if(ind1 == string::npos)
			ind1 = ind2;
		else if(ind2 != string::npos)
			ind1 = max(ind1, ind2);

	//	assign the filename
		if(ind1 == string::npos)
			mtrlFilename = mtrlFullFilename;
		else
			mtrlFilename.assign(mtrlFullFilename, ind1, mtrlFullFilename.size() - ind1);
	}

//	write the .obj file
	{
		ofstream out(filename);

	//	write the header
		out << "# exported from libGrid" << endl;

	//	write material file
		if(pvMaterials)
			out << "mtllib " << mtrlFilename << endl;

		int indexDimension = 1;

	//	store indices in the vertices
		AInt aInt;
		grid.attach_to_vertices(aInt);
		Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);

	//	write vertex data
		{
			int counter = 1;
			Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
			for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter, ++counter)
			{
				out << "v " << aaPos[*iter].x() << " " << aaPos[*iter].y() << " " <<
								aaPos[*iter].z() << endl;
				aaInt[*iter] = counter;
			}
		}

	//	if texture data is supplied we have to write it too
		if(paTexCoord != nullptr)
		{
			indexDimension++;
			Grid::VertexAttachmentAccessor<AVector2> aaTex(grid, *paTexCoord);
			for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter)
			{
				out << "vt " << aaTex[*iter].x() << " " <<
								aaTex[*iter].y() << " 0.0" << endl;
			}
		}

	//	if a subset handler is supplied, we'll export each subset as a separate object.
		if(pSubsetHandler)
		{
			for(int i = 0; i < pSubsetHandler->num_subsets(); ++i)
			{
			//	write object
				if(pSubsetHandler->subset_info(i).name.size() > 0)
					out << "o " << pSubsetHandler->subset_info(i).name << endl;
				else
					out << "o " << "sub_" << i << endl;

			//	write material reference
				if((pvMaterials != nullptr) && (pSubsetHandler->subset_info(i).materialIndex != -1))
					out << "usemtl " << (*pvMaterials)[pSubsetHandler->subset_info(i).materialIndex].m_strName << endl;
				else
					out << "usemtl (null)" << endl;

			//	get the goc
				GridObjectCollection goc = pSubsetHandler->
						get_grid_objects_in_subset(i);

				for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
				//	write the edges of this subset
					WriteEdges(out, goc.begin<Edge>(lvl),
								goc.end<Edge>(lvl), indexDimension, aaInt);

				//	write the faces of this subset
					WriteFaces(out, goc.begin<Face>(lvl), goc.end<Face>(lvl),
							   indexDimension, aaInt);
				}
			}
		}
		else
		{
		//	since no subset handler was supplied, we'll treat the whole grid as on object.
			out << "o nonameknown" << endl;
			out << "usemtl (null)" << endl;
		//	write the edges
			WriteEdges(out, grid.edges_begin(), grid.edges_end(), indexDimension, aaInt);

		//	write the faces
			WriteFaces(out, grid.faces_begin(), grid.faces_end(), indexDimension, aaInt);
		}

	//	clean up
		grid.detach_from_vertices(aInt);
	//	.obj done
		out.close();
	}

//	write the material file
	if(pvMaterials)
	{
		ofstream out(mtrlFullFilename.c_str());
		out << "# exported from libMesh" << endl;
		out << "# Material Count: " << pvMaterials->size() << endl;

		for(uint i = 0; i < pvMaterials->size(); ++i)
		{
			vector4& vD = (*pvMaterials)[i].m_vDiffuse;
			out << "newmtl " << (*pvMaterials)[i].m_strName << endl;
			out << "Kd " << vD.x() << " " << vD.y() << " " << vD.z() << endl;
			out << "d " << vD.w() << endl;
			if((*pvMaterials)[i].m_strTextureDiffuse.size() > 0)
				out << "map_Kd " << (*pvMaterials)[i].m_strTextureDiffuse << endl;
		}

	//	done
		out.close();
	}

	return true;
}

}//	end of namespace
