// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m05 d28

#include "extrude.h"
#include "common/util/new_hash.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_grid/algorithms/geom_obj_util/volume_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	Extrude
void Extrude(Grid& grid,
			std::vector<VertexBase*>* pvVerticesInOut,
			std::vector<EdgeBase*>* pvEdgesInOut,
			std::vector<Face*>* pvFacesInOut,
			const vector3& direction,
			uint extrusionOptions,
			APosition& aPos)
{
	UG_DLOG(LIB_GRID, 0, "extruding...\n");
	
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//TODO: find a better guess.
//	first we'll determine the rough size of the hash.
//	we'll simply add all vertices, edges and faces.
	uint hashSize = 10;
	uint numNewFaces = 0;//	the number of faces that are generated directly in this method (autogeneration during volume-extrude is not regarded).
	if(pvVerticesInOut)
		hashSize += pvVerticesInOut->size();
	if(pvEdgesInOut){
		hashSize += pvEdgesInOut->size();
		numNewFaces += pvEdgesInOut->size();
	}
	if(pvFacesInOut){
		hashSize += pvFacesInOut->size();
		numNewFaces += pvFacesInOut->size();
	}

	if(hashSize == 0)
		return;

//	the hash:
	typedef NewHash<uint, VertexBase*> VertexHash;
	VertexHash vrtHash(hashSize);
	vrtHash.reserve(hashSize);

//	we'll record created faces in this vector, since we have to fix the
//	orientation later on (only if pvEdgesInOut has been specified).
	vector<Face*> vNewFaces;
	bool bRecordNewFaces = false;
	if((pvEdgesInOut != NULL) && (extrusionOptions & EO_CREATE_FACES))
	{
		if(!pvEdgesInOut->empty()){
			bRecordNewFaces = true;
			vNewFaces.resize(numNewFaces);
		}
	}
	
	size_t newFaceCount = 0;
//	if faces are extruded we want them in the first section of vNewFaces (they shall define the orientation).
//	Thats why we start in the second section in this case.
	if(pvFacesInOut)
		newFaceCount = pvFacesInOut->size();

//	first we'll extrude all vertices. For each a new edge will be created.
	if(pvVerticesInOut)
	{
		UG_DLOG(LIB_GRID, 1, "  extruding vertices: " << pvVerticesInOut->size() << endl);
		vector<VertexBase*>& vVertices = *pvVerticesInOut;
		for(uint i = 0; i < vVertices.size(); ++i)
		{
			VertexBase* vOld = vVertices[i];
		//	create a new vertex and store it in the hash.
		//	use the attachment_data_index of the old one as key.
		//	WARNING: this is only secure as long as nobody calls defragment while the hash is active!
			VertexBase* v = *grid.create<Vertex>(vOld);
			vrtHash.insert(grid.get_attachment_data_index(vOld), v);

		//	calculate new position
			aaPos[v] = aaPos[vOld];
			VecAdd(aaPos[v], aaPos[v], direction);

		//	create an edge between both vertices.
			grid.create<Edge>(EdgeDescriptor(vOld, v), vOld);

		//	overwrite the vertex in pvVerticesInOut
			vVertices[i] = v;
		}
		UG_DLOG(LIB_GRID, 1, "  extrunding vertices done.\n");
	}

//	now extrude edges.
	if(pvEdgesInOut)
	{
		UG_DLOG(LIB_GRID, 1, "  extruding edges: " << pvEdgesInOut->size() << endl);
		vector<EdgeBase*>& vEdges = *pvEdgesInOut;
		for(uint i = 0; i < vEdges.size(); ++i)
		{
			EdgeBase* e = vEdges[i];

		//	check for both boundary points whether the new vertices have already been created.
		//	if not then create them here and store them.
			VertexBase* v[2];

			for(uint j = 0; j < 2; ++j)
			{
				if(!vrtHash.get_entry(v[j], grid.get_attachment_data_index(e->vertex(j))))
				{
					v[j] = *grid.create<Vertex>(e->vertex(j));
					vrtHash.insert(grid.get_attachment_data_index(e->vertex(j)), v[j]);
				//	calculate new position
					aaPos[v[j]] = aaPos[e->vertex(j)];
					VecAdd(aaPos[v[j]], aaPos[v[j]], direction);
				}
			}

		//	both new vertices exist now.
		//	create the new edge
			Edge* eNew = *grid.create<Edge>(EdgeDescriptor(v[0], v[1]), e);

		//	overwrite the edge in pvEdgesInOut
			vEdges[i] = eNew;

		//	finally create the face.
			if(extrusionOptions & EO_CREATE_FACES)
			{
			//	create a quadrilateral from the four points.
			//	the orientation will be done later on.
				Face* f = *grid.create<Quadrilateral>(QuadrilateralDescriptor(v[0], v[1], e->vertex(1), e->vertex(0)), e);
				vNewFaces[newFaceCount++] = f;
			}
		}
		UG_DLOG(LIB_GRID, 1, "  extrunding edges done.\n");
	}

//	now extrude faces.
	if(pvFacesInOut)
	{
		UG_DLOG(LIB_GRID, 1, "  extruding faces: " << pvFacesInOut->size() << endl);
	//	fill the first section of vNewFaces
		newFaceCount = 0;

		vector<Face*>& vFaces = *pvFacesInOut;
		for(uint i = 0; i < vFaces.size(); ++i)
		{
			Face* f = vFaces[i];
			assert(f->num_vertices() < 5 && "can't deal with faces that have more than 4 vertices!");

		//	check for all points whether the new vertices have already been created.
		//	if not then create them here and store them.
			VertexBase* v[4];

			uint numVrts = f->num_vertices();

			for(uint j = 0; j < numVrts; ++j)
			{
				uint oldVrtInd = grid.get_attachment_data_index(f->vertex(j));
				if(!vrtHash.get_entry(v[j], oldVrtInd)){
					v[j] = *grid.create<Vertex>(f->vertex(j));
					vrtHash.insert(oldVrtInd, v[j]);
				//	calculate new position
					aaPos[v[j]] = aaPos[f->vertex(j)];
					VecAdd(aaPos[v[j]], aaPos[v[j]], direction);
				}
			}

		//	all new vertices exist now.
		//	create the new face
			Face* fNew = NULL;
			if(numVrts == 3)
			{
				UG_DLOG(LIB_GRID, 2, "    " << i << ": creating tri...\n");
				fNew = *grid.create<Triangle>(TriangleDescriptor(v[0], v[1], v[2]), f);

			//	create the volume
				if(extrusionOptions & EO_CREATE_VOLUMES)
				{
					Prism prism(f->vertex(0), f->vertex(1), f->vertex(2),
								fNew->vertex(0), fNew->vertex(1), fNew->vertex(2));
					
				//	check the orientation and create the prism
					if(!CheckOrientation(&prism, aaPos)){
						VolumeDescriptor vd;
						prism.get_flipped_orientation(vd);
						grid.create_by_cloning(&prism, vd, f);	
					}
					else
						grid.create_by_cloning(&prism, prism, f);					

				}
			}
			else if(numVrts == 4)
			{
				UG_DLOG(LIB_GRID, 2, "    " << i << ": creating quad...\n");
				fNew = *grid.create<Quadrilateral>(QuadrilateralDescriptor(v[0], v[1], v[2], v[3]), f);

			//	create the volume
				if(extrusionOptions & EO_CREATE_VOLUMES)
				{
					Hexahedron hex(f->vertex(0), f->vertex(1), f->vertex(2), f->vertex(3),
								   fNew->vertex(0), fNew->vertex(1), fNew->vertex(2), fNew->vertex(3));
					
				//	check the orientation and create the hexahedron
					if(!CheckOrientation(&hex, aaPos)){
						VolumeDescriptor vd;
						hex.get_flipped_orientation(vd);
						grid.create_by_cloning(&hex, vd, f);	
					}
					else
						grid.create_by_cloning(&hex, hex, f);					
				}
			}
			else{
				assert(!"face has bad number of vertices!");
			}

			if(fNew){
			//	overwrite the face in pvFacesInOut
				vFaces[i] = fNew;
			//	store the face in vNewFaces - but only if we have
			//	to re-orientate them later on.
				if(bRecordNewFaces){
					UG_DLOG(LIB_GRID, 2, "    storing face for reordering...\n");
					vNewFaces[newFaceCount++] = fNew;
				}
			}
		}
		UG_DLOG(LIB_GRID, 1, "  extruding faces done.\n");
	}

//	if faces were extruded from edges, we have to fix the orientation now
	if(bRecordNewFaces){
		UG_DLOG(LIB_GRID, 1, "  reordering faces...\n");
		FixFaceOrientation(grid, vNewFaces.begin(), vNewFaces.end());
		UG_DLOG(LIB_GRID, 1, "  reordering faces done.\n");
	}
}

}//	end of namespace
