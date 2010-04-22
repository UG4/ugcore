// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m05 d28

#include "extrude.h"
#include "common/util/hash.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"

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
			APosition& aPos,
			bool invertOrientation)
{
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
	typedef Hash<VertexBase*, uint> VertexHash;
	VertexHash vrtHash(hashSize);

//	we'll record created faces in this vector, since we have to fix the
//	orientation later on (only if pvEdgesInOut has been specified).
	vector<Face*> vNewFaces;
	bool bRecordNewFaces = (pvEdgesInOut != NULL) && (extrusionOptions & EO_CREATE_FACES);
	if(bRecordNewFaces)
		vNewFaces.resize(numNewFaces);
	size_t newFaceCount = 0;
//	if faces are extruded we want them in the first section of vNewFaces (they shall define the orientation).
//	Thats why we start in the second section in this case.
	if(pvFacesInOut)
		newFaceCount = pvFacesInOut->size();

//	first we'll extrude all vertices. For each a new edge will be created.
	if(pvVerticesInOut)
	{
		vector<VertexBase*>& vVertices = *pvVerticesInOut;
		for(uint i = 0; i < vVertices.size(); ++i)
		{
			VertexBase* vOld = vVertices[i];
		//	create a new vertex and store it in the hash.
		//	use the attachment_data_index of the old one as key.
		//	WARNING: this is only secure as long as nobody calls defragment while the hash is active!
			VertexBase* v = *grid.create<Vertex>(vOld);
			vrtHash.add(v, grid.get_attachment_data_index(vOld));

		//	calculate new position
			aaPos[v] = aaPos[vOld];
			VecAdd(aaPos[v], aaPos[v], direction);

		//	create an edge between both vertices.
			grid.create<Edge>(EdgeDescriptor(vOld, v), vOld);

		//	overwrite the vertex in pvVerticesInOut
			vVertices[i] = v;
		}
	}

//	now extrude edges.
	if(pvEdgesInOut)
	{
		vector<EdgeBase*>& vEdges = *pvEdgesInOut;
		for(uint i = 0; i < vEdges.size(); ++i)
		{
			EdgeBase* e = vEdges[i];

		//	check for both boundary points whether the new vertices have already been created.
		//	if not then create them here and store them.
			VertexBase* v[2];

			for(uint j = 0; j < 2; ++j)
			{
				if(vrtHash.has_entries(grid.get_attachment_data_index(e->vertex(j))))
					v[j] = vrtHash.first(grid.get_attachment_data_index(e->vertex(j)));
				else
				{
					v[j] = *grid.create<Vertex>(e->vertex(j));
					vrtHash.add(v[j], grid.get_attachment_data_index(e->vertex(j)));
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
	}

//	now extrude faces.
	if(pvFacesInOut)
	{
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
				if(vrtHash.has_entries(oldVrtInd))
				{
					v[j] = vrtHash.first(oldVrtInd);
				}
				else
				{
					v[j] = *grid.create<Vertex>(f->vertex(j));
					vrtHash.add(v[j], oldVrtInd);
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
				fNew = *grid.create<Triangle>(TriangleDescriptor(v[0], v[1], v[2]), f);

			//	create the volume
				if(extrusionOptions & EO_CREATE_VOLUMES)
				{
				//	create a hexahedron from the four points.
					if(invertOrientation)
						grid.create<Prism>(PrismDescriptor(f->vertex(2), f->vertex(1), f->vertex(0),
														fNew->vertex(2), fNew->vertex(1), fNew->vertex(0)), f);
					else
						grid.create<Prism>(PrismDescriptor(f->vertex(0), f->vertex(1), f->vertex(2),
														fNew->vertex(0), fNew->vertex(1), fNew->vertex(2)), f);
				}
			}
			else if(numVrts == 4)
			{
				fNew = *grid.create<Quadrilateral>(QuadrilateralDescriptor(v[0], v[1], v[2], v[3]), f);

			//	create the volume
				if(extrusionOptions & EO_CREATE_VOLUMES)
				{
				//	create a hexahedron from the four points.
					if(invertOrientation)
						grid.create<Hexahedron>(HexahedronDescriptor(f->vertex(3), f->vertex(2), f->vertex(1), f->vertex(0),
														fNew->vertex(3), fNew->vertex(2), fNew->vertex(1), fNew->vertex(0)), f);
					else
						grid.create<Hexahedron>(HexahedronDescriptor(f->vertex(0), f->vertex(1), f->vertex(2), f->vertex(3),
														fNew->vertex(0), fNew->vertex(1), fNew->vertex(2), fNew->vertex(3)), f);
				}
			}
			else
			{
				assert(!"face has bad number of vertices!");
			}

			if(fNew)
			{
			//	overwrite the face in pvFacesInOut
				vFaces[i] = fNew;
			//	store the face in vNewFaces - but only if we have
			//	to re-orientate them later on.
				if(bRecordNewFaces)
					vNewFaces[newFaceCount++] = fNew;
			}
		}
	}

//	if faces were extruded from edges, we have to fix the orientation now
	if(bRecordNewFaces){
		FixOrientation(grid, vNewFaces.begin(), vNewFaces.end());
	}
}

}//	end of namespace
