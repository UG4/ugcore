// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m05 d28

#include "extrude.h"
#include "common/util/hash.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_grid/algorithms/geom_obj_util/volume_util.h"

using namespace std;

namespace ug
{

static bool ExtrusionHelper_CheckOrientation(Volume* v, Grid::VertexAttachmentAccessor<Attachment<vector1> >& aaPos)
{
    UG_THROW("Can't check orientation of a degenerated volume element!");
    return false;
}

static bool ExtrusionHelper_CheckOrientation(Volume* v, Grid::VertexAttachmentAccessor<Attachment<vector2> >& aaPos)
{
    UG_THROW("Can't check orientation of a degenerated volume element!");
    return false;
}

static bool ExtrusionHelper_CheckOrientation(Volume* v, Grid::VertexAttachmentAccessor<Attachment<vector3> >& aaPos)
{
    return CheckOrientation(v, aaPos);
}

////////////////////////////////////////////////////////////////////////
//	Extrude
template <class vector_t>
void Extrude(Grid& grid,
			std::vector<Vertex*>* pvVerticesInOut,
			std::vector<Edge*>* pvEdgesInOut,
			std::vector<Face*>* pvFacesInOut,
			const vector_t& direction,
			uint extrusionOptions,
			Attachment<vector_t>& aPos,
			std::vector<Volume*>* pvVolsOut)
{
	UG_DLOG(LIB_GRID, 0, "extruding...\n");
	
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	Grid::VertexAttachmentAccessor<Attachment<vector_t> > aaPos(grid, aPos);

	Extrude(grid, pvVerticesInOut, pvEdgesInOut, pvFacesInOut,
			direction, aaPos, extrusionOptions, pvVolsOut);
}

template <class TAAPos>
void Extrude(Grid& grid,
			std::vector<Vertex*>* pvVerticesInOut,
			std::vector<Edge*>* pvEdgesInOut,
			std::vector<Face*>* pvFacesInOut,
			const typename TAAPos::ValueType& direction,
			TAAPos aaPos,
			uint extrusionOptions,
			std::vector<Volume*>* pvVolsOut)
{
	if(pvVolsOut)
		pvVolsOut->clear();

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
	typedef Hash<uint, Vertex*> VertexHash;
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
		vector<Vertex*>& vVertices = *pvVerticesInOut;
		for(uint i = 0; i < vVertices.size(); ++i)
		{
			Vertex* vOld = vVertices[i];
		//	create a new vertex and store it in the hash.
		//	use the attachment_data_index of the old one as key.
		//	WARNING: this is only secure as long as nobody calls defragment while the hash is active!
			Vertex* v = *grid.create<RegularVertex>(vOld);
			vrtHash.insert(grid.get_attachment_data_index(vOld), v);

		//	calculate new position
			aaPos[v] = aaPos[vOld];
			VecAdd(aaPos[v], aaPos[v], direction);

		//	create an edge between both vertices.
			grid.create<RegularEdge>(EdgeDescriptor(vOld, v), vOld);

		//	overwrite the vertex in pvVerticesInOut
			vVertices[i] = v;
		}
		UG_DLOG(LIB_GRID, 1, "  extruding vertices done.\n");
	}

//	now extrude edges.
	if(pvEdgesInOut)
	{
		UG_DLOG(LIB_GRID, 1, "  extruding edges: " << pvEdgesInOut->size() << endl);
		vector<Edge*>& vEdges = *pvEdgesInOut;
		for(uint i = 0; i < vEdges.size(); ++i)
		{
			Edge* e = vEdges[i];

		//	check for both boundary points whether the new vertices have already been created.
		//	if not then create them here and store them.
			Vertex* v[2];

			for(uint j = 0; j < 2; ++j)
			{
				if(!vrtHash.get_entry(v[j], grid.get_attachment_data_index(e->vertex(j))))
				{
					v[j] = *grid.create<RegularVertex>(e->vertex(j));
					vrtHash.insert(grid.get_attachment_data_index(e->vertex(j)), v[j]);
				//	calculate new position
					aaPos[v[j]] = aaPos[e->vertex(j)];
					VecAdd(aaPos[v[j]], aaPos[v[j]], direction);
				}
			}

		//	both new vertices exist now.
		//	create the new edge
			RegularEdge* eNew = *grid.create<RegularEdge>(EdgeDescriptor(v[0], v[1]), e);

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
			Vertex* v[4];

			uint numVrts = f->num_vertices();

			for(uint j = 0; j < numVrts; ++j)
			{
				uint oldVrtInd = grid.get_attachment_data_index(f->vertex(j));
				if(!vrtHash.get_entry(v[j], oldVrtInd)){
					v[j] = *grid.create<RegularVertex>(f->vertex(j));
					vrtHash.insert(oldVrtInd, v[j]);
				//	calculate new position
					aaPos[v[j]] = aaPos[f->vertex(j)];
					VecAdd(aaPos[v[j]], aaPos[v[j]], direction);
				}
			}

		//	all new vertices exist now.
		//	create the new face
			Face* fNew = NULL;
			Volume* vol = NULL;
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
					if (!ExtrusionHelper_CheckOrientation(&prism, aaPos)){
						VolumeDescriptor vd;
						prism.get_flipped_orientation(vd);
						vol = *grid.create_by_cloning(&prism, vd, f);	
					}
					else
						vol = *grid.create_by_cloning(&prism, prism, f);					

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
					if(!ExtrusionHelper_CheckOrientation(&hex, aaPos)){
						VolumeDescriptor vd;
						hex.get_flipped_orientation(vd);
						vol = *grid.create_by_cloning(&hex, vd, f);	
					}
					else
						vol = *grid.create_by_cloning(&hex, hex, f);					
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

			if(pvVolsOut && vol)
				pvVolsOut->push_back(vol);
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

template void Extrude<vector1>(Grid&,
								std::vector<Vertex*>*,
								std::vector<Edge*>*,
								std::vector<Face*>*,
								const vector1&,
								uint,
								Attachment<vector1>&,
								std::vector<Volume*>*);

template void Extrude<vector2>(Grid&,
								std::vector<Vertex*>*,
								std::vector<Edge*>*,
								std::vector<Face*>*,
								const vector2&,
								uint,
								Attachment<vector2>&,
								std::vector<Volume*>*);

template void Extrude<vector3>(Grid&,
								std::vector<Vertex*>*,
								std::vector<Edge*>*,
								std::vector<Face*>*,
								const vector3&,
								uint,
								Attachment<vector3>&,
								std::vector<Volume*>*);

template void Extrude<Grid::VertexAttachmentAccessor<Attachment<vector1> > >(
								Grid&,
								std::vector<Vertex*>*,
								std::vector<Edge*>*,
								std::vector<Face*>*,
								const vector1&,
								Grid::VertexAttachmentAccessor<Attachment<vector1> >,
								uint,
								std::vector<Volume*>*);

template void Extrude<Grid::VertexAttachmentAccessor<Attachment<vector2> > >(
								Grid&,
								std::vector<Vertex*>*,
								std::vector<Edge*>*,
								std::vector<Face*>*,
								const vector2&,
								Grid::VertexAttachmentAccessor<Attachment<vector2> >,
								uint,
								std::vector<Volume*>*);

template void Extrude<Grid::VertexAttachmentAccessor<Attachment<vector3> > >(
								Grid&,
								std::vector<Vertex*>*,
								std::vector<Edge*>*,
								std::vector<Face*>*,
								const vector3&,
								Grid::VertexAttachmentAccessor<Attachment<vector3> >,
								uint,
								std::vector<Volume*>*);
}//	end of namespace
