//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d12

#include <queue>
#include "lib_grid/lg_base.h"
#include "../geom_obj_util/geom_obj_util.h"
#include "../refinement/regular_refinement.h"
#include "grid_adaption.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
					   VertexBase* vrtCenter, const vector3& normal,
					   number radius, APosition& aPos)
{
	AInt aInt;
	grid.attach_to_edges(aInt);
	bool retVal = AdaptSurfaceGridToCylinder(selOut, grid, vrtCenter,
											normal, radius, aInt, aPos);
	grid.detach_from_edges(aInt);
	return retVal;
}

////////////////////////////////////////////////////////////////////////
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
						   VertexBase* vrtCenter, const vector3& normal,
						   number radius, AInt& aInt, APosition& aPos)
{
//	defines when a point is considered to be close to the rim
//	CLOSE_TO_RIM has to be higher than 0!
	const number CLOSE_TO_RIM = 0.8;

//	the quality factor is a lower border for quality decrease
	const number QUALITY_FACTOR = 0.6;

//	the radius has to be big enough
	if(radius < SMALL)
		return false;

//	get the position attachment accessor
	if(!grid.has_vertex_attachment(aPos))
		return false;

//	make sure all required options are enabled in grid
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		LOG("  INFO in AdaptSurfaceGridToCylinder: auto-enabling FACEOPT_AUTOGENERATE_EDGES.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	position accessor
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	position of the center
	vector3 center = aaPos[vrtCenter];

//	alias for the selector
	Selector& sel = selOut;
//	clear the selector
	sel.clear();
//	store whether selection_inheritance was enabled - if not we'll enable it now
	bool selInheritanceWasEnabled = sel.selection_inheritance_enabled();
	sel.enable_selection_inheritance(true);

//	we'll store faces that lie inside the cylinder in this vector
	vector<Face*> vCylinderFaces;

//	we'll collect faces in this vector
	vector<Face*> vFaces;

//	we'll collect edges in this vector
	vector<EdgeBase*> vEdges;

//	we'll store normals of the associated faces in this vector
	vector<vector3>	vNormals;

//	begin marking
	grid.begin_marking();

/*
//	this queue holds vertices which are checked in the next iteration
	queue<VertexBase*> qVrts;
	qVrts.push(vrtCenter);
*/
//	this vectot contains the verices on which the algorithm works.
//	it is extended during the algo.
	vector<VertexBase*> vVrts;
	vVrts.push_back(vrtCenter);
	grid.mark(vrtCenter);

//	those indices define the area in vVrts on which we work
	size_t iFirst = 0;
	size_t iEnd = 1;

	while(iFirst != iEnd)
	{
	//	as long as there are vertices in the queue we have to iterate
		//while(!qVrts.empty())
		for(int vrtInd = iFirst; vrtInd < iEnd; ++ vrtInd)
		{
		//	clear the selector
			sel.clear();

		//	the vertex that we'll check
			//VertexBase* vrt = qVrts.front();
			//qVrts.pop();
			VertexBase* vrt = vVrts[vrtInd];

		//	collect associated faces - they will be needed twice
			vFaces.clear();
			for(FaceIterator iter = grid.associated_faces_begin(vrt);
				iter != grid.associated_faces_end(vrt); ++iter)
				vFaces.push_back(*iter);

		//	calculate normals
			vNormals.resize(vFaces.size());
			for(size_t i = 0; i < vFaces.size(); ++i)
				CalculateNormal(vNormals[i], vFaces[i], aaPos);

		//	iterate over all associated edges of vrt
			EdgeBaseIterator iterEnd = grid.associated_edges_end(vrt);
			for(EdgeBaseIterator iter = grid.associated_edges_begin(vrt);
				iter != iterEnd; ++iter)
			{
				EdgeBase* e = *iter;
			//	get the connected vertex and mark it - if it was alreay marked, then ignore it
				VertexBase* cv = GetConnectedVertex(e, vrt);
				if(grid.is_marked(cv))
					continue;
				grid.mark(cv);

			//	the position of the connected vertex
				vector3 cpos = aaPos[cv];
			//	distance of the conencted vertex to the central axis of the cylinder.
				vector3 projPos;
				ProjectPointToRay(projPos, cpos, center, normal);
				vector3 dir;
				VecSubtract(dir, cpos, projPos);
				number dist = VecLength(dir);

			//	if the connected vertex is close to the boundary of the cylinder
			//	or outside, we'll try to project it onto the cylinder.
			//	We'll use two measures to avoid geometry-conflicts.
			//	the first one guarantees that the quality of the area will not get
			//	too bad. The second makes sure that no triangles are flipped.
				if(dist / radius > CLOSE_TO_RIM){
					number initialQuality = AreaFaceQuality(vFaces.begin(), vFaces.end(), aaPos);
				//	move the vertex to the rim
					VecScale(dir, dir, radius / dist);
					VecAdd(aaPos[cv], projPos, dir);

					bool passedTests = true;

				//	compare qualities
					if(AreaFaceQuality(vFaces.begin(), vFaces.end(), aaPos)
						< QUALITY_FACTOR * initialQuality)
						passedTests = false;

				//	compare normals
					for(size_t i = 0; i < vFaces.size(); ++i){
						vector3 n;
						CalculateNormal(n, vFaces[i], aaPos);
						if(VecDot(n, vNormals[i]) < QUALITY_FACTOR){
							passedTests = false;
							break;
						}
					}
					
				//	if one of the tests failed, undo the movement
					if(!passedTests){
						aaPos[cv] = cpos;
					//	if the vertex lies outside, we'll have to split the edge.
					//	if not we'll check it in the next iteration.
						if(dist > radius){
							sel.select(e);
							LOG("1");
						}
						else{
							vVrts.push_back(cv);
							LOG("2");
						}
					}
				}
				else{
				//	check the vertex in the next iteration
					vVrts.push_back(cv);
					LOG("3");
				}
			}

		//	if there are very large edges that cut the cylinder,
		//	we'll split them.
		//	iterate over associated faces
			for(size_t i = 0; i < vFaces.size(); ++i){
				Face* f = vFaces[i];
				if(!grid.is_marked(f)){
					grid.mark(f);
					vCylinderFaces.push_back(f);
				//	collect associated edges
					CollectEdges(vEdges, grid, f);
					for(size_t j = 0; j < vEdges.size(); ++j){
						EdgeBase* e = vEdges[j];
						if(!grid.is_marked(e)){
						//	if both endpoints are marked
							if(grid.is_marked(e->vertex(0)) && grid.is_marked(e->vertex(1)))
							{
								grid.mark(e);
							//	check whether the edge should be splitted
								if(VecDistance(aaPos[e->vertex(0)], aaPos[e->vertex(1)]) > radius * 1.05)
								{
									sel.select(e);
									LOG("4");
								}
							}
						}
					}
				}
			}

		//	refine marked edges
			if(sel.num<EdgeBase>() > 0){
				if(!Refine(grid, sel, aInt)){
				//	clean up
					LOG("Refine failed in AdaptSurfaceGridToCylinder\n");
					grid.end_marking();
					sel.enable_selection_inheritance(selInheritanceWasEnabled);
					return false;
				}
			//	if new vertices have been created, they are now selected
				for(VertexBaseIterator iter = sel.begin<VertexBase>();
					iter != sel.end<VertexBase>(); ++iter)
				{
				//	if the vertex is outside of the cylinder, we'll project it
				//	back on the cylinder.
				//	this operation is save (although it might produce triangles
				//	of bad quality).
					vector3 projPos;
					ProjectPointToRay(projPos, aaPos[*iter], center, normal);
					vector3 dir;
					VecSubtract(dir, aaPos[*iter], projPos);
					number dist = VecLength(dir);
					if(dist > radius){
						VecScale(dir, dir, radius / dist);
						VecAdd(aaPos[*iter], projPos, dir);
					}

				//	add the new vertex to vVrts and mark it
					vVrts.push_back(*iter);
					grid.mark(*iter);
				}
			}

		//	adjust iFirst and iEnd
			iFirst = iEnd;
			iEnd = vVrts.size();
		}
	}

//	select all faces in vCylinderFaces
	sel.clear();
	sel.select(vCylinderFaces.begin(), vCylinderFaces.end());

//	clean up
	grid.end_marking();
	sel.enable_selection_inheritance(selInheritanceWasEnabled);

//	done
	return true;
}

}//	end of namespace
