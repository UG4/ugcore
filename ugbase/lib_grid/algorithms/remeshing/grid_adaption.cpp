#include <queue>
#include <stack>
#include "lib_grid/lg_base.h"
#include "../geom_obj_util/geom_obj_util.h"
#include "../refinement/regular_refinement.h"
#include "../refinement/refinement_projectors/misc_refinement_projectors.h"
#include "grid_adaption.h"

using namespace std;

namespace ug
{

bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
							   Vertex* vrtCenter, const vector3& normal,
							   number radius, number rimSnapThreshold,  AInt& aInt,
							   APosition& aPos)
{
	if(!grid.has_vertex_attachment(aPos)){
		UG_THROW("Position attachment required!");
	}

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

	if(rimSnapThreshold < 0)
		rimSnapThreshold = 0;

	if(rimSnapThreshold > (radius - SMALL))
		rimSnapThreshold = radius - SMALL;

	const number smallRadius = radius - rimSnapThreshold;
	const number smallRadiusSq = smallRadius * smallRadius;
	const number largeRadius = radius + rimSnapThreshold;
	const number largeRadiusSq = largeRadius * largeRadius;

//	the cylinder geometry
	vector3 axis;
	VecNormalize(axis, normal);
	vector3 center = aaPos[vrtCenter];

//	recursively select all vertices in the cylinder which can be reached from a
//	selected vertex by following an edge. Start with the given one.
//	We'll also select edges which connect inner with outer vertices. Note that
//	some vertices are considered rim-vertices (those with a distance between
//	smallRadius and largeRadius). Those are neither considered inner nor outer.
	Selector& sel = selOut;
	sel.clear();
	sel.select(vrtCenter);

	stack<Vertex*> vrtStack;
	vrtStack.push(vrtCenter);

	Grid::edge_traits::secure_container edges;
	Grid::face_traits::secure_container faces;
	vector<Quadrilateral*> quads;

	while(!vrtStack.empty()){
		Vertex* curVrt = vrtStack.top();
		vrtStack.pop();

	//	we have to convert associated quadrilaterals to triangles.
	//	Be careful not to alter the array of associated elements while we iterate
	//	over it...
		quads.clear();
		grid.associated_elements(faces, curVrt);
		for(size_t i = 0; i < faces.size(); ++i){
			if(faces[i]->num_vertices() == 4){
				Quadrilateral* q = dynamic_cast<Quadrilateral*>(faces[i]);
				if(q)
					quads.push_back(q);
			}
		}

		for(size_t i = 0; i < quads.size(); ++i){
			Triangulate(grid, quads[i], &aaPos);
		}

	//	now check whether edges leave the cylinder and mark them accordingly.
	//	Perform projection of vertices to the cylinder rim for vertices which
	//	lie in the threshold area.
		grid.associated_elements(edges, curVrt);

		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
			Edge* e = edges[i_edge];
			Vertex* vrt = GetConnectedVertex(e, curVrt);

			if(sel.is_selected(vrt))
				continue;

			vector3 p = aaPos[vrt];
			vector3 proj;
			ProjectPointToRay(proj, p, center, axis);
			number distSq = VecDistanceSq(p, proj);

			if(distSq < smallRadiusSq){
				sel.select(vrt);
				vrtStack.push(vrt);
			}
			else if(distSq < largeRadiusSq){
				sel.select(vrt);
			//	cut the ray from center through p with the cylinder hull to calculate
			//	the new position of vrt.
				vector3 dir;
				VecSubtract(dir, p, center);
				number t0, t1;
				if(RayCylinderIntersection(t0, t1, center, dir, center, axis, radius))
					VecScaleAdd(aaPos[vrt], 1., center, t1, dir);
			}
			else{
			//	the edge will be refined later on
				sel.select(e);
			}
		}
	}

//	refine selected edges and use a special refinement callback, which places
//	new vertices on edges which intersect a cylinder on the cylinders hull.
	RefinementCallback_IntersectCylinder refCallback(grid, center, axis, radius, aPos);
	Refine(grid, sel, aInt, &refCallback);

//	finally select all triangles which lie in the cylinder
	sel.clear();
	vrtStack.push(vrtCenter);
	sel.select(vrtCenter);

	while(!vrtStack.empty()){
		Vertex* curVrt = vrtStack.top();
		vrtStack.pop();
		grid.associated_elements(faces, curVrt);

		for(size_t i_face = 0; i_face < faces.size(); ++i_face){
			Face* f = faces[i_face];
			if(sel.is_selected(f))
				continue;

			sel.select(f);

			for(size_t i = 0; i < f->num_vertices(); ++i){
				Vertex* vrt = f->vertex(i);
				if(!sel.is_selected(vrt)){
					number dist = DistancePointToRay(aaPos[vrt], center, axis);
					if(dist < (radius - SMALL)){
						sel.select(vrt);
						vrtStack.push(vrt);
					}
				}
			}
		}
	}

	sel.clear<Vertex>();

	return true;
}

////////////////////////////////////////////////////////////////////////
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
						   Vertex* vrtCenter, const vector3& normal,
						   number radius, number rimSnapThreshold, APosition& aPos)
{
	AInt aInt;
	grid.attach_to_edges(aInt);
	bool retVal = AdaptSurfaceGridToCylinder(selOut, grid, vrtCenter, normal,
											radius, rimSnapThreshold, aInt, aPos);
	grid.detach_from_edges(aInt);
	return retVal;
}

}//	end of namespace
