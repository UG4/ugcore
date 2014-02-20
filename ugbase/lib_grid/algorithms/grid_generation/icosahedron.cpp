// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 29.02.2012 (m,d,y)
 
#include "icosahedron.h"
#include "lib_grid/algorithms/refinement/regular_refinement.h"
#include "../refinement/refinement_projectors/standard_refinement_projectors.h"

namespace ug{

///	Creates an Icosahedron
void GenerateIcosahedron(Grid& grid, const vector3& center,
						 number radius, AVector3& aPos)
{
	const number A = 0.85065080835204;
	const number B = 0.525731112119134;

//	create the vertices
	const number coords[12][3] = {	{-B, A, 0}, {0, B, A}, {B, A, 0}, {0, B, -A},
									{-A, 0, B}, {A, 0, B}, {A, 0, -B}, {-A, 0, -B},
									{-B, -A, 0}, {0, -B, A}, {B, -A, 0}, {0, -B, -A}};

	const int inds[20][3] = {	{0, 1, 2}, {0, 2, 3},
								{3, 7, 0}, {7, 4, 0}, {0, 4, 1},
								{1, 5, 2}, {5, 6, 2}, {2, 6, 3},
								{4, 9, 1}, {1, 9, 5}, {7, 3, 11}, {3, 6, 11},
								{11, 8, 7}, {7, 8, 4}, {8, 9, 4},
								{9, 10, 5}, {10, 6, 5}, {10, 11, 6},
								{8, 10, 9}, {8, 11, 10}};

	Vertex* vrts[12];

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::AttachmentAccessor<Vertex, AVector3> aaPos(grid, aPos);

	for(size_t i = 0; i < 12; ++i){
		vrts[i] = *grid.create<RegularVertex>();
		aaPos[vrts[i]] = vector3(coords[i][0], coords[i][1], coords[i][2]);
		VecScaleAdd(aaPos[vrts[i]], 1.0, center, radius, aaPos[vrts[i]]);
	}

//	construct triangles
	for(size_t i = 0; i < 20; ++i){
		grid.create<Triangle>(TriangleDescriptor(vrts[inds[i][0]], vrts[inds[i][1]],
												 vrts[inds[i][2]]));
	}
}

void GenerateIcosphere(Grid& grid, const vector3& center, number radius,
					   int numRefinements, AVector3& aPos, Selector* psel)
{
	Selector defSel;
	if(!psel){
		defSel.assign_grid(grid);
		psel = &defSel;
	}

	Selector& sel = *psel;

//	enable autoselection
	bool autoselectionEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	clear the selector
	sel.clear();

//	create an icosahedron. All elements of the sphere will be selected, since we
//	enabled autoselection
	GenerateIcosahedron(grid, center, radius, aPos);

//	perform refinement steps
//	we need a temporary int attachment for the edges
	AInt aInt;
	grid.attach_to_edges(aInt);

//	use a refinement callback to project the new vertices to the sphere
	RefinementCallbackSphere<APosition> sphereProjecton(grid, aPos, center, radius);

	for(int i = 0; i < numRefinements; ++i)
		Refine(grid, sel, aInt, &sphereProjecton);

//	remove attachments
	grid.detach_from_edges(aInt);

//	restore autoselection
	sel.enable_autoselection(autoselectionEnabled);
}

}// end of namespace
