// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.07.2011 (m,d,y)
 
#include <vector>
#include <string>
#include <sstream>
#include "ug_bridge.h"
#include "lib_discretization/domain.h"
#include "lib_grid/lib_grid.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
///	Creates a global domain refiner.
/**	Automatically chooses whether a parallel refiner is required.*/
template <typename TDomain>
static IRefiner* GlobalDomainRefiner(TDomain* dom)
{
//todo: support normal grids, too!
	#ifdef UG_PARALLEL
		if(pcl::GetNumProcesses() > 1){
			return new ParallelGlobalRefiner_MultiGrid(*dom->get_distributed_grid_manager());
		}
	#endif

	return new GlobalMultiGridRefiner(dom->get_grid());
}

////////////////////////////////////////////////////////////////////////////////
///	Creates an adaptive hanging node domain refiner.
/**	Automatically chooses whether a parallel refiner is required.*/
template <typename TDomain>
static IRefiner* HangingNodeDomainRefiner(TDomain* dom)
{
//todo: support normal grids, too!
	#ifdef UG_PARALLEL
		if(pcl::GetNumProcesses() > 1){
			return new ParallelHangingNodeRefiner_MultiGrid(*dom->get_distributed_grid_manager());
		}
	#endif

	return new HangingNodeRefiner_MultiGrid(dom->get_grid());
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements from refinement.
/**	If used in a parallel environment only elements on the calling procs
 * are marked.
 */
static void MarkForRefinement_All(IRefiner& ref)
{
	Grid* g = ref.get_associated_grid();
	if(!g){
		UG_LOG("Refiner is not registered at a grid. Aborting.\n");
		return;
	}
	ref.mark(g->vertices_begin(), g->vertices_end());
	ref.mark(g->edges_begin(), g->edges_end());
	ref.mark(g->faces_begin(), g->faces_end());
	ref.mark(g->volumes_begin(), g->volumes_end());
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all vertices in the given d-dimensional sphere.
template <class TDomain>
void MarkForRefinement_VerticesInSphere(TDomain& dom, IRefiner& refiner,
									const typename TDomain::position_type& center,
									number radius)
{
	typedef typename TDomain::position_type 			position_type;
	typedef typename TDomain::position_accessor_type	position_accessor_type;

//	make sure that the refiner was created for the given domain
	if(refiner.get_associated_grid() != &dom.get_grid()){
		throw(UGError("ERROR in MarkForRefinement_VerticesInSphere: "
					"Refiner was not created for the specified domain. Aborting."));
	}

	Grid& grid = *refiner.get_associated_grid();
	position_accessor_type& aaPos = dom.get_position_accessor();

//	we'll compare against the square radius.
	number radiusSq = radius * radius;

//	we'll store associated edges, faces and volumes in those containers
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all vertices of the grid. If a vertex is inside the given sphere,
//	then we'll mark all associated elements.
	for(VertexBaseIterator iter = grid.begin<VertexBase>();
		iter != grid.end<VertexBase>(); ++iter)
	{
		if(VecDistanceSq(center, aaPos[*iter]) <= radiusSq){
			CollectAssociated(vEdges, grid, *iter);
			CollectAssociated(vFaces, grid, *iter);
			CollectAssociated(vVols, grid, *iter);

			refiner.mark(vEdges.begin(), vEdges.end());
			refiner.mark(vFaces.begin(), vFaces.end());
			refiner.mark(vVols.begin(), vVols.end());
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements which lie completely in the given d-dimensional sphere.
/**	Valid parameters for TElem are EdgeBase, Face, Volume.*/
template <class TDomain, class TElem>
void MarkForRefinement_ElementsInSphere(TDomain& dom, IRefiner& refiner,
									const typename TDomain::position_type& center,
									number radius)
{
	typedef typename TDomain::position_type 			position_type;
	typedef typename TDomain::position_accessor_type	position_accessor_type;
	typedef typename geometry_traits<TElem>::iterator	ElemIter;

//	make sure that the refiner was created for the given domain
	if(refiner.get_associated_grid() != &dom.get_grid()){
		throw(UGError("ERROR in MarkForRefinement_VerticesInCube: "
					"Refiner was not created for the specified domain. Aborting."));
	}

	Grid& grid = *refiner.get_associated_grid();
	position_accessor_type& aaPos = dom.get_position_accessor();

//	we'll compare against the square radius.
	number radiusSq = radius * radius;

//	we'll store associated edges, faces and volumes in those containers
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all elements of the grid. If all vertices of an element are inside
//	the given sphere, then we'll mark those elements.
	for(ElemIter iter = grid.begin<TElem>();
		iter != grid.end<TElem>(); ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	bool flag to check whether all vertices are in the sphere
		bool bInSphere = true;

	//	loop all vertices
		for(size_t i = 0; i < elem->num_vertices(); ++i)
		{
		//	check if vertex is in sphere
			if(VecDistanceSq(center, aaPos[elem->vertex(i)]) > radiusSq)
				bInSphere = false;
		}

	//	mark the element
		if(bInSphere)
			refiner.mark(elem);
	}
}

////////////////////////////////////////////////////////////////////////////////
///	Marks all elements which have vertices in the given d-dimensional cube.
/**	Make sure that TAPos is an attachment of vector_t position types.*/
template <class TDomain>
void MarkForRefinement_VerticesInCube(TDomain& dom, IRefiner& refiner,
									const typename TDomain::position_type& min,
									const typename TDomain::position_type& max)
{
	typedef typename TDomain::position_type 			position_type;
	typedef typename TDomain::position_accessor_type	position_accessor_type;

//	make sure that the refiner was created for the given domain
	if(refiner.get_associated_grid() != &dom.get_grid()){
		throw(UGError("ERROR in MarkForRefinement_VerticesInCube: "
					"Refiner was not created for the specified domain. Aborting."));
	}

	Grid& grid = *refiner.get_associated_grid();
	position_accessor_type& aaPos = dom.get_position_accessor();

//	we'll store associated edges, faces and volumes in those containers
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	iterate over all vertices of the grid. If a vertex is inside the given cube,
//	then we'll mark all associated elements.
	for(VertexBaseIterator iter = grid.begin<VertexBase>();
		iter != grid.end<VertexBase>(); ++iter)
	{
	//	Position
		position_type& pos = aaPos[*iter];

	//	check flag
		bool bRefine = true;

	//	check node
		for(size_t d = 0; d < pos.size(); ++d)
			if(pos[d] < min[d] || max[d] < pos[d])
				bRefine = false;

		if(bRefine)
		{
			CollectAssociated(vEdges, grid, *iter);
			CollectAssociated(vFaces, grid, *iter);
			CollectAssociated(vVols, grid, *iter);

			refiner.mark(vEdges.begin(), vEdges.end());
			refiner.mark(vFaces.begin(), vFaces.end());
			refiner.mark(vVols.begin(), vVols.end());
		}
	}
}



namespace bridge{

static bool RegisterRefinementBridge_DomIndep(Registry& reg, string parentGroup)
{
	try
	{
	//	get group string
		stringstream groupString; groupString << parentGroup << "/Refinement";
		string grp = groupString.str();

	//	register domain independent mark methods
		reg.add_function("MarkForRefinement_All", &MarkForRefinement_All, grp);
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterRefinementBridge_DomIndep: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////
template <class TDomain>
static bool RegisterRefinementBridge_DomDep(Registry& reg, string parentGroup)
{
	typedef TDomain 							domain_type;
	typedef typename TDomain::position_type		pos_type;

	try
	{
	//	get group string
		stringstream groupString; groupString << parentGroup << "/Refinement";
		string grp = groupString.str();

	//	refiner factory-method registration
	//	Note that the refiners themselfs have already been registered in lib_grid_bridge.
		reg.add_function("GlobalDomainRefiner",
						 &GlobalDomainRefiner<domain_type>, grp);
		reg.add_function("HangingNodeDomainRefiner",
						 &HangingNodeDomainRefiner<domain_type>, grp);

	//	register domain dependent mark methods
		reg.add_function("MarkForRefinement_VerticesInSphere",
					&MarkForRefinement_VerticesInSphere<domain_type>, grp)
			.add_function("MarkForRefinement_EdgesInSphere",
					&MarkForRefinement_ElementsInSphere<domain_type, EdgeBase>, grp)
			.add_function("MarkForRefinement_FacesInSphere",
					&MarkForRefinement_ElementsInSphere<domain_type, Face>, grp)
			.add_function("MarkForRefinement_VolumesInSphere",
					&MarkForRefinement_ElementsInSphere<domain_type, Volume>, grp)
			.add_function("MarkForRefinement_VerticesInCube",
					&MarkForRefinement_VerticesInCube<domain_type>, grp);

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterRefinementBridge_DomDep: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////
bool RegisterRefinementBridge(Registry& reg, string parentGroup)
{
	bool bSuccess = true;

	bSuccess &= RegisterRefinementBridge_DomIndep(reg, parentGroup);

#ifdef UG_DIM_1
	bSuccess &= RegisterRefinementBridge_DomDep<Domain<1, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	bSuccess &= RegisterRefinementBridge_DomDep<Domain<2, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	bSuccess &= RegisterRefinementBridge_DomDep<Domain<3, MultiGrid, MGSubsetHandler> >(reg, parentGroup);
#endif
	return bSuccess;
}

}// end of namespace
}// end of namespace
