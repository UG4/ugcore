// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 10.02.2011 (m,d,y)

#include <iostream>
#include <sstream>
#include <vector>

#include "registry/registry.h"
#include "bridge.h"

#include "common/profiler/profiler.h"
#include "lib_grid/lib_grid.h"

#include "lib_disc/domain.h"
#include "lib_disc/domain_util.h"

#include "lib_disc/parallelization/domain_distribution.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/load_balancing.h"
	#include "lib_grid/parallelization/parallelization.h"
#endif

// ONLY TEMPORARY
#include "lib_grid/visualization/grid_visualization.h"

using namespace std;

namespace ug{

//	This method is only a temporary test method and will be replaced by
//	a more sophisticated approach
template <typename TDomain>
static void MinimizeMemoryFootprint(TDomain& dom)
{
	dom.grid()->set_options(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION
							 | GRIDOPT_AUTOGENERATE_SIDES);
}

template <typename TDomain>
static void LoadAndRefineDomain(TDomain& domain, const char* filename,
								int numRefs)
{
	try{
		LoadDomain(domain, filename);
	}
	UG_CATCH_THROW("LoadAndRefineDomain: Could not load domain at file: "<<filename);

	GlobalMultiGridRefiner ref(*domain.grid());
	for(int i = 0; i < numRefs; ++i)
		ref.refine();
}

template <typename TDomain>
static bool SavePartitionMap(PartitionMap& pmap, TDomain& domain,
							 const char* filename)
{
	if(domain.grid().get() != pmap.get_partition_handler().grid())
	{
		UG_LOG("WARNING in SavePartitionMap: The given partition map was not"
				" created for the given domain. Aborting...\n");
		return false;
	}

	return SavePartitionMapToFile(pmap, filename, domain.position_attachment());
}


template <typename TDomain>
static bool TestDomainInterfaces(TDomain* dom)
{
	#ifdef UG_PARALLEL
		return TestGridLayoutMap(*dom->grid(),
					dom->distributed_grid_manager()->grid_layout_map());
	#endif
	return true;
}

//	ONLY TEMPORARY!!! DO NOT USE THIS METHOD!!!
template <typename TDomain>
static void TestDomainVisualization(TDomain& dom)
{
	GridVisualization<number, int, typename TDomain::position_attachment_type> gridVis;
	gridVis.set_grid(*dom.grid(), dom.position_attachment());
	gridVis.update_visuals();

//	log vertex positions
	UG_LOG("GridVisualization:\n");
	UG_LOG("Vertices (" << gridVis.num_vertices() << "):");
	const number* pos = gridVis.vertex_positions();
	for(int i = 0; i < gridVis.num_vertices(); ++i){
		UG_LOG(" (" <<  pos[3*i] << ", " << pos[3*i+1] << ", " << pos[3*i+2] << ")");
	}
	UG_LOG("\n\n");

//	iterate over all visuals
	for(int i_vis = 0; i_vis < gridVis.num_visuals(); ++i_vis){
		UG_LOG("Visual " << i_vis << ":\n");
	//	log the triangles of this visual (triple of vertex indices.)
		UG_LOG("Triangles:");
		const int* tris = gridVis.triangle_list(i_vis);
		for(int i = 0; i < gridVis.num_triangles(i_vis); ++i){
			UG_LOG(" [" <<  tris[3*i] << ", " << tris[3*i+1] << ", " << tris[3*i+2] << "]");
		}
		UG_LOG("\n\n");

	//	log the triangle normals of this visual
		UG_LOG("Tri-Normals:");
		const number* norms = gridVis.face_normals(i_vis);

		for(int i = 0; i < gridVis.num_triangles(i_vis); ++i){
			UG_LOG(" (" <<  norms[3*i] << ", " << norms[3*i+1] << ", " << norms[3*i+2] << ")");
		}
		UG_LOG("\n\n");
	}
}


template <typename TDomain>
static void ScaleDomain(TDomain& dom, number sx, number sy, number sz)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	typename TDomain::grid_type& g = *dom.grid();
	vector3 s(sx, sy, sz);

	const int numCoords = TDomain::position_type::Size;
	UG_ASSERT(numCoords <= 3, "too many coordinates.");

	for(VertexBaseIterator iter = g.vertices_begin();
		iter != g.vertices_end(); ++iter)
	{
		for(int i = 0; i < numCoords; ++i)
			aaPos[*iter][i] *= s[i];
	}
}


template <typename TDomain>
static void TranslateDomain(TDomain& dom, number tx, number ty, number tz)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	typename TDomain::grid_type& g = *dom.grid();
	vector3 t(tx, ty, tz);

	const int numCoords = TDomain::position_type::Size;
	UG_ASSERT(numCoords <= 3, "too many coordinates.");

	for(VertexBaseIterator iter = g.vertices_begin();
		iter != g.vertices_end(); ++iter)
	{
		for(int i = 0; i < numCoords; ++i)
			aaPos[*iter][i] += t[i];
	}
}



namespace bridge{

template <typename TDomain>
static bool RegisterDomainInterface_(Registry& reg, string grp)
{
//	the dimension suffix
	string dimSuffix = GetDomainSuffix<TDomain>();

//	the dimension tag
	string dimTag = GetDomainTag<TDomain>();

//	Domain
	{
		string name = string("Domain").append(dimSuffix);
		reg.add_class_<TDomain>(name, grp)
			.add_constructor()
			.add_method("subset_handler", static_cast<SmartPtr<MGSubsetHandler> (TDomain::*)()>(&TDomain::subset_handler))
			.add_method("grid", static_cast<SmartPtr<MultiGrid> (TDomain::*)()>(&TDomain::grid))
			.add_method("get_dim", static_cast<int (TDomain::*)() const>(&TDomain::get_dim))
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "Domain", dimTag);
	}

// 	LoadDomain
	reg.add_function("LoadDomain", static_cast<void (*)(TDomain&, const char*)>(
					 &LoadDomain<TDomain>), grp,
					"", "Domain # Filename | load-dialog | endings=[\"ugx\"]; description=\"*.ugx-Files\" # Number Refinements",
					"Loads a domain", "No help");

//	LoadAndRefineDomain
	reg.add_function("LoadAndRefineDomain", &LoadAndRefineDomain<TDomain>, grp,
					"", "Domain # Filename # NumRefines | load-dialog | endings=[\"ugx\"]; description=\"*.ugx-Files\" # Number Refinements",
					"Loads a domain and performs global refinement", "No help");
//	SaveDomain
	reg.add_function("SaveDomain", &SaveDomain<TDomain>, grp,
					"", "Domain # Filename|save-dialog",
					"Saves a domain", "No help");

//	SavePartitionMap
	reg.add_function("SavePartitionMap", &SavePartitionMap<TDomain>, grp,
					"Success", "PartitionMap # Domain # Filename|save-dialog",
					"Saves a partition map", "No help");

//	DistributeDomain
	reg.add_function("DistributeDomain", static_cast<bool (*)(TDomain&)>(
					 &DistributeDomain<TDomain>), grp);

	reg.add_function("DistributeDomain", static_cast<bool (*)(TDomain&, PartitionMap&)>(
					 &DistributeDomain<TDomain>), grp);

//	todo: remove this
	{
		string name = string("DistributeDomain").append(dimSuffix);
		reg.add_function(name.c_str(), static_cast<bool (*)(TDomain&)>(
						 &DistributeDomain<TDomain>), grp);
	}

	reg.add_function("PartitionDomain_Bisection",
					 &PartitionDomain_Bisection<TDomain>, grp);

	reg.add_function("PartitionDomain_MetisKWay",
					 &PartitionDomain_MetisKWay<TDomain>, grp);

	reg.add_function("RedistributeDomain",
					 &RedistributeDomain<TDomain>, grp);

//	transform the domain
	reg.add_function("ScaleDomain", &ScaleDomain<TDomain>, grp);
	reg.add_function("TranslateDomain", &TranslateDomain<TDomain>, grp);

//	debugging
	reg.add_function("TestDomainInterfaces", &TestDomainInterfaces<TDomain>, grp);

//	ONLY TEMPORARY
	reg.add_function("TestDomainVisualization", &TestDomainVisualization<TDomain>, grp);
	reg.add_function("MinimizeMemoryFootprint", &MinimizeMemoryFootprint<TDomain>, grp);

	return true;
}

///	methods that are only available for 2d and 3d are registered here
template <typename TDomain>
static bool RegisterDomainInterface_2d_3d(Registry& reg, string grp)
{
	reg.add_function("PartitionDomain_RegularGrid",
					 &PartitionDomain_RegularGrid<TDomain>, grp);

	return true;
}

bool RegisterDomainInterface(Registry& reg, string parentGroup)
{
	bool bSuccess = true;

	string grp = parentGroup; grp.append("/Domain");

#ifdef UG_DIM_1
	bSuccess &= RegisterDomainInterface_<Domain<1, MultiGrid, MGSubsetHandler> >(reg, grp);
#endif
#ifdef UG_DIM_2
	bSuccess &= RegisterDomainInterface_<Domain<2, MultiGrid, MGSubsetHandler> >(reg, grp);
	bSuccess &= RegisterDomainInterface_2d_3d<Domain<2, MultiGrid, MGSubsetHandler> >(reg, grp);
#endif
#ifdef UG_DIM_3
	bSuccess &= RegisterDomainInterface_<Domain<3, MultiGrid, MGSubsetHandler> >(reg, grp);
	bSuccess &= RegisterDomainInterface_2d_3d<Domain<3, MultiGrid, MGSubsetHandler> >(reg, grp);
#endif
	return bSuccess;
}

}// end of namespace
}// end of namespace
