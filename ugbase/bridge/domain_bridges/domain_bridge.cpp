// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 10.02.2011 (m,d,y)

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

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
namespace Domain{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Domain
	{
		typedef IDomain<> TBase;
		string name = string("Domain").append(suffix);
		reg.add_class_<TDomain, TBase>(name, grp)
			.add_constructor()
//			.add_method("subset_handler", static_cast<SmartPtr<MGSubsetHandler> (TDomain::*)()>(&TDomain::subset_handler))
//			.add_method("grid", static_cast<SmartPtr<MultiGrid> (TDomain::*)()>(&TDomain::grid))
//			.add_method("get_dim", static_cast<int (TDomain::*)() const>(&TDomain::get_dim))
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "Domain", tag);
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
		string name = string("DistributeDomain").append(suffix);
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
}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
//	IDomain
	{
		typedef IDomain<> T;
		reg.add_class_<T>("IDomain", grp)
			.add_method("subset_handler", static_cast<SmartPtr<MGSubsetHandler> (T::*)()>(&T::subset_handler))
			.add_method("grid", static_cast<SmartPtr<MultiGrid> (T::*)()>(&T::grid))
			.add_method("get_dim", static_cast<int (T::*)() const>(&T::get_dim))
			.set_construct_as_smart_pointer(true);
	}

}

}; // end Functionality


///	methods that are only available for 2d and 3d are registered here
struct Functionality2d3d
{
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	reg.add_function("PartitionDomain_RegularGrid",
					 &PartitionDomain_RegularGrid<TDomain>, grp);
}
}; // end Functionality2d3d

}// end Domain

void RegisterBridge_Domain(Registry& reg, string grp)
{
	grp.append("/Domain");

	typedef Domain::Functionality Functionality;
	typedef boost::mpl::list<
	#ifdef UG_DIM_2
			Domain2d
	#endif
	#if defined UG_DIM_2 && defined UG_DIM_3
			,
	#endif
	#ifdef UG_DIM_3
			Domain3d
	#endif
	> CompileDomain2d3dList;

	try{
		RegisterDomainDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Domain::Functionality2d3d, CompileDomain2d3dList>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace
}// end of namespace
