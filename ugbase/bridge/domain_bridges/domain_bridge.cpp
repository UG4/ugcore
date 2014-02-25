// created by Andreas Vogel, Sebastian Reiter
// s.b.reiter@googlemail.com
// 10.02.2011 (m,d,y)

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#include "common/profiler/profiler.h"

#include "lib_disc/domain.h"
#include "lib_disc/domain_util.h"

#include "lib_disc/parallelization/domain_distribution.h"

#include "lib_grid/algorithms/refinement/global_multi_grid_refiner.h"
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"

#include "lib_grid/algorithms/subset_util.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/load_balancing.h"
	#include "lib_grid/parallelization/parallelization.h"
	#include "lib_disc/parallelization/domain_load_balancer.h"
#endif

#include "lib_grid/parallelization/util/partition_weighting_callbacks.h"

using namespace std;

namespace ug{

/**
 * \defgroup domain_bridge Domain Bridge
 * \ingroup bridge
 * \{
 */

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
	PROFILE_FUNC_GROUP("grid");
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
	PROFILE_FUNC_GROUP("grid");
	if(domain.grid().get() != pmap.get_partition_handler()->grid())
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


template <typename TDomain>
static void ScaleDomain(TDomain& dom, number sx, number sy, number sz)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	typename TDomain::grid_type& g = *dom.grid();
	vector3 s(sx, sy, sz);

	const int numCoords = TDomain::position_type::Size;
	UG_ASSERT(numCoords <= 3, "too many coordinates.");

	for(VertexIterator iter = g.vertices_begin();
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

	for(VertexIterator iter = g.vertices_begin();
		iter != g.vertices_end(); ++iter)
	{
		for(int i = 0; i < numCoords; ++i)
			aaPos[*iter][i] += t[i];
	}
}

/**
 * Calculates the area sum of faces in given domain and subset handler sh.
 * Please note that the sum is returned for faces in subset with index si and
 * on grid level lvl in the domain dom.
 *
 * \param dom domain
 * \param sh subset handler
 * \param si subset index
 * \param lvl grid level
 *
 * \return \c number area sum
 */
template <typename TDomain>
static number FaceArea(TDomain& dom, ISubsetHandler& sh, int si, size_t lvl)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	UG_ASSERT(TDomain::position_type::Size <= 3, "too many coordinates.");

	return FaceArea(sh, si, lvl, aaPos);
}

/**
 * Calculates the area sum of faces in given domain with implicit given
 * subset handler. Please note the sum is returned for faces in subset with
 * index si and on grid level lvl in the domain dom.
 *
 * \param dom domain
 * \param si subset index
 * \param lvl grid level
 *
 * \return \c number area sum
 */
template <typename TDomain>
static number FaceArea(TDomain& dom, int si, size_t lvl)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	UG_ASSERT(TDomain::position_type::Size <= 3, "too many coordinates.");

return FaceArea(*dom.subset_handler(), si, lvl, aaPos);
}

/**
 * Calculates the area sum of faces in a given domain with implicit given
 * subset handler. Please note the sum is returned for faces in subset with
 * index si on grid level 0 in the domain dom.
 *
 * \param dom domain
 * \param si subset index
 *
 * \return \c number area sum
 */
template <typename TDomain>
static number FaceArea(TDomain& dom, int si)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	UG_ASSERT(TDomain::position_type::Size <= 3, "too many coordinates.");

	return FaceArea(*dom.subset_handler(), si, 0, aaPos);
}

/**
 * Calculates the area sum of faces in a given domain. Please note the sum
 * is returned for all faces selected by the selector sel in the domain dom.
 *
 * \param dom domain
 * \param sel selector
 *
 * \return \c number area sum
 *
 */
template <typename TDomain>
static number FaceArea(TDomain& dom, ISelector& sel)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	UG_ASSERT(TDomain::position_type::Size <= 3, "too many coordinates.");

	return FaceArea(sel, aaPos);
}


template <class TDomain, class TElem>
static TElem* GetElementByCoordinate(TDomain& dom, number x, number y, number z)
{
	vector3 tv(x, y, z);
	typename TDomain::position_type v;
	VecCopy(v, tv, 0);

	typename TDomain::grid_type& g = *dom.grid();
	return FindClosestByCoordinate<TElem>(v, g.template begin<TElem>(), g.template end<TElem>(), dom.position_accessor());
}


template <typename TDomain>
static number GetMaxEdgeLength(TDomain& dom)
{
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	typename TDomain::grid_type& g = *dom.grid();

	number maxLenSq = 0;
	for(EdgeIterator eiter = g.template begin<Edge>();
		eiter != g.template end<Edge>(); ++eiter)
	{
		maxLenSq = max(maxLenSq, EdgeLengthSq(*eiter, aaPos));
	}

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator com;
	number gMaxLenSq = com.allreduce(maxLenSq, PCL_RO_MAX);
	return sqrt(gMaxLenSq);
#else
	return sqrt(maxLenSq);
#endif
}

// end group domain_bridge
/// \}

namespace bridge{
namespace Domain{

/// \addtogroup domain_bridge
/// \{

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
			.add_method("empty", &TDomain::empty)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "Domain", tag);
	}


// 	MaxElementDiameter
	reg.add_function("MaxElementDiameter", static_cast<number (*)(TDomain&, int)>(
					 &MaxElementDiameter<TDomain>), grp);

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
					"", "Domain # Filename|save-dialog| endings=[\"ugx\"]",
					"Saves a domain", "No help");

//	SavePartitionMap
	reg.add_function("SavePartitionMap", &SavePartitionMap<TDomain>, grp,
					"Success", "PartitionMap # Domain # Filename|save-dialog",
					"Saves a partition map", "No help");

//	Domain Distribution
// note that an overload of this method exists, registered in disc_bridges/domain_disc_bridge
	reg.add_function("DistributeDomain",
		static_cast<bool (*)(TDomain&, PartitionMap&, bool)>(&DistributeDomain<TDomain>),
		grp);

//	PartitionDomain
	reg.add_function("PartitionDomain_MetisKWay",
					 static_cast<bool (*)(TDomain&, PartitionMap&, int, size_t, int, int)>(&PartitionDomain_MetisKWay<TDomain>), grp);
	reg.add_function("PartitionDomain_MetisKWay",
					 static_cast<bool (*)(TDomain&, PartitionMap&, int, size_t, SmartPtr<PartitionWeighting>)>(&PartitionDomain_MetisKWay<TDomain>), grp);

	reg.add_function("PartitionDomain_LevelBased",
					 &PartitionDomain_LevelBased<TDomain>, grp, "", "domain#partitionMap#numPartitions#level");

	reg.add_function("PartitionDistributedDomain_LevelBased",
					 &PartitionDistributedDomain_LevelBased<TDomain>, grp, "", "domain#partitionMap#numPartitions#level");

//	transform the domain
	reg.add_function("ScaleDomain", &ScaleDomain<TDomain>, grp, "", "dom#sx#sy#sz");
	reg.add_function("TranslateDomain", &TranslateDomain<TDomain>, grp, "", "dom#tx#ty#tz");

//  calculate area covered by faces
	reg.add_function("FaceArea", static_cast<number (*)(TDomain&, ISubsetHandler&, int, size_t)>(&FaceArea<TDomain>), grp, "Area sum", "Domain#Subset handler#Subset index#Grid level");
	reg.add_function("FaceArea", static_cast<number (*)(TDomain&, int, size_t)>(&FaceArea<TDomain>), grp, "Area sum", "Domaim#Subset index#Grid level");
	reg.add_function("FaceArea", static_cast<number (*)(TDomain&, int)>(&FaceArea<TDomain>), grp, "Area sum", "Domain#Subset index");
	reg.add_function("FaceArea", static_cast<number (*)(TDomain&, ISelector&)>(&FaceArea<TDomain>), grp, "Area sum", "Domain#Selector");

//	element access
	reg.add_function("GetVertexByCoordinate", &GetElementByCoordinate<TDomain, Vertex>, grp);
	reg.add_function("GetEdgeByCoordinate", &GetElementByCoordinate<TDomain, Edge>, grp);
	reg.add_function("GetFaceByCoordinate", &GetElementByCoordinate<TDomain, Face>, grp);
	reg.add_function("GetVolumeByCoordinate", &GetElementByCoordinate<TDomain, Volume>, grp);

//	geometry information
	reg.add_function("GetMaxEdgeLength", &GetMaxEdgeLength<TDomain>, grp);

//	debugging
	reg.add_function("TestDomainInterfaces", &TestDomainInterfaces<TDomain>, grp);

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
//	DomainInfo
	{
		reg.add_class_<DomainInfo>("DomainInfo", grp)
			.add_constructor()
			.add_method("element_type", &DomainInfo::element_type)
			.add_method("num_levels", &DomainInfo::num_levels)
			.add_method("num_elements_on_level", &DomainInfo::num_elements_on_level)
			.add_method("num_local_elements_on_level", &DomainInfo::num_local_elements_on_level)
			.add_method("num_local_ghosts_on_level", &DomainInfo::num_local_ghosts_on_level)
			.add_method("num_subsets", &DomainInfo::num_subsets)
			.add_method("subset_dim", &DomainInfo::subset_dim)
			.add_method("to_string", &DomainInfo::to_string)
			.set_construct_as_smart_pointer(true);
	}

//	IDomain
	{
		typedef IDomain<> T;
		reg.add_class_<T>("IDomain", grp)
			.add_method("domain_info", &T::domain_info, "DomainInfo")
			.add_method("get_dim", static_cast<int (T::*)() const>(&T::get_dim))
			.add_method("grid", static_cast<SmartPtr<MultiGrid> (T::*)()>(&T::grid), "grid")
			.add_method("subset_handler", static_cast<SmartPtr<MGSubsetHandler> (T::*)()>(&T::subset_handler))
			.add_method("create_additional_subset_handler", &T::create_additional_subset_handler, "bool")
			.add_method("additional_subset_handler_names", &T::additional_subset_handler_names, "vector<string>")
			.add_method("additional_subset_handler",
					static_cast<SmartPtr<MGSubsetHandler> (T::*)(string)>(&T::additional_subset_handler),
					"SubsetHandler")
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

// end group domain_bridge
/// \}

}// end Domain

/// \addtogroup domain_bridge
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
		RegisterCommon<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Domain::Functionality2d3d, CompileDomain2d3dList>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace
}// end of namespace
