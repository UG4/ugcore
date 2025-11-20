/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <vector>
#include "registry/registry.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_grid/tools/subset_group.h"

#include <lib_grid/algorithms/projections/overlying_subset_finder.hpp>

using namespace std;

namespace ug{

/// TEMPORARY QUICK HACK! DON'T USE! WILL BE REPLACED SOON!
class DomainRayTracer {
	public:
		using vector_t = vector3;

		DomainRayTracer(Domain3d& dom) :
			m_tree(*dom.grid(), dom.position_attachment()),
			m_small(SMALL),
			m_dom(&dom)
		{}

		void set_small(number small)	{m_small = small;}

		void init (const std::vector<int>& subsetIndices)
		{
			MultiGrid& mg = *m_dom->grid();
			MGSubsetHandler& sh = *m_dom->subset_handler();

			std::vector<Triangle*>	tris;
			for(size_t is = 0; is < subsetIndices.size(); ++is){
				int si = subsetIndices[is];
				for(int lvl = 0; lvl < (int)sh.num_levels(); ++lvl){
					for(TriangleIterator it = sh.begin<Triangle>(si, lvl);
						it != sh.end<Triangle>(si, lvl); ++it)
					{
						Triangle* t = *it;
						if(!mg.has_children(t))
							tris.push_back(t);
					}
				}
			}

			m_tree.create_tree(tris.begin(), tris.end());
		}
		
		void init (const char* subsetNames)
		{
			SubsetGroup ssGrp(m_dom->subset_handler());
			ssGrp.add(TokenizeString(subsetNames));
			init(ssGrp.index_vector());
		}

		size_t trace_ray(const vector_t& from, const vector_t& dir)
		{
			m_tracePoints.clear();
			RayElementIntersections(m_intersectionRecords, m_tree, from, dir, m_small);
			for(size_t _vfeI = 0; _vfeI < m_intersectionRecords.size(); ++_vfeI){ intersection_record_t& r = m_intersectionRecords[_vfeI];{
				m_tracePoints.push_back(PointOnRay(from, dir, r.smin));
			}};
			return m_tracePoints.size();
		}

		size_t num_trace_points() const				{return m_tracePoints.size();}
		const vector_t& trace_point(size_t i) const	{return m_tracePoints[i];}
		number trace_point_x(size_t i) const		{return m_tracePoints[i].x();}
		number trace_point_y(size_t i) const		{return m_tracePoints[i].y();}
		number trace_point_z(size_t i) const		{return m_tracePoints[i].z();}

		// size_t closest_point_index() const
		// {
		// 	if(m_tracePoints.empty())
		// 		return 0;

		// //	todo: use local-ray-coord
		// 	size_t index = 0;
		// 	number minDistSq = VecDistanceSq(m_from, m_tracePoints[0]);
		// 	for(size_t i = 1; i < m_tracePoints.size(); ++i){
		// 		number distSq = VecDistanceSq(m_from, m_tracePoints[i]);
		// 		if(distSq < minDistSq){
		// 			minDistSq = distSq;
		// 			index = i;
		// 		}
		// 	}

		// 	return index;
		// }

	///	local coordinate regarding trace-ray
		// number local_trace_point_coordinate(size_t i) const;

		// int trace_point_subset_index(size_t i) const;

	private:
		using tree_t = lg_ntree<3, 3, Triangle>;
		using intersection_record_t = RayElemIntersectionRecord<Triangle*>;

		std::vector<vector_t>	m_tracePoints;
		std::vector<intersection_record_t>	m_intersectionRecords;
		tree_t					m_tree;
		number					m_small;
		Domain3d*				m_dom;
};



namespace bridge {
namespace domain_ray_tracing {

struct Functionality {

	static void Common(Registry& reg, string grp)
	{
		using T = DomainRayTracer;
		
		reg.add_class_<DomainRayTracer>("DomainRayTracer", grp)
				.add_constructor<void (*)(Domain3d&)> ()
				.add_method("set_small", &T::set_small, "", "small", "")
				.add_method("init", static_cast<void (T::*) (const std::vector<int>&)>(&T::init), "", "subsetIndices", "")
				.add_method("init", static_cast<void (T::*) (const char*)>(&T::init), "", "subsetNames", "")
				.add_method("trace_ray", &T::trace_ray, "", "rayFrom # rayTo", "")
				.add_method("num_trace_points", &T::num_trace_points, "numPoints", "", "")
				.add_method("trace_point", &T::trace_point, "point", "index", "")
				.add_method("trace_point_x", &T::trace_point_x, "xCoord", "index", "")
				.add_method("trace_point_y", &T::trace_point_y, "yCoord", "index", "")
				.add_method("trace_point_z", &T::trace_point_z, "zCoord", "index", "");
	}



	/**
	 * Function called for the registration of Domain dependent parts.
	 * All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param grp				group for sorting of functionality
	 */
	template <typename TDomain>
	static void Domain(Registry& reg, string grp)
	{
		using domain_type = TDomain;

		using T = OverlyingSubsetFinder<TDomain>;

		string suffix = GetDomainSuffix<TDomain>();
		string tag = GetDomainTag<TDomain>();


		string name = string("OverlyingSubsetFinder").append(suffix);
		
		reg.add_class_<T>(name, grp)
				.template add_constructor<void (*)(SmartPtr<domain_type>, const std::string& subsets)> ()
				.add_method("findOverlyingSubset", &T::findOverlyingSubset, "", "point to search over", "")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "OverlyingSubsetFinder", tag);
						
	}

}; // end Functionality

}  // end of namespace domain_ray_tracing

void RegisterBridge_DomainRayTracing(Registry& reg, string grp) {
	grp.append("/RayTracing");

	try {
		RegisterDomain2d3dDependent<domain_ray_tracing::Functionality>(reg, grp);
#ifdef UG_DIM_3
		RegisterCommon<domain_ray_tracing::Functionality>(reg, grp);
#endif
	} UG_REGISTRY_CATCH_THROW(grp);
}

} // end of namespace bridge
}//	end of namespace ug
