// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include <vector>
#include "registry/registry.h"
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_grid/tools/subset_group.h"


using namespace std;

namespace ug{

/// TEMPORARY QUICK HACK! DON'T USE! WILL BE REPLACED SOON!
class DomainRayTracer {
	public:
		typedef vector3 vector_t;

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
			for_each_in_vec(intersection_record_t& r, m_intersectionRecords){
				m_tracePoints.push_back(PointOnRay(from, dir, r.smin));
			}end_for;
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
		typedef lg_ntree<3, 3, Triangle>	tree_t;
		typedef RayElemIntersectionRecord<Triangle*>	intersection_record_t;

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
		typedef DomainRayTracer T;
		
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
}; // end Functionality

}  // end of namespace domain_ray_tracing

void RegisterBridge_DomainRayTracing(Registry& reg, string grp) {
	grp.append("/RayTracing");

	try {
#ifdef UG_DIM_3
		RegisterCommon<domain_ray_tracing::Functionality>(reg, grp);
#endif
	} UG_REGISTRY_CATCH_THROW(grp);
}

} // end of namespace bridge
}//	end of namespace ug
