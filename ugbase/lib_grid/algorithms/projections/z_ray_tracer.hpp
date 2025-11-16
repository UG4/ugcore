#include <vector>
#include <string>
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include <lib_grid/grid_objects/grid_dim_traits.h>
#ifdef UG_PARALLEL
#include "lib_grid/parallelization/gather_grid.h"
#include "lib_grid/parallelization/distributed_grid.h"
#endif
#include "common/common.h"

namespace ug
{

template <typename TDomain>
class ZRayTracer
{

    ///	grid function type
    using domain_type = TDomain;
			
///	world dimension
	static constexpr int dim = domain_type::dim;
	
/// type of the position accessor
    using position_attachment_type = typename domain_type::position_attachment_type;
	
///	side type
    using side_t = typename grid_dim_traits<dim-1>::element_type;

    using top_tracer_tree_t = lg_ntree<dim-1, dim, side_t>;

    using top_intersection_record_t = RayElemIntersectionRecord<side_t*>;
    
public:

///	class constructor
    ZRayTracer
    (
        SmartPtr<domain_type> sp_domain, ///< the domain
        const std::string & top_ss_names, ///< top surface subset names
        int grid_level,
        bool useLocalTopFacesOnly = false
    )
    :	m_sp_domain (sp_domain),
        m_top_ss_grp (sp_domain->subset_handler (), TokenizeString (top_ss_names))
    {
        #ifndef UG_PARALLEL
            useLocalTopFacesOnly = true;
        #endif

        using SideIterator = typename Grid::traits<side_t>::iterator;

        MultiGrid & mg = * m_sp_domain->grid ();
        MGSubsetHandler & sh = * m_sp_domain->subset_handler ();
        std::vector<side_t*> topSides;
        
    //	get all the elements to collect
        if(useLocalTopFacesOnly){
            m_top_tracer_tree.set_grid(*m_sp_domain->grid (), m_sp_domain->position_attachment ());
            for (size_t i = 0; i < m_top_ss_grp.size (); i++)
            {
                int si = m_top_ss_grp [i];
                
                if (grid_level >= 0){ // if the grid level for the top is specified
                    for (SideIterator it = sh.begin<side_t> (si, grid_level);
                                            it != sh.end<side_t> (si, grid_level); ++it)
                        topSides.push_back (*it);
                }
                else{
                    for (int lvl = 0; lvl < (int) sh.num_levels(); lvl++){
                        for (SideIterator it = sh.begin<side_t> (si, lvl);
                                                it != sh.end<side_t> (si, lvl); ++it)
                        {
                            side_t* t = *it;
                            if (! mg.has_children (t))
                                topSides.push_back (t);
                        }
                    }
                }
            }
        }
        else {
            #ifdef UG_PARALLEL
                m_top_tracer_tree.set_grid(m_top_grid, m_sp_domain->position_attachment ());
                
                DistributedGridManager* dgm = mg.distributed_grid_manager();

            //	select the Sides on the top
                Selector sel (mg);
                for (size_t i = 0; i < m_top_ss_grp.size (); i++)
                {
                    int si = m_top_ss_grp [i];
                    
                    if (grid_level >= 0) // if the grid level for the top is specified
                        for (SideIterator it = sh.begin<side_t> (si, grid_level);
                                                it != sh.end<side_t> (si, grid_level); ++it)
                            sel.select (*it);
                    else
                        for (int lvl = 0; lvl < (int) sh.num_levels(); lvl++)
                            for (SideIterator it = sh.begin<side_t> (si, lvl);
                                                    it != sh.end<side_t> (si, lvl); ++it)
                            {
                                side_t* t = *it;
                                if (! (mg.has_children (t) || (dgm && dgm->is_ghost(t))))
                                    sel.select (t);
                            }
                }
                
            //	copy the top Sides into a new grid
                GridDataSerializationHandler serializer;
                serializer.add
                    (GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (mg, m_sp_domain->position_attachment ()));
                    
                GridDataSerializationHandler deserializer;
                deserializer.add
                    (GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (m_top_grid, m_sp_domain->position_attachment ()));

                AllGatherGrid (m_top_grid, sel, serializer, deserializer);

                for (SideIterator it = m_top_grid.begin<side_t> (); it != m_top_grid.end<side_t> (); ++it)
                    topSides.push_back (*it);
                
                // UG_LOG("DEBUG: SAVING allgathered m_top_grid to file in ls_init...\n");
                // SaveGridToFile(m_top_grid, mkstr("top_grid_p" << pcl::ProcRank() << ".ugx").c_str(),
                //                m_sp_domain->position_attachment());
            #endif
        }
    //	compose the tree
        m_top_tracer_tree.create_tree (topSides.begin (), topSides.end ());
    }
    
///	computes the minimum z-coordiate of the top (returns true if the top found, false otherwise)
    bool get_min_at
    (
        const MathVector<dim> & over, ///< to look over this point
        number tolerance, ///< the tolerance of the ray tracer
        number & z ///< the result
    )
    {
    //	find all the intersections
        MathVector<dim> up_dir;
        up_dir = 0;
        up_dir [dim - 1] = 1;
        m_top_intersection_records.clear ();
        RayElementIntersections (m_top_intersection_records, m_top_tracer_tree, over, up_dir, tolerance);
        
    //	check if there are intersections at all
        if (m_top_intersection_records.size () == 0)
            return false;
        
    //	find the lowest point
        MathVector<dim> x = PointOnRay (over, up_dir, m_top_intersection_records[0].smin);
        number z_min = x [dim - 1];
        // UG_LOG("<dbg>  x: " << x << std::endl);
        for (size_t i = 1; i < m_top_intersection_records.size (); i++)
        {
            top_intersection_record_t & r = m_top_intersection_records [i];
            x = PointOnRay (over, up_dir, r.smin);
            if (x [dim - 1] < z_min)
                z_min = x [dim - 1];
        }
        
        z = z_min;
        return true;
    }

private:

///	multigrid of the domain
    SmartPtr<domain_type> m_sp_domain;

///	subset group of the top faces
    SubsetGroup m_top_ss_grp;

#	ifdef UG_PARALLEL
///	auxiliary grid of the top faces
    Grid m_top_grid;
#	endif
///	tracer tree of the top faces
    top_tracer_tree_t m_top_tracer_tree;
    
///	array to store all the intersections
    std::vector<top_intersection_record_t> m_top_intersection_records;
};

}
