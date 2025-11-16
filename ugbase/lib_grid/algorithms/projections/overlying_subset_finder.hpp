#include <common/util/string_util.h>
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include <lib_grid/grid_objects/grid_dim_traits.h>
#ifdef UG_PARALLEL
#include "lib_grid/parallelization/gather_grid.h"
#include "lib_grid/parallelization/broadcast.h"
#include "lib_grid/parallelization/distributed_grid.h"
#endif
#include "lib_grid/file_io/file_io.h"

namespace ug {

	template <typename TDomain>
	class OverlyingSubsetFinder
	{
		static constexpr int dim = TDomain::dim;
		using domain_type = TDomain;

		using side_t = typename grid_dim_traits<dim-1>::element_type;

		using top_tracer_tree_t = lg_ntree<dim-1, dim, side_t>;

		using top_intersection_record_t = RayElemIntersectionRecord<side_t*>;

		using position_attachment_type = typename domain_type::position_attachment_type;
		using position_accessor_type = typename domain_type::position_accessor_type;

		using SideIterator = typename Grid::traits<side_t>::iterator;

		public:

			OverlyingSubsetFinder(        
				SmartPtr<TDomain> sp_domain, ///< the domain
        		const std::string & top_ss_names) ///< top surface subset names
				: m_sp_domain (sp_domain),
				m_top_ss_grp (sp_domain->subset_handler (), TokenizeString (top_ss_names))
			{

				MultiGrid & mg = * m_sp_domain->grid ();
        		MGSubsetHandler & sh = * m_sp_domain->subset_handler ();
        		std::vector<side_t*> topSides;			

				m_top_tracer_tree.set_grid(m_top_grid, m_sp_domain->position_attachment ());
                
				mg.attach_to<side_t>(subset_index_attachment);
				m_top_grid.attach_to<side_t>(subset_index_attachment);

				Grid::AttachmentAccessor<side_t, Attachment<int> > subset_index_accessor;
				subset_index_accessor.access(mg, subset_index_attachment);

            //	select the Sides on the top
                Selector sel (mg);
                for (size_t i = 0; i < m_top_ss_grp.size (); i++)
                {
                    int si = m_top_ss_grp [i];                   
                    
					for (SideIterator it = sh.begin<side_t> (si, 0); it != sh.end<side_t> (si, 0); ++it)
					{
						sel.select (*it);
						subset_index_accessor[*it] = si;          
					}          
                }

				position_accessor_type accessor;
				accessor.access(mg, m_sp_domain->position_attachment());


				//UG_LOG_ALL_PROCS("Total count of Edges on this node before distributing: " << sel.num<side_t>() << std::endl);

				// for (SideIterator it = sel.begin<side_t>(); it != sel.end<side_t>(); ++it)
				// {
				// 	UG_LOG_ALL_PROCS(subset_index_accessor[(*it)] << ": " << accessor[(*it)->vertex(0)] << "-" << accessor[(*it)->vertex(1)] << std::endl);
				// }
                
            	//	copy the top Sides into a new grid

				#ifdef UG_PARALLEL

                GridDataSerializationHandler serializer;
                serializer.add
                    (GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (mg, m_sp_domain->position_attachment ()));
                serializer.add
					(GeomObjAttachmentSerializer<side_t, Attachment<int> >::create (mg, subset_index_attachment));   
                
				GridDataSerializationHandler deserializer;
                deserializer.add
                    (GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (m_top_grid, m_sp_domain->position_attachment ()));
				deserializer.add
                    (GeomObjAttachmentSerializer<side_t, Attachment<int> >::create (m_top_grid, subset_index_attachment));

                BroadcastGrid (m_top_grid, sel, serializer, deserializer, 0);

				#endif

                for (SideIterator it = m_top_grid.begin<side_t> (); it != m_top_grid.end<side_t> (); ++it)
                    topSides.push_back (*it);

       			m_top_tracer_tree.create_tree (topSides.begin (), topSides.end ());

				// debug();

				// SaveGridToFile(m_top_grid, mkstr("top_grid_p" << pcl::ProcRank() << ".ugx").c_str(),
                //                m_sp_domain->position_attachment());
			}

			std::string findOverlyingSubset(const std::vector<number> & point)
			{

				if (point.size () != dim)
					UG_THROW ("OverlyingSubsetFinder::findOverlyingSubset: Exactly " << dim << " coordinates have to be specified!");
		
				MathVector<dim> measure_point;

				for(size_t d = 0; d < dim; d++)
				{
					measure_point[d] = point[d];
				}

				//	find all the intersections
				MathVector<dim> up_dir;
				up_dir = 0;
				up_dir [dim - 1] = 1;
				m_top_intersection_records.clear ();
				RayElementIntersections (m_top_intersection_records, m_top_tracer_tree, measure_point, up_dir, 0.001);
				
			//	check if there are intersections at all
				if (m_top_intersection_records.size () == 0)
					return "";
				
				Grid::AttachmentAccessor<side_t, Attachment<int> > subset_index_accessor;
				subset_index_accessor.access(m_top_grid, subset_index_attachment);

				//debug();

			//	find the lowest point
				MathVector<dim> x = PointOnRay (measure_point, up_dir, m_top_intersection_records[0].smin);
				number z_min = x [dim - 1];
				int si = subset_index_accessor[m_top_intersection_records[0].elem];				

				// UG_LOG("<dbg>  x: " << x << std::endl);
				for (size_t i = 1; i < m_top_intersection_records.size (); i++)
				{
					top_intersection_record_t & r = m_top_intersection_records [i];
					x = PointOnRay (measure_point, up_dir, r.smin);
					if (x [dim - 1] < z_min)
					{
						z_min = x [dim - 1];
						si = subset_index_accessor[r.elem];
					}
				}
				
				return m_sp_domain->subset_handler ()->get_subset_name(si);
			}

		void debug()
		{
			position_accessor_type accessor;
			accessor.access(m_top_grid, m_sp_domain->position_attachment());

			Grid::AttachmentAccessor<side_t, Attachment<int> > subset_index_accessor;
			subset_index_accessor.access(m_top_grid, subset_index_attachment);

			UG_LOG_ALL_PROCS("Total count of Edges: " << m_top_grid.num<side_t>() << std::endl);

			for (SideIterator it = m_top_grid.begin<side_t>(); it != m_top_grid.end<side_t>(); ++it)
			{
				UG_LOG_ALL_PROCS(subset_index_accessor[(*it)] << ": ");
				for (size_t i = 0; i < dim; i++)
				{
					UG_LOG_ALL_PROCS(accessor[(*it)->vertex(i)] << "-");
				}
				UG_LOG_ALL_PROCS(std::endl);
			}

		}

		private:
			SmartPtr<domain_type> m_sp_domain;

    		SubsetGroup m_top_ss_grp;
    		Grid m_top_grid;
			Attachment<int> subset_index_attachment;

			top_tracer_tree_t m_top_tracer_tree;
			
		///	array to store all the intersections
			std::vector<top_intersection_record_t> m_top_intersection_records;

	};

}
