// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.08.2012 (m,d,y)

#include "check_associated_elements.h"
#include "lib_grid/grid/grid_util.h"

namespace ug{
namespace grid_unit_tests{

void CheckAssociatedEdgesOfVolumes(Grid& g)
{
//	VOLOPT_STORE_ASSOCIATED_EDGES has to be enabled, so that this method makes sense...
	if(!g.option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES)){
		UG_LOG("WARNING: Autoenabling VOLOPT_STORE_ASSOCIATED_EDGES in CheckAssociatedEdgesOfVolumes.\n");
		g.enable_options(VOLOPT_STORE_ASSOCIATED_EDGES);
	}

	g.begin_marking();

//	iterate over all volumes
	for(VolumeIterator iter = g.volumes_begin(); iter != g.volumes_end(); ++iter){
		Volume* vol = *iter;
		g.clear_marks();

	//	get all associated edges
		std::vector<EdgeBase*> edges;
		CollectAssociated(edges, g, vol);

	//	iterate over them and mark them. Make sure none is marked twice.
		for(size_t i = 0; i < edges.size(); ++i){
			if(g.is_marked(edges[i])){
				UG_THROW("Edge is contained in associated edges of volume twice!");
			}
			g.mark(edges[i]);
		}

	//	make sure that the elements have the right order, if VOLOPT_AUTOGENERATE_EDGES
	//	or GRIDOPT_AUTOGENERATE_SIDES is active
		if(g.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES) ||
			g.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
		{
		//	edges should be sorted...
		//	also check get_edge(vol, i)...
			EdgeDescriptor ed;
			for(size_t i_ed = 0; i_ed < vol->num_edges(); ++i_ed){
				vol->edge_desc(i_ed, ed);
				EdgeBase* e = g.get_edge(ed);

				if(e != g.get_edge(vol, i_ed)){
					UG_THROW("Grid::get_edge(vol, i) does not return the i-th edge of vol!");
				}

				if(e != edges[i_ed]){
					UG_THROW("AssociatedEdges should contain edges in the correct order when autogeneration is active!");
				}
			}
		}

	//	now iterate over associated edges of all corner vertices of vol
	//	and make sure that each is marked
		Volume::ConstVertexArray vrts = vol->vertices();
		for(size_t i_vrt = 0; i_vrt < vol->num_vertices(); ++i_vrt){
			CollectAssociated(edges, g, vrts[i_vrt]);
			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				if(VolumeContains(vol, edges[i_edge])){
					if(!g.is_marked(edges[i_edge])){
						UG_THROW("Edge is contained in volume but not in volume's associated-edge-container!");
					}
				}
			}
		}
	}

	g.end_marking();
}

void CheckAssociatedVolumesOfEdges(Grid& g)
{
//	VOLOPT_STORE_ASSOCIATED_EDGES has to be enabled, so that this method makes sense...
	if(!g.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES)){
		UG_LOG("WARNING: Autoenabling EDGEOPT_STORE_ASSOCIATED_VOLUMES in CheckAssociatedVolumesOfEdges.\n");
		g.enable_options(EDGEOPT_STORE_ASSOCIATED_VOLUMES);
	}

	g.begin_marking();

//	iterate over all volumes
	for(EdgeBaseIterator iter = g.edges_begin(); iter != g.edges_end(); ++iter){
		EdgeBase* e = *iter;
		g.clear_marks();

	//	get all associated volumes
		std::vector<Volume*> vols;
		CollectAssociated(vols, g, e);

	//	iterate over them and mark them. Make sure none is marked twice.
		for(size_t i = 0; i < vols.size(); ++i){
			if(g.is_marked(vols[i])){
				UG_THROW("Volume is contained in associated volumes of edge twice!");
			}
			g.mark(vols[i]);
		}

	//	now iterate over associated volumes of all corner vertices of e
	//	and make sure that each is marked
		EdgeBase::ConstVertexArray vrts = e->vertices();
		for(size_t i_vrt = 0; i_vrt < e->num_vertices(); ++i_vrt){
			CollectAssociated(vols, g, vrts[i_vrt]);
			for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
				if(VolumeContains(vols[i_vol], e)){
					if(!g.is_marked(vols[i_vol])){
						UG_THROW("Volume is contained in edge but not in edge's associated-volume-container!");
					}
				}
			}
		}
	}

	g.end_marking();
}

}//	end of namespace
}//	end of namespace
