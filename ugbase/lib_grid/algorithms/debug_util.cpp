// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 13.01.2011 (m,d,y)
 
#include "debug_util.h"
using namespace std;

namespace ug{

void PrintElementNumbers(const GeometricObjectCollection& goc)
{
	UG_LOG("grid element numbers:\n");
	for(size_t i = 0; i < goc.num_levels(); ++i)
	{
		if(goc.num_levels() > 1){
			UG_LOG("level " << i << endl);
		}
		UG_LOG("  vertices total:\t" << goc.num<VertexBase>(i) << endl);
		UG_LOG("    normal vrts:\t" << goc.num<Vertex>(i) << endl);
		UG_LOG("    hanging vrts:\t" << goc.num<HangingVertex>(i) << endl);

		UG_LOG("  edges total:\t\t" << goc.num<EdgeBase>(i) << endl);
		UG_LOG("    normal edges:\t" << goc.num<Edge>(i) << endl);
		UG_LOG("    constraining edges:\t" << goc.num<ConstrainingEdge>(i) << endl);
		UG_LOG("    constrained edges:\t" << goc.num<ConstrainedEdge>(i) << endl);

		UG_LOG("  faces total:\t\t" << goc.num<Face>(i) << endl);
		UG_LOG("    normal triangles:\t" << goc.num<Triangle>(i) << endl);
		UG_LOG("    constraining tris:\t" << goc.num<ConstrainingTriangle>(i) << endl);
		UG_LOG("    constrained tris:\t" << goc.num<ConstrainedTriangle>(i) << endl);

		UG_LOG("    normal quads:\t" << goc.num<Quadrilateral>(i) << endl);
		UG_LOG("    constraining quads:\t" << goc.num<ConstrainingQuadrilateral>(i) << endl);
		UG_LOG("    constrained quads:\t" << goc.num<ConstrainedQuadrilateral>(i) << endl);

		UG_LOG("  volumes total:\t" << goc.num<Volume>(i) << endl);

		UG_LOG(endl);
	}
}

void PrintGridElementNumbers(Grid& grid)
{
	PrintElementNumbers(grid.get_geometric_object_collection());
}

void PrintGridElementNumbers(MultiGrid& mg)
{
	PrintElementNumbers(mg.get_geometric_object_collection());
}

void PrintGridElementNumbers(GridSubsetHandler& sh)
{
	PrintElementNumbers(sh.get_geometric_object_collection());
}

template <class TGeomObj>
void PrintAttachmentInfo(Grid& grid)
{
	typedef Grid::AttachmentPipe::ConstAttachmentEntryIterator AttIter;

//	iterate over all attachments of the grid
	Grid::AttachmentPipe& pipe = grid.get_attachment_pipe<TGeomObj>();

	int counter = 1;
	size_t totalSize = 0;
	for(AttIter iter = pipe.attachments_begin();
		iter != pipe.attachments_end(); ++iter, ++counter)
	{
	//	name
		IAttachment* att = iter->m_pAttachment;
		UG_LOG("Attachment " << counter << " (" << att->get_name() << "): ");

	//	size
		IAttachmentDataContainer* con = iter->m_pContainer;
		UG_LOG(con->occupied_memory() << " bytes\n");
		totalSize += con->occupied_memory();
	}

	UG_LOG(counter - 1 << " attachments with a total size of "
			<< totalSize << " bytes.\n");
}

void PrintAttachmentInfo(Grid& grid)
{
	UG_LOG("Vertex Attachments:\n");
	PrintAttachmentInfo<VertexBase>(grid);

	UG_LOG("\nEdge Attachments:\n");
	PrintAttachmentInfo<EdgeBase>(grid);

	UG_LOG("\nFace Attachments:\n");
	PrintAttachmentInfo<Face>(grid);

	UG_LOG("\nVolume Attachments:\n");
	PrintAttachmentInfo<Volume>(grid);
}

}// end of namespace
