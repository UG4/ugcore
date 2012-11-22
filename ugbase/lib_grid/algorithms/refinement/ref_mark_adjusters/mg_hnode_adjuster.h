// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#ifndef __H__UG__mg_hnode_adjuster__
#define __H__UG__mg_hnode_adjuster__

#include "../ref_mark_adjuster_interface.h"

namespace ug{

class MGHNodeAdjuster;
typedef SmartPtr<MGHNodeAdjuster> SPMGHNodeAdjuster;

///	Makes sure that that hierarchical dependencies are considered during hnode refinement.
/**	This adjuster regards the grid as a serial grid. If the grid represents a part
 * of a distributed grid, then the additional use of a parallel adjuster is required.*/
class MGHNodeAdjuster : public IRefMarkAdjuster
{
	public:
		static SPMGHNodeAdjuster create()		{return SPMGHNodeAdjuster(new MGHNodeAdjuster);}

		virtual ~MGHNodeAdjuster()	{}

		virtual AdjustRetVal ref_marks_changed(IRefiner& ref,
											const std::vector<VertexBase*>& vrts,
											const std::vector<EdgeBase*>& edges,
											const std::vector<Face*>& faces,
											const std::vector<Volume*>& vols);
};

}// end of namespace

#endif
