#ifndef __H__UG__std_hnode_adjuster__
#define __H__UG__std_hnode_adjuster__

#include "../ref_mark_adjuster_interface.h"

namespace ug{

class StdHNodeAdjuster;
typedef SmartPtr<StdHNodeAdjuster> SPStdHNodeAdjuster;

///	Makes sure that elements are marked correctly so that hnode-refinement produces a valid grid.
/**	This adjuster regards the grid as a serial grid. If the grid represents a part
 * of a distributed grid, then the additional use of a parallel adjuster is required.*/
class StdHNodeAdjuster : public IRefMarkAdjuster
{
	public:
		static SPStdHNodeAdjuster create()	{return SPStdHNodeAdjuster(new StdHNodeAdjuster);}

		virtual ~StdHNodeAdjuster()	{}

		virtual void ref_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols);
};

}// end of namespace

#endif
