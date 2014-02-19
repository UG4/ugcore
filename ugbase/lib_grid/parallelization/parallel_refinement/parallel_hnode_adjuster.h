// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#ifndef __H__UG__parallel_hnode_adjuster__
#define __H__UG__parallel_hnode_adjuster__

#include "lib_grid/algorithms/refinement/ref_mark_adjuster_interface.h"
#include "../distributed_grid.h"
#include "pcl/pcl_interface_communicator.h"
#include "pcl/pcl_process_communicator.h"

namespace ug{

class ParallelHNodeAdjuster;
typedef SmartPtr<ParallelHNodeAdjuster> SPParallelHNodeAdjuster;

///	Makes sure that that marks are propagated over process interfaces
class ParallelHNodeAdjuster : public IRefMarkAdjuster
{
	public:
		static SPParallelHNodeAdjuster create()		{return SPParallelHNodeAdjuster(new ParallelHNodeAdjuster);}

		virtual ~ParallelHNodeAdjuster()	{}

		virtual void ref_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols);

	private:
		pcl::ProcessCommunicator m_procCom;
		pcl::InterfaceCommunicator<VertexLayout> m_intfComVRT;
		pcl::InterfaceCommunicator<EdgeLayout> m_intfComEDGE;
		pcl::InterfaceCommunicator<FaceLayout> m_intfComFACE;
};

}// end of namespace

#endif
