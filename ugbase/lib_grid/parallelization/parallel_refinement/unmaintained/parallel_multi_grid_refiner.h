#ifndef __H__LIB_GRID__PARALLEL_MULTI_GRID_REFINER__
#define __H__LIB_GRID__PARALLEL_MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/multi_grid_refiner.h"
#include "../distributed_grid.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_refinement
/// @{

//	predeclarations
template <class TLayout>
class RefinementMarkDistributor;

///	DEPRECIATED
/**	DEPRECIATED!
 * This class is only intended for some test. Shouldn't be used.
 */
class ParallelMultiGridRefiner : public MultiGridRefiner
{
	friend class RefinementMarkDistributor<VertexLayout>;
	friend class RefinementMarkDistributor<EdgeLayout>;
	friend class RefinementMarkDistributor<FaceLayout>;
	friend class RefinementMarkDistributor<VolumeLayout>;
	
	public:
		//ParallelMultiGridRefiner();
		ParallelMultiGridRefiner(DistributedGridManager& distGridMgr);
		virtual ~ParallelMultiGridRefiner();

	///	BE CAREFUL! only for debugging purposes. Enabled by default.
	/**	Disable communication only if you now what you are doing.
	 *	Disabled communication can lead to asynchronous pcl-interfaces
	 *	and thus to severe errors during later communication.*/
		inline void enable_communication(bool bEnable)		{m_bCommunicationEnabled = bEnable;}
	///	only for debugging purposes
		inline bool communication_is_enabled()				{return m_bCommunicationEnabled;}
		
	protected:
		virtual void collect_objects_for_refine();

		virtual void refinement_step_begins();
		virtual void refinement_step_ends();
		
		virtual void set_rule(Vertex* e, RefinementMark mark);
		virtual void set_rule(Edge* e, RefinementMark mark);
		virtual void set_rule(Face* e, RefinementMark mark);
		virtual void set_rule(Volume* e, RefinementMark mark);

	/**	Distributes marks for all elements that are stored in
	 *	m_vNewlyMarkedInterface...
	 *	you may optionally pass a vector to which elements will
	 *	be written that received a mark during communication.
	 *	Elements will be appended to the existing ones.*/
		template <class TDistributor, class TCommunicator>
		void
		exchange_data(TDistributor& distributor,
						TCommunicator& communicator,
						std::vector<typename TDistributor::Element>* pvReceivedElemsOut = NULL);
					
		template <class TMarkDistributor>
		void mark_received_elements(TMarkDistributor& distributor);
		
		void clear_newly_marked_element_buffers();
		
	///	adjust selection based on received elements
		void adjust_parallel_selection(const std::vector<Vertex*>* pvVrts,
										const std::vector<Edge*>* pvEdges,
										const std::vector<Face*>* pvFaces,
										const std::vector<Volume*>* pvVolumes);
		
	protected:
		DistributedGridManager& m_distGridMgr;
		
		std::vector<Vertex*>	m_vNewlyMarkedInterfaceVrts;
		std::vector<Edge*>		m_vNewlyMarkedInterfaceEdges;
		std::vector<Face*>			m_vNewlyMarkedInterfaceFaces;
		std::vector<Volume*>		m_vNewlyMarkedInterfaceVols;
		
		bool m_bCommunicationEnabled;	///< only for debugging purposes
};

/// @}
}//	end of namespace

#endif
