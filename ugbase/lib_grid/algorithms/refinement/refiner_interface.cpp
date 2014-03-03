// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 02.03.2012 (m,d,y)
 
#include "refiner_interface.h"
#include "lib_grid/lib_grid_messages.h"
#include "common/catch_std.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif

namespace ug{

bool IRefiner::mark(GridObject* o, RefinementMark refMark)
{
	switch(o->base_object_id()){
		case VERTEX:	return mark(static_cast<Vertex*>(o), refMark);
		case EDGE:		return mark(static_cast<Edge*>(o), refMark);
		case FACE:		return mark(static_cast<Face*>(o), refMark);
		case VOLUME:	return mark(static_cast<Volume*>(o), refMark);
	}
	return false;
}

RefinementMark IRefiner::get_mark(GridObject* o)
{
	switch(o->base_object_id()){
		case VERTEX:	return get_mark(static_cast<Vertex*>(o));
		case EDGE:		return get_mark(static_cast<Edge*>(o));
		case FACE:		return get_mark(static_cast<Face*>(o));
		case VOLUME:	return get_mark(static_cast<Volume*>(o));
	}
	return RM_NONE;
}

void IRefiner::adaption_begins()
{
	if(!m_messageHub.valid()){
		UG_THROW("A message-hub has to be assigned to IRefiner before adaption_begins may be called. "
				"Make sure that you assigned a grid to the refiner you're using.");
	}

//	we'll schedule an adaption-begins message
	if(adaptivity_supported())
		m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_ADAPTION_BEGINS));
	else
		m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_ADAPTION_BEGINS));

	m_adaptionIsActive = true;
}

void IRefiner::adaption_ends()
{
	if(!m_messageHub.valid()){
		UG_THROW("A message-hub has to be assigned to IRefiner before adaption_ends may be called. "
				"Make sure that you assigned a grid to the refiner you're using.");
	}

	if(adaptivity_supported())
		m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_ADAPTION_ENDS));
	else
		m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_ADAPTION_ENDS));

	m_adaptionIsActive = false;
}

void IRefiner::refine()
{
	try
	{
		if(!m_messageHub.valid()){
			UG_THROW("A message-hub has to be assigned to IRefiner before refine may be called. "
					"Make sure that you assigned a grid to the refiner you're using.");
		}

	//	we'll schedule an adaption-begins message, if adaption is not yet enabled.
		bool locallyActivatedAdaption = false;
		if(!m_adaptionIsActive){
			locallyActivatedAdaption = true;
			adaption_begins();
		}

	//	now post a message, which informs that refinement begins
//		if(adaptivity_supported())
//			m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_REFINEMENT_BEGINS));
//		else
//			m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_REFINEMENT_BEGINS));

	//	now perform refinement
		perform_refinement();

	//	post a message that refinement has been finished
//		if(adaptivity_supported())
//			m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_REFINEMENT_ENDS));
//		else
//			m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_REFINEMENT_ENDS));

	//	and finally - if we posted an adaption-begins message, then we'll post
	//	an adaption ends message, too.
		if(locallyActivatedAdaption)
			adaption_ends();
	}
	CATCH_STD_EXCEPTIONS();
}


bool IRefiner::coarsen()
{
//	if coarsen isn't supported, we'll leave right away
	if(!coarsening_supported())
		return false;

	if(!m_messageHub.valid()){
		UG_THROW("A message-hub has to be assigned to IRefiner before coarsen may be called. "
				"Make sure that you assigned a grid to the refiner you're using.");
	}

//	we'll schedule an adaption-begins message, if adaption is not yet enabled.
	bool locallyActivatedAdaption = false;
	if(!m_adaptionIsActive){
		locallyActivatedAdaption = true;
		adaption_begins();
	}

//	now post a message, which informs that coarsening begins
//	if(adaptivity_supported())
//		m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_COARSENING_BEGINS));
//	else
//		m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_COARSENING_BEGINS));

//	now perform coarsening
	bool retVal = perform_coarsening();

//	post a message that coarsening has been finished
//	if(adaptivity_supported())
//		m_messageHub->post_message(GridMessage_Adaption(GMAT_HNODE_COARSENING_ENDS));
//	else
//		m_messageHub->post_message(GridMessage_Adaption(GMAT_GLOBAL_COARSENING_ENDS));

//	and finally - if we posted an adaption-begins message, then we'll post
//	an adaption ends message, too.
	if(locallyActivatedAdaption)
		adaption_ends();

//	done
	return retVal;
}


void IRefiner::set_message_hub(SPMessageHub msgHub)
{
	m_messageHub = msgHub;
}

void IRefiner::set_adjusted_marks_debug_filename(const char* filename)
{
	if(!filename)
		m_adjustedMarksDebugFilename = "";
	else
		m_adjustedMarksDebugFilename = filename;
}

void IRefiner::num_marked_edges(std::vector<int>& numMarkedEdgesOut)
{
	#ifdef UG_PARALLEL
		std::vector<int> numLocal;
		num_marked_edges_local(numLocal);
		pcl::ProcessCommunicator com;
		com.allreduce(numLocal, numMarkedEdgesOut, PCL_RO_SUM);
	#else
		num_marked_edges_local(numMarkedEdgesOut);
	#endif
}

void IRefiner::num_marked_faces(std::vector<int>& numMarkedFacesOut)
{
	#ifdef UG_PARALLEL
		std::vector<int> numLocal;
		num_marked_faces_local(numLocal);
		pcl::ProcessCommunicator com;
		com.allreduce(numLocal, numMarkedFacesOut, PCL_RO_SUM);
	#else
		num_marked_faces_local(numMarkedFacesOut);
	#endif
}

void IRefiner::num_marked_volumes(std::vector<int>& numMarkedVolsOut)
{
	#ifdef UG_PARALLEL
		std::vector<int> numLocal;
		num_marked_volumes_local(numLocal);
		pcl::ProcessCommunicator com;
		com.allreduce(numLocal, numMarkedVolsOut, PCL_RO_SUM);
	#else
		num_marked_volumes_local(numMarkedVolsOut);
	#endif
}

}// end of namespace
