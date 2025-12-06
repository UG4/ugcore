/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "refiner_interface.h"
#include "lib_grid/lib_grid_messages.h"
#include "common/catch_std.h"

//	FOR DEBUGGING ONLY:
#include "lib_grid/tools/periodic_boundary_manager.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
	#include "pcl/pcl_util.h"
#endif

namespace ug{

bool IRefiner::mark(GridObject* o, RefinementMark refMark)
{
	switch(o->base_object_id()){
		case VERTEX: return mark(static_cast<Vertex*>(o), refMark);
		case EDGE: return mark(static_cast<Edge*>(o), refMark);
		case FACE: return mark(static_cast<Face*>(o), refMark);
		case VOLUME: return mark(static_cast<Volume*>(o), refMark);
	}
	return false;
}

RefinementMark IRefiner::get_mark(GridObject* o) const
{
	switch(o->base_object_id()){
		case VERTEX: return get_mark(static_cast<Vertex*>(o));
		case EDGE: return get_mark(static_cast<Edge*>(o));
		case FACE: return get_mark(static_cast<Face*>(o));
		case VOLUME: return get_mark(static_cast<Volume*>(o));
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
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
	PROFILE_BEGIN_GROUP(IRefiner_refine, "grid");
	try
	{
		if(m_projector.invalid() && grid()){
			Grid& g = *grid();
			if(g.has_vertex_attachment(aPosition)){
				m_projector = make_sp(new RefinementProjector(MakeGeometry3d(g, aPosition)));
			}
			else if(g.has_vertex_attachment(aPosition2)){
				m_projector = make_sp(new RefinementProjector(MakeGeometry3d(g, aPosition2)));
			}
			else if(g.has_vertex_attachment(aPosition1)){
				m_projector = make_sp(new RefinementProjector(MakeGeometry3d(g, aPosition1)));
			}
		}

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
		if(locallyActivatedAdaption){
			adaption_ends();
		}
	}
	CATCH_STD_EXCEPTIONS();
	
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
}


bool IRefiner::coarsen()
{
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
	PROFILE_BEGIN_GROUP(IRefiner_coarsen, "grid");
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
	#ifdef UG_PARALLEL
		PCL_DEBUG_BARRIER_ALL();
	#endif
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


int IRefiner::
get_local_edge_mark(Face* f, Edge* e) const
{
	const int edgeInd = GetEdgeIndex(f, e);
	UG_COND_THROW(edgeInd == -1, "Given edge is not an edge of the given face.");

	if(marked_local(f)){
		const int faceLocalMark = get_local_mark(f);
		return (faceLocalMark >> edgeInd) & 1;
	}
	else if(marked_full(f)){
		return 1;
	}
	else if(marked_closure(f)){

		return static_cast<int>(marked_full(e));
	}
	return 0;
}

int IRefiner::
get_local_edge_mark(Volume* vol, Edge* e) const
{
	const int edgeInd = GetEdgeIndex(vol, e);
	UG_COND_THROW(edgeInd == -1, "Given edge is not an edge of the given volume.");

	if(marked_local(vol)){
		const int volLocalMark = get_local_mark(vol);
		return (volLocalMark >> edgeInd) & 1;
	}
	else if(marked_full(vol)){
		return 1;
	}
	else if(marked_closure(vol)){
		return static_cast<int>(marked_full(e));
	}
	return 0;
}

int IRefiner::
get_local_face_mark(Volume* vol, Face* f) const
{
	int sideInd = GetFaceIndex(vol, f);
	UG_COND_THROW(sideInd == -1, "Given face is not a side of the given volume.");

	if(marked_local(vol)){
		const int volLocalMark = get_local_mark(vol);
		const size_t numFaceVrts = f->num_vertices();
		Face::ConstVertexArray vrts = f->vertices();

		int vinds[MAX_FACE_VERTICES];
		for(size_t i = 0; i < numFaceVrts; ++i){
			vinds[i] = GetVertexIndex(vol, vrts[i]);
		}

		int sideMark = 0;

		for(size_t i = 0; i < numFaceVrts; ++i){
			const int edgeInd =
					vol->get_edge_index_from_vertices(
							vinds[i], vinds[(i+1)%numFaceVrts]);
			
			sideMark |= ((volLocalMark >> edgeInd) & 1) << i;
		}

		return sideMark;
	}
	else if(marked_full(vol)){
		int sideMark = 0;
		const size_t numFaceVrts = f->num_vertices();
		for(size_t i = 0; i < numFaceVrts; ++i)
			sideMark |= 1 << i;
		return sideMark;
	}
	else if(marked_closure(vol)){
		int sideMark = 0;
		const size_t numFaceEdges = f->num_edges();
		for(size_t iedge = 0; iedge < numFaceEdges; ++iedge){
			if(marked_full(const_cast<IRefiner*>(this)->grid()->get_edge(f, iedge)))
				sideMark |= 1 << iedge;
		}

		return sideMark;
	}
	return 0;
}


size_t IRefiner::num_marked_edges(std::vector<int>& numMarkedEdgesOut)
{
	#ifdef UG_PARALLEL
		std::vector<int> numLocal;
		num_marked_edges_local(numLocal);
		pcl::ProcessCommunicator com;
		com.allreduce(numLocal, numMarkedEdgesOut, PCL_RO_SUM);
	#else
		num_marked_edges_local(numMarkedEdgesOut);
	#endif

	size_t total = 0;
	for(size_t i = 0; i < numMarkedEdgesOut.size(); ++i){
		total += numMarkedEdgesOut[i];
	}
	return total;
}

size_t IRefiner::num_marked_faces(std::vector<int>& numMarkedFacesOut)
{
	#ifdef UG_PARALLEL
		std::vector<int> numLocal;
		num_marked_faces_local(numLocal);
		pcl::ProcessCommunicator com;
		com.allreduce(numLocal, numMarkedFacesOut, PCL_RO_SUM);
	#else
		num_marked_faces_local(numMarkedFacesOut);
	#endif

	size_t total = 0;
	for(size_t i = 0; i < numMarkedFacesOut.size(); ++i){
		total += numMarkedFacesOut[i];
	}
	return total;
}

size_t IRefiner::num_marked_volumes(std::vector<int>& numMarkedVolsOut)
{
	#ifdef UG_PARALLEL
		std::vector<int> numLocal;
		num_marked_volumes_local(numLocal);
		pcl::ProcessCommunicator com;
		com.allreduce(numLocal, numMarkedVolsOut, PCL_RO_SUM);
	#else
		num_marked_volumes_local(numMarkedVolsOut);
	#endif

	size_t total = 0;
	for(size_t i = 0; i < numMarkedVolsOut.size(); ++i){
		total += numMarkedVolsOut[i];
	}
	return total;
}

size_t IRefiner::num_marked_elements(std::vector<int>& numMarkedElemsOut)
{
	if(grid()){
		Grid& g = *grid();

		size_t numVolumes = g.num<Volume>();
		size_t numFaces = g.num<Face>();
		size_t numEdges = g.num<Edge>();

		#ifdef UG_PARALLEL
			pcl::ProcessCommunicator com;
			numVolumes = com.allreduce(numVolumes, PCL_RO_SUM);
			numFaces = com.allreduce(numFaces, PCL_RO_SUM);
			numEdges = com.allreduce(numEdges, PCL_RO_SUM);
		#endif

		if(numVolumes > 0)
			return num_marked_volumes(numMarkedElemsOut);
		else if(numFaces > 0)
			return num_marked_faces(numMarkedElemsOut);
		else if(numEdges > 0)
			return num_marked_edges(numMarkedElemsOut);
		else
			return 0;
	}
	return 0;
}

}// end of namespace
