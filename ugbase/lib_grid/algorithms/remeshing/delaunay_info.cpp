/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "delaunay_info.h"

namespace ug{

template <class TAAPos>
DelaunayInfo<TAAPos>::
DelaunayInfo(Grid& g, TAAPos& aaPos,
			 Grid::edge_traits::callback cbConstrainedEdge)
	: m_grid(g), m_aaPos(aaPos), m_minAngle(0),
	  m_maxDot(1), m_cbConstrainedEdge(cbConstrainedEdge),
	  m_candidateRecordingEnabled(false)
{
	g.attach_to_vertices_dv(m_aMark, NONE);
	g.attach_to_edges_dv(m_aMark, NONE);
	g.attach_to_faces_dv(m_aMark, NONE);
	m_aaMark.access(g, m_aMark, true, true, true, false);
	
	g.attach_to_edges_dv(m_aCandidate, false);
	m_aaCandidateEDGE.access(g, m_aCandidate);

	if(!g.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in DelaunayInfo: enabling FACEOPT_AUTOGENERATE_EDGES.\n");
		g.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	g.register_observer(this, OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER | OT_FACE_OBSERVER);
}

template <class TAAPos>
DelaunayInfo<TAAPos>::
~DelaunayInfo()
{
	m_grid.unregister_observer(this);
	enable_face_classification(0);
	m_grid.detach_from_vertices(m_aMark);
	m_grid.detach_from_edges(m_aMark);
	m_grid.detach_from_faces(m_aMark);
	m_grid.detach_from_edges(m_aCandidate);
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
set_mark(Face* f, typename DelaunayInfo<TAAPos>::Mark mark)
{
	m_aaMark[f] = mark;
	if(face_classification_enabled()){
		FaceInfo* fi = m_aaFaceInfo[f];
		if((mark == INNER) || (mark == NEW_INNER)){
			if(!fi){
				m_aaFaceInfo[f] = fi = new FaceInfo();
				fi->f = f;
			}

			if(!fi->classified)
				classify_face(f);
		}
		else if(fi){
			if(fi->classified){
			//	we can't delete the face-info now. instead we'll only set fi->f to nullptr
			//	so that it can be identified as invalid in m_faceQueue. Clean-up is performed later.
				fi->f = nullptr;
			}
			else{
			//	delete the associated face-info
				delete fi;
				m_aaFaceInfo[f] = nullptr;
			}
		}
	}
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
is_dart_segment(Edge* e)
{
	Edge::ConstVertexArray vrts = e->vertices();
	return is_segment(e) && (mark(vrts[0]) == DART || mark(vrts[1]) == DART);
}

template <class TAAPos>
bool DelaunayInfo<TAAPos>::
is_new_dart_segment(Edge* e)
{
	Edge::ConstVertexArray vrts = e->vertices();
	return mark(e) == NEW_SEGMENT && (mark(vrts[0]) == DART || mark(vrts[1]) == DART);
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
is_dart_shell_segment(Edge* e)
{
	Edge::ConstVertexArray vrts = e->vertices();
	return is_segment(e)
			&& ((mark(vrts[0]) == DART && mark(vrts[1]) == SHELL)
			    || (mark(vrts[0]) == SHELL && mark(vrts[1]) == DART));
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
push_candidate(Edge* e)
{
	if(!is_candidate(e)){
		m_aaCandidateEDGE[e] = true;
		m_qEdgeCandidates.push(e);
	}
}


template <class TAAPos>
Edge* DelaunayInfo<TAAPos>::
pop_candidate()
{
	Edge* e = m_qEdgeCandidates.front();
	m_qEdgeCandidates.pop();
	m_aaCandidateEDGE[e] = false;
	return e;
}

template <class TAAPos>
void DelaunayInfo<TAAPos>::
start_candidate_recording()
{
	UG_COND_THROW(m_candidateRecordingEnabled,
				  "Call stop_candidate_recording before recording new candidates!");
	m_recordedCandidates.clear();
	m_candidateRecordingEnabled = true;
}

template <class TAAPos>
void DelaunayInfo<TAAPos>::
stop_candidate_recording()
{
	if(m_candidateRecordingEnabled){
		for(size_t i = 0; i < m_recordedCandidates.size(); ++i){
			if(m_recordedCandidates[i]){
				push_candidate(m_recordedCandidates[i]);
			}
		}
		m_recordedCandidates.clear();
		m_candidateRecordingEnabled = false;
	}
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
enable_face_classification(number minAngle)
{
	bool faceClassificationWasEnabled = face_classification_enabled();

	m_maxDot = fabs(cos(deg_to_rad(minAngle)));
	if(minAngle > 0){
		if(faceClassificationWasEnabled){
		//	careful - faceInfo->classified has to be set to false!
			m_faceQueue = FacePriorityQueue();

		//	classification was already enabled.
			for(FaceIterator iter = m_grid.begin<Face>();
				iter != m_grid.end<Face>(); ++iter)
			{
				if(is_inner(*iter)){
					m_aaFaceInfo[*iter]->classified = false;
					classify_face(*iter);
				}
			}
		}
		else{
			m_grid.attach_to_faces_dv(m_aFaceInfo, nullptr);
			m_aaFaceInfo.access(m_grid, m_aFaceInfo);
		//	we have to create face-infos before classification
			for(FaceIterator iter = m_grid.begin<Face>();
				iter != m_grid.end<Face>(); ++iter)
			{
				if(is_inner(*iter)){
					m_aaFaceInfo[*iter] = new FaceInfo;
					m_aaFaceInfo[*iter]->f = *iter;
					classify_face(*iter);
				}
			}
		}
	}
	else if(faceClassificationWasEnabled){
	// unclassify faces.
		for(FaceIterator iter = m_grid.begin<Face>();
			iter != m_grid.end<Face>(); ++iter)
		{
			if(m_aaFaceInfo[*iter]){
				delete m_aaFaceInfo[*iter];
			}
		}
		m_faceQueue = FacePriorityQueue();
		m_grid.detach_from_faces(m_aFaceInfo);
		m_aaFaceInfo.invalidate();
	}
	m_minAngle = minAngle;
}

template <class TAAPos>
bool DelaunayInfo<TAAPos>::
classified_faces_left()
{
//	pop illegal entries
	while(!m_faceQueue.empty()){
		FaceInfo* fi = m_faceQueue.top();
		if(fi->f == nullptr){
			m_faceQueue.pop();
			delete fi;
		}
		else
			break;
	}
	return !m_faceQueue.empty();
}


template <class TAAPos>
Face* DelaunayInfo<TAAPos>::
pop_classified_face()
{
//	pop illegal entries
	while(!m_faceQueue.empty()){
		FaceInfo* fi = m_faceQueue.top();
		if(fi->f == nullptr){
			m_faceQueue.pop();
			delete fi;
		}
		else
			break;
	}

	if(m_faceQueue.empty())
		return nullptr;
	FaceInfo* fi = m_faceQueue.top();
	m_faceQueue.pop();
	fi->classified = false;
	return fi->f;
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
is_classified(Face* f)
{
	if(face_classification_enabled()){
		return m_aaFaceInfo[f]->classified;
	}
	return false;
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
is_classifiable(Face* f)
{
	int subtended = -1;
	int numShell = 0;
	for(size_t i = 0; i < 3; ++i){
		if(mark(f->vertex(i)) == SHELL)
			++numShell;
		else
			subtended = i;
	}

	vector_t v[3] = {m_aaPos[f->vertex(0)], m_aaPos[f->vertex(1)], m_aaPos[f->vertex(2)]};

	if(numShell == 3){
		number d1sq = VecDistanceSq(v[0], v[1]);
		number d2sq = VecDistanceSq(v[1], v[2]);
		number d3sq = VecDistanceSq(v[2], v[0]);
		if(d1sq < d2sq){
			if(d1sq < d3sq)
				subtended = 2;
			else
				subtended = 1;
		}
		else if(d2sq < d3sq)
			subtended = 0;
		else
			subtended = 1;
	}

	if(numShell >= 2){
		vector_t dir1, dir2;
		VecSubtract(dir1, v[(subtended+1)%3], v[subtended]);
		VecSubtract(dir2, v[(subtended+2)%3], v[subtended]);
		if(VecAngle(dir1, dir2) < PI / 3. + SMALL){
			return false;
		}
	}
	return true;
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
classify_face(Face* f)
{
	UG_COND_THROW(!face_classification_enabled(),
				  "Face classification has to be enabled before calling 'classify_face'");
	FaceInfo* fi = m_aaFaceInfo[f];
	UG_COND_THROW(!fi, "Invalid face-info attached");
	UG_COND_THROW(fi->classified, "Only unclassified faces may be classified");

//	only triangles can be classified
	if(f->num_vertices() != 3){
		fi->classified = false;
		return false;
	}

	if(!is_classifiable(f))
		return false;
	
	vector_t& v1 = m_aaPos[f->vertex(0)];
	vector_t& v2 = m_aaPos[f->vertex(1)];
	vector_t& v3 = m_aaPos[f->vertex(2)];
	
//	if at least two of the vertices are SHELL vertices, we won't classify the triangle if
//	the subtended angle to the shortest edge between shell vertices is smaller than PI/3
	// {
	// 	int subtended = -1;
	// 	int numShell = 0;
	// 	for(size_t i = 0; i < 3; ++i){
	// 		if(mark(f->vertex(i)) == SHELL)
	// 			++numShell;
	// 		else
	// 			subtended = i;
	// 	}

	// 	if(numShell == 3){
	// 		number d1sq = VecDistanceSq(v1, v2);
	// 		number d2sq = VecDistanceSq(v2, v3);
	// 		number d3sq = VecDistanceSq(v3, v1);
	// 		if(d1sq < d2sq){
	// 			if(d1sq < d3sq)
	// 				subtended = 2;
	// 			else
	// 				subtended = 1;
	// 		}
	// 		else if(d2sq < d3sq)
	// 			subtended = 0;
	// 		else
	// 			subtended = 1;
	// 	}

	// 	if(numShell >= 2){
	// 		vector_t dir1, dir2;
	// 		VecSubtract(dir1, m_aaPos[f->vertex((subtended+1)%3)], m_aaPos[f->vertex(subtended)]);
	// 		VecSubtract(dir2, m_aaPos[f->vertex((subtended+2)%3)], m_aaPos[f->vertex(subtended)]);
	// 		if(VecAngle(dir1, dir2) < PI / 3. + SMALL){
	// 			return false;
	// 		}
	// 	}
	// }

//	perform classification
//	calculate min angle
	vector_t v12, v13, v23;
	VecSubtract(v12, v2, v1);	VecNormalize(v12, v12);
	VecSubtract(v13, v3, v1);	VecNormalize(v13, v13);
	VecSubtract(v23, v3, v2);	VecNormalize(v23, v23);

	number d1 = fabs(VecDot(v12, v13));
	number d2 = fabs(VecDot(v12, v23));
	number d3 = fabs(VecDot(v13, v23));

	number highestDot = 0;
	if(d1 > m_maxDot){
	//	check edges
		Edge* e1 = m_grid.get_edge(f, 0);
		Edge* e2 = m_grid.get_edge(f, 2);
		if(!(is_segment(e1) && is_segment(e2)))
			highestDot = d1;
	}

	if((highestDot < m_maxDot) && (d2 > m_maxDot)){
	//	check edges
		Edge* e1 = m_grid.get_edge(f, 0);
		Edge* e2 = m_grid.get_edge(f, 1);
		if(!(is_segment(e1) && is_segment(e2)))
			highestDot = d2;
	}

	if((highestDot < m_maxDot) && (d3 > m_maxDot)){
	//	check edges
		Edge* e1 = m_grid.get_edge(f, 1);
		Edge* e2 = m_grid.get_edge(f, 2);
		if(!(is_segment(e1) && is_segment(e2)))
			highestDot = d3;
	}

	if(highestDot > m_maxDot){
	//	calculate the radius of the circumcenter for the priority
		vector_t cc;
		if(TriangleCircumcenter(cc, v1, v2, v3)){
			fi->priority = VecDistanceSq(cc, v1);
			fi->classified = true;
			m_faceQueue.push(fi);
		}
		else{
			// UG_LOG("Couldn't calculate triangle-circumcenter. ");
			// UG_LOG("Ignoring triangle with center at: " << CalculateCenter(f, m_aaPos) << "\n");
			return false;
		}
	}

//todo	if the face-queue gets too large, we should do some clean up
	return true;
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
vertex_created(Grid* grid, Vertex* vrt, GridObject* pParent, bool replacesParent)
{

	set_mark(vrt, NEW_INNER);
	if(pParent && (pParent->base_object_id() == EDGE)
		&& is_segment(static_cast<Edge*>(pParent)))
	{
		set_mark(vrt, NEW_SEGMENT);
	}
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
edge_created(Grid* grid, Edge* e, GridObject* pParent, bool replacesParent)
{

	set_mark(e, NEW_INNER);
	if(pParent && (pParent->base_object_id() == EDGE)
		&& is_segment(static_cast<Edge*>(pParent)))
	{
		set_mark(e, NEW_SEGMENT);
	}

	if(m_candidateRecordingEnabled)
		m_recordedCandidates.push_back(e);
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
face_created(Grid* grid, Face* f, GridObject* pParent, bool replacesParent)
{
	if(pParent && (pParent->base_object_id() == FACE)
		&& is_inner(static_cast<Face*>(pParent)))
	{
		set_mark(f, NEW_INNER);
	}
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
edge_to_be_erased(Grid* grid, Edge* e, Edge* replacedBy)
{
	if(m_candidateRecordingEnabled){
	//	if e is a recorded candidate, we'll simply overwrite it
	//	with the last valid candidate and shorten the candidates array.
		for(size_t i = 0; i < m_recordedCandidates.size(); ++i){
			if(m_recordedCandidates[i] == e){
				m_recordedCandidates[i] = m_recordedCandidates.back();
				m_recordedCandidates.pop_back();
				break;
			}
		}
	}
}

template <class TAAPos>
void DelaunayInfo<TAAPos>::
face_to_be_erased(Grid* grid, Face* f, Face* replacedBy)
{
//	unmark the face.
	set_mark(f, NONE);
}


template class DelaunayInfo<Grid::VertexAttachmentAccessor<AVector3> >;

}//	end of namespace
