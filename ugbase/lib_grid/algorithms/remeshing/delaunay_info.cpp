// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "delaunay_info.h"

namespace ug{

template <class TAAPos>
DelaunayInfo<TAAPos>::
DelaunayInfo(Grid& g, TAAPos& aaPos,
			 Grid::edge_traits::callback cbConstrainedEdge)
	: m_grid(g), m_aaPos(aaPos), m_numMarkedFaces(0), m_minAngle(0),
	  m_maxDot(1), m_cbConstrainedEdge(cbConstrainedEdge),
	  m_faceClassificationEnabled(false),
	  m_candidateRecordingEnabled(false)
{
	g.attach_to_vertices_dv(m_aCandidateMark, 0);
	m_aaMarkedVRT.access(g, m_aCandidateMark);

	g.attach_to_edges_dv(m_aCandidateMark, 0);
	m_aaMarkedEDGE.access(g, m_aCandidateMark);

	g.attach_to_faces_dv(m_aFaceInfo, NULL);
	m_aaFaceInfo.access(g, m_aFaceInfo);

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
	enable_face_classification(0);
	m_grid.unregister_observer(this);
	m_grid.detach_from_vertices(m_aCandidateMark);
	m_grid.detach_from_edges(m_aCandidateMark);
	m_grid.detach_from_faces(m_aCandidateMark);
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
mark(Vertex* vrt, bool mark)
{
	if(mark)	m_aaMarkedVRT[vrt] = 1;
	else		m_aaMarkedVRT[vrt] = 0;
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
push_candidate(Edge* e)
{
	if(!is_candidate(e)){
		m_aaMarkedEDGE[e] = 1;
		m_qEdgeCandidates.push(e);
	}
}


template <class TAAPos>
Edge* DelaunayInfo<TAAPos>::
pop_candidate()
{
	Edge* e = m_qEdgeCandidates.front();
	m_qEdgeCandidates.pop();
	m_aaMarkedEDGE[e] = 0;
	return e;
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
mark(Face* f, bool mark)
{
	FaceInfo* fi = m_aaFaceInfo[f];
	if(mark && (fi == NULL)){
		++m_numMarkedFaces;
		if(face_classification_enabled()){
			fi = m_aaFaceInfo[f] = new FaceInfo;
			fi->f = f;
			classify_face(f);
		}
		else{
			m_aaFaceInfo[f] = &m_faceMark;
		}
	}
	else if((!mark) && (fi != NULL)){
		--m_numMarkedFaces;
		if(face_classification_enabled()){
			if(is_classified(f)){
			//	we can't delete the face-info now. instead we'll only set
			//	info->f to NULL
				fi->f = NULL;
			}
			else{
			//	delete the associated face-info
				delete fi;
			}
		}
		m_aaFaceInfo[f] = NULL;
	}
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
classified_faces_left()
{
//	pop illegal entries
	while(!m_faceQueue.empty()){
		FaceInfo* fi = m_faceQueue.top();
		if(fi->f == NULL){
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
		if(fi->f == NULL){
			m_faceQueue.pop();
			delete fi;
		}
		else
			break;
	}

//	one should always call classified_faces_left before calling this method!
	assert(!m_faceQueue.empty());
	FaceInfo* fi = m_faceQueue.top();
	m_faceQueue.pop();
	fi->classified = false;
	return fi->f;
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
enable_face_classification(number minAngle)
{
	bool faceClassificationWasEnabled = face_classification_enabled();
	m_faceClassificationEnabled = (minAngle > 0);

	m_maxDot = fabs(cos(deg_to_rad(minAngle)));
	if(minAngle > 0){
		if(faceClassificationWasEnabled){
		//	careful - faceInfo->classified has to be set to false!
			m_faceQueue = FacePriorityQueue();

		//	classification was already enabled.
			for(FaceIterator iter = m_grid.begin<Face>();
				iter != m_grid.end<Face>(); ++iter)
			{
				if(is_marked(*iter)){
					m_aaFaceInfo[*iter]->classified = false;
					classify_face(*iter);
				}
			}
		}
		else{
		//	we have to create face-infos before classification
			for(FaceIterator iter = m_grid.begin<Face>();
				iter != m_grid.end<Face>(); ++iter)
			{
				if(is_marked(*iter)){
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
			if(is_marked(*iter)){
				delete m_aaFaceInfo[*iter];
				m_aaFaceInfo[*iter] = &m_faceMark;
			}
		}
		m_faceQueue = FacePriorityQueue();
	}
	m_minAngle = minAngle;
	m_faceClassificationEnabled = (minAngle > 0);
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
vertex_created(Grid* grid, Vertex* vrt, GridObject* pParent, bool replacesParent)
{
	m_aaMarkedVRT[vrt] = 1;

//	new vertices created on constrained edges shall not be marked
	if(pParent){
		if(pParent->base_object_id() == EDGE){
			if(is_constrained(static_cast<Edge*>(pParent)))
				m_aaMarkedVRT[vrt] = 0;
		}
	}
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
edge_created(Grid* grid, Edge* e, GridObject* pParent, bool replacesParent)
{
	m_aaMarkedEDGE[e] = 0;

	if(pParent){
		if(pParent->base_object_id() == EDGE){
			if(is_constrained(static_cast<Edge*>(pParent)))
				mark_as_constrained(e);
		}
	}

	if(m_candidateRecordingEnabled){
		m_recordedCandidates.push_back(e);
	}
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
face_created(Grid* grid, Face* f, GridObject* pParent, bool replacesParent)
{
//	if the new face has a parent (it should always have one if a split
//	was performed), and if that parent is marked, then we'll mark the new
//	face.
	if(pParent){
		if(pParent->base_object_id() == FACE){
			if(is_marked(static_cast<Face*>(pParent)))
				mark(f);
		}
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
	mark(f, false);
}


template <class TAAPos>
bool DelaunayInfo<TAAPos>::
is_classified(Face* f)
{
	if(face_classification_enabled()){
		if(is_marked(f)){
			return m_aaFaceInfo[f]->classified;
		}
	}
	return false;
}


template <class TAAPos>
void DelaunayInfo<TAAPos>::
classify_face(Face* f)
{
	FaceInfo* fi = m_aaFaceInfo[f];
	assert(is_marked(f));
	assert(fi != NULL);
	assert(fi != &m_faceMark);
	assert(face_classification_enabled());
	assert(fi->classified == false);

//	only triangles can be classified
	if(f->num_vertices() != 3){
		fi->classified = false;
		return;
	}

//	calculate min angle
	vector_t& v1 = m_aaPos[f->vertex(0)];
	vector_t& v2 = m_aaPos[f->vertex(1)];
	vector_t& v3 = m_aaPos[f->vertex(2)];

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
		if(!(is_constrained(e1) && is_constrained(e2)))
			highestDot = d1;
	}

	if((highestDot < m_maxDot) && (d2 > m_maxDot)){
	//	check edges
		Edge* e1 = m_grid.get_edge(f, 0);
		Edge* e2 = m_grid.get_edge(f, 1);
		if(!(is_constrained(e1) && is_constrained(e2)))
			highestDot = d2;
	}

	if((highestDot < m_maxDot) && (d3 > m_maxDot)){
	//	check edges
		Edge* e1 = m_grid.get_edge(f, 1);
		Edge* e2 = m_grid.get_edge(f, 2);
		if(!(is_constrained(e1) && is_constrained(e2)))
			highestDot = d3;
	}

	if(highestDot > m_maxDot){
		fi->classified = true;
	//	calculate the radius of the circumcenter for the priority
		vector_t cc;
		if(!TriangleCircumcenter(cc, v1, v2, v3)){
			fi->priority = 1.e12;//	just a number to increase priority
		}
		else
			fi->priority = VecDistanceSq(cc, v1);
		m_faceQueue.push(fi);
	}

//todo	if the face-queue gets too large, we have to clean it up
//		compare m_numMarkedFaces and queue->size...

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

template class DelaunayInfo<Grid::VertexAttachmentAccessor<AVector3> >;

}//	end of namespace
