// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_delaunay_info
#define __H__UG_delaunay_info

#include <queue>
#include <vector>
#include <sstream>
#include "lib_grid/lg_base.h"
#include "common/ug_config.h"

namespace ug{

/** This class intended for internal use in delaunay related algorithms.
 *
 * This helper class describes the set of triangles and edges, on which
 * delaunay triangulation shall be performed. All triangles for which the
 * delaunay triangulation shall hold, have to be marked through the mark
 * method.
 * All edges, which are candidates for a swap, should be declared as such
 * through the push_candidate method.
 *
 * Vertices created while the DelaunayInfo exists, are automatically marked, if
 * they are not created on a constrained edge.
 */
template <class TAAPos>
class UG_API DelaunayInfo : public GridObserver
{
	typedef TAAPos	AAPos;
	typedef typename TAAPos::ValueType	vector_t;

	public:
		DelaunayInfo(Grid& g, TAAPos& aaPos,
					 Grid::edge_traits::callback cbConstrainedEdge);

		~DelaunayInfo();

		Grid& grid()				{return m_grid;}
		AAPos& position_accessor()	{return m_aaPos;}

		void mark(Vertex* vrt, bool mark = true);

		bool is_marked(Vertex* vrt)		{return m_aaMarkedVRT[vrt] != 0;}

		void mark_as_constrained(Edge* e)	{m_aaMarkedEDGE[e] = 2;}

		bool is_constrained(Edge* e)	{return (m_aaMarkedEDGE[e] == 2) || m_cbConstrainedEdge(e);}

		bool is_candidate(Edge* e)		{return m_aaMarkedEDGE[e] == 1;}

		void push_candidate(Edge* e);
		Edge* pop_candidate();
		bool candidates_left()				{return !m_qEdgeCandidates.empty();}

		void mark(Face* f, bool mark = true);
		bool is_marked(Face* f)					{return m_aaFaceInfo[f] != NULL;}

		template <class TIter>
		void mark(TIter begin, TIter end)
		{
			while(begin != end) {mark(*begin); ++begin;}
		}

		bool classified_faces_left();
		Face* pop_classified_face();
		bool face_classification_enabled()		{return m_faceClassificationEnabled;}
		void enable_face_classification(number minAngle);

		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent,
									bool replacesParent);

		virtual void edge_created(Grid* grid, Edge* e,
									GridObject* pParent,
									bool replacesParent);

		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent,
									bool replacesParent);

		
		virtual void edge_to_be_erased(Grid* grid, Edge* e, Edge* replacedBy);

		virtual void face_to_be_erased(Grid* grid, Face* f, Face* replacedBy);

	///	newly created edges will be recorded as possible new candidates
	/**	All newly created edges will be added to a list of possible candidates,
	 * however, they are not added to the list of candidates until stop_candidate_recording()
	 * has been called. This is important since recorded possible candidates may be erased
	 * again from the grid (opposed to real candidates, which may not be erased).*/
		void start_candidate_recording();

	///	stops candidate recording and pushes all valid recorded edges to the list of real candidates
		void stop_candidate_recording();

	private:
		struct FaceInfo{
			FaceInfo() : f(NULL), priority(0), classified(false)	{}
			Face* f;
			number priority;
			bool classified;
		};

		struct CompareFaceInfo{
			bool operator()(const FaceInfo* fi1, const FaceInfo* fi2) const
			{
				return fi1->priority < fi2->priority;
			}
		};

		typedef Attachment<FaceInfo*> AFaceInfo;

		typedef std::priority_queue<FaceInfo*, std::vector<FaceInfo*>,
									CompareFaceInfo>
			FacePriorityQueue;


		bool is_classified(Face* f);
		void classify_face(Face* f);

	private:
		Grid& 	m_grid;
		AAPos	m_aaPos;

		AByte 	m_aCandidateMark;
		Grid::AttachmentAccessor<Vertex, AByte>	m_aaMarkedVRT;
		Grid::AttachmentAccessor<Edge, AByte>	m_aaMarkedEDGE;
		std::queue<Edge*>	m_qEdgeCandidates;

	//	pointer to m_faceMark is used to mark faces, if face-classification is disabled.
		int			m_numMarkedFaces;
		FaceInfo	m_faceMark;
		AFaceInfo	m_aFaceInfo;
		Grid::AttachmentAccessor<Face, AFaceInfo>	m_aaFaceInfo;
		FacePriorityQueue	m_faceQueue;

		number							m_minAngle;
		number							m_maxDot;
		Grid::edge_traits::callback 	m_cbConstrainedEdge;

		bool	m_faceClassificationEnabled;

	//	helper vector used during start/stop_candidate_recording
		bool	m_candidateRecordingEnabled;
		std::vector<Edge*>	m_recordedCandidates;
};

}//	end of namespace

#endif	//__H__UG_delaunay_info
