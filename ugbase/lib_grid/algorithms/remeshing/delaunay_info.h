// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_delaunay_info
#define __H__UG_delaunay_info

#include <queue>
#include <vector>
#include <sstream>
#include "common/ug_config.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/attachment_util.h"

namespace ug{

/** This class intended for internal use in delaunay related algorithms.
 */
template <class TAAPos>
class UG_API DelaunayInfo : public GridObserver
{
	typedef TAAPos	AAPos;
	typedef typename TAAPos::ValueType	vector_t;

	public:
		enum Mark{
			NONE,
			INNER,
			NEW_INNER,
			SEGMENT,
			NEW_SEGMENT,
			DART,
			SHELL
		};

		DelaunayInfo(Grid& g, TAAPos& aaPos,
					 Grid::edge_traits::callback cbConstrainedEdge);

		~DelaunayInfo();

		Grid& grid()				{return m_grid;}
		AAPos& position_accessor()	{return m_aaPos;}


	/**	\warning	init_marks may currently only be called once! Undefined behaviour
	 *				if called multiple times or after marks have already been set.
	 *	\todo	think about adding a cleanup at the beginning of init_marks (costly)*/
		template <class TIter>
		void init_marks(TIter trisBegin, TIter trisEnd, bool pushFlipCandidates);

	////////////////////////////////////////////////////////////////
	//	MARKS
		void set_mark(Vertex* v, Mark m)		{m_aaMark[v] = m;}
		Mark mark(Vertex* v) const				{return static_cast<Mark>(m_aaMark[v]);}

		void set_mark(Edge* e, Mark m)			{m_aaMark[e] = m;}
		Mark mark(Edge* e) const				{return static_cast<Mark>(m_aaMark[e]);}

		void set_mark(Face* f, Mark m);
		Mark mark(Face* f) const				{return static_cast<Mark>(m_aaMark[f]);}

		template <class TIter>
		void set_mark(TIter begin, TIter end, Mark m)
		{
			while(begin != end) {set_mark(*begin, m); ++begin;}
		}

	///	returns true if the specified element is either marked as INNER or as NEW_INNER
		template <class TElem>
		bool is_inner(TElem* e);

	///	returns true if the specified element is either marked as SEGMENT, NEW_SEGMENT, DART or SHELL
		template <class TElem>
		bool is_segment(TElem* e);

	///	returns true if the specified edge is a segment and is connected to a DART vertex
		bool is_dart_segment(Edge* e);

	///	returns true if the specified edge is a new segment and is connected to a DART vertex
		bool is_new_dart_segment(Edge* e);

	///	returns true if the specified edge is a segment and connects a DART and a SHELL vertex
		bool is_dart_shell_segment(Edge* e);

	///	returns true if the face can be classified
	/**	Faces which contain two shell vertices whose subtended angle is small
	 * are not classifiable.*/
		bool is_classifiable(Face* f);

	////////////////////////////////////////////////////////////////
	//	CANDIDATES
	/**	candidates are used during MakeDelaunay to define the set of edges which
	 * may have to be swapped to obtain a delaunay triangulation.
	 * \{ */
		bool is_candidate(Edge* e)		{return m_aaCandidateEDGE[e];}
		void push_candidate(Edge* e);
		Edge* pop_candidate();
		bool candidates_left()			{return !m_qEdgeCandidates.empty();}
	/** \} */

	///	newly created edges will be recorded as possible new candidates
	/**	All newly created edges will be added to a list of possible candidates,
	 * however, they are not added to the list of candidates until stop_candidate_recording()
	 * has been called. This is important since recorded possible candidates may be erased
	 * again from the grid (opposed to real candidates, which may not be erased).*/
		void start_candidate_recording();

	///	stops candidate recording and pushes all valid recorded edges to the list of real candidates
		void stop_candidate_recording();


	////////////////////////////////////////////////////////////////
	//	FACE-CLASSIFICATION
	/**	Face classification is used to define and order the set of faces which
	 * whose circumcenter may have to be inserted in order to obtain a delaunay-
	 * triangulation with prescribed minimal angle during QualityGridGeneration.
	 * Face classification is not required for MakeDelaunay.
	 * \{ */
		bool classified_faces_left();
		Face* pop_classified_face();
		bool face_classification_enabled()		{return m_aaFaceInfo.valid();}
		void enable_face_classification(number minAngle);

		number min_angle() const				{return m_minAngle;}
		number max_dot() const					{return m_maxDot;}
	/** \} */

	////////////////////////////////////////////////////////////////
	//	GRID-OBERSERVER CALLBACKS
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
		bool classify_face(Face* f);

	private:
		Grid& 	m_grid;
		AAPos	m_aaPos;

		AByte 	m_aMark;
		MultiElementAttachmentAccessor<AByte>		m_aaMark;

		AFaceInfo	m_aFaceInfo;
		Grid::AttachmentAccessor<Face, AFaceInfo>	m_aaFaceInfo;
		FacePriorityQueue	m_faceQueue;

		number							m_minAngle;
		number							m_maxDot;
		Grid::edge_traits::callback 	m_cbConstrainedEdge;

	//	candidate edges
		ABool									m_aCandidate;
		Grid::AttachmentAccessor<Edge, ABool>	m_aaCandidateEDGE;
		std::queue<Edge*>						m_qEdgeCandidates;

	//	helper vector used during start/stop_candidate_recording
		bool	m_candidateRecordingEnabled;
		std::vector<Edge*>	m_recordedCandidates;
};

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "delaunay_info_impl.h"


#endif	//__H__UG_delaunay_info
