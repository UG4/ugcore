// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.09.2011 (m,d,y)

#ifndef __H__UG__delaunay_triangulation__
#define __H__UG__delaunay_triangulation__

#include <queue>
#include <vector>
#include <sstream>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"
#include "lib_grid/algorithms/polychain_util.h"
#include "common/ug_config.h"

//temporary
#include "lib_grid/file_io/file_io.h"

namespace ug
{

class UG_API DelaunayDebugSaver
{
	public:
		static DelaunayDebugSaver& inst()
		{
			static DelaunayDebugSaver dds;
			return dds;
		}

		void save(Grid& g, const char* msg)
		{
			if(m_saveEnabled){
				std::stringstream ss;
				ss << "dbg" << m_counter++ << ".ugx";

				UG_LOG(msg << ": debug-save to " << ss.str() << std::endl);
				SubsetHandler sh(g);
				AssignGridToSubset(g, sh, 0);
				SaveGridToFile(g, sh, ss.str().c_str());
			}
		}

		void enable_save(bool enable = true)	{m_saveEnabled = enable;}

	private:
		DelaunayDebugSaver() : m_saveEnabled(false), m_counter(0) {}

		bool	m_saveEnabled;
		int		m_counter;
};

inline void EnableDelaunayDebugSave(bool enable = true)
{
	DelaunayDebugSaver::inst().enable_save(enable);
}

inline void DelaunayDebugSave(Grid& g, const char* msg)
{
	DelaunayDebugSaver::inst().save(g, msg);
}


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
class DelaunayInfo : public GridObserver
{
	typedef TAAPos	AAPos;
	typedef typename TAAPos::ValueType	vector_t;

	public:
		DelaunayInfo(Grid& g, TAAPos& aaPos,
					 Grid::edge_traits::callback cbConstrainedEdge)
			: m_grid(g), m_aaPos(aaPos), m_numMarkedFaces(0), m_minAngle(0),
			  m_maxDot(1), m_cbConstrainedEdge(cbConstrainedEdge),
			  m_faceClassificationEnabled(false)
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

		~DelaunayInfo()
		{
			enable_face_classification(0);
			m_grid.unregister_observer(this);
			m_grid.detach_from_vertices(m_aCandidateMark);
			m_grid.detach_from_edges(m_aCandidateMark);
			m_grid.detach_from_faces(m_aCandidateMark);
		}

		Grid& grid()				{return m_grid;}
		AAPos& position_accessor()	{return m_aaPos;}

		void mark(VertexBase* vrt, bool mark = true)
		{
			if(mark)	m_aaMarkedVRT[vrt] = 1;
			else		m_aaMarkedVRT[vrt] = 0;
		}

		bool is_marked(VertexBase* vrt)		{return m_aaMarkedVRT[vrt] != 0;}

		void mark_as_constrained(EdgeBase* e)	{m_aaMarkedEDGE[e] = 2;}

		bool is_constrained(EdgeBase* e)	{return (m_aaMarkedEDGE[e] == 2) || m_cbConstrainedEdge(e);}

		bool is_candidate(EdgeBase* e)		{return m_aaMarkedEDGE[e] == 1;}

		void push_candidate(EdgeBase* e)
		{
			if(!is_candidate(e)){
				m_aaMarkedEDGE[e] = 1;
				m_qEdgeCandidates.push(e);
			}
		}

		EdgeBase* pop_candidate()
		{
			EdgeBase* e = m_qEdgeCandidates.front();
			m_qEdgeCandidates.pop();
			m_aaMarkedEDGE[e] = 0;
			return e;
		}

		bool candidates_left()				{return !m_qEdgeCandidates.empty();}

		bool is_marked(Face* f)					{return m_aaFaceInfo[f] != NULL;}

		void mark(Face* f, bool mark = true)
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

		template <class TIter>
		void mark(TIter begin, TIter end)
		{
			while(begin != end) {mark(*begin); ++begin;}
		}

		bool classified_faces_left()
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

		Face* pop_classified_face()
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

		bool face_classification_enabled()		{return m_faceClassificationEnabled;}

		void enable_face_classification(number minAngle)
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

		virtual void vertex_created(Grid* grid, VertexBase* vrt,
									GeometricObject* pParent,
									bool replacesParent)
		{
			m_aaMarkedVRT[vrt] = 1;

		//	new vertices created on constrained edges shall not be marked
			if(pParent){
				if(pParent->base_object_id() == EDGE){
					if(is_constrained(static_cast<EdgeBase*>(pParent)))
						m_aaMarkedVRT[vrt] = 0;
				}
			}
		}

		virtual void edge_created(Grid* grid, EdgeBase* e,
									GeometricObject* pParent,
									bool replacesParent)
		{
			m_aaMarkedEDGE[e] = 0;

			if(pParent){
				if(pParent->base_object_id() == EDGE){
					if(is_constrained(static_cast<EdgeBase*>(pParent)))
						mark_as_constrained(e);
				}
			}
		}

		virtual void face_created(Grid* grid, Face* f,
									GeometricObject* pParent,
									bool replacesParent)
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

		virtual void face_to_be_erased(Grid* grid, Face* f, Face* replacedBy)
		{
		//	unmark the face.
			mark(f, false);
		}

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


		bool is_classified(Face* f)
		{
			if(face_classification_enabled()){
				if(is_marked(f)){
					return m_aaFaceInfo[f]->classified;
				}
			}
			return false;
		}

		void classify_face(Face* f)
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
				EdgeBase* e1 = m_grid.get_edge(f, 0);
				EdgeBase* e2 = m_grid.get_edge(f, 2);
				if(!(is_constrained(e1) && is_constrained(e2)))
					highestDot = d1;
			}

			if((highestDot < m_maxDot) && (d2 > m_maxDot)){
			//	check edges
				EdgeBase* e1 = m_grid.get_edge(f, 0);
				EdgeBase* e2 = m_grid.get_edge(f, 1);
				if(!(is_constrained(e1) && is_constrained(e2)))
					highestDot = d2;
			}

			if((highestDot < m_maxDot) && (d3 > m_maxDot)){
			//	check edges
				EdgeBase* e1 = m_grid.get_edge(f, 1);
				EdgeBase* e2 = m_grid.get_edge(f, 2);
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

	private:
		Grid& 	m_grid;
		AAPos	m_aaPos;

		AByte 	m_aCandidateMark;
		Grid::AttachmentAccessor<VertexBase, AByte>	m_aaMarkedVRT;
		Grid::AttachmentAccessor<EdgeBase, AByte>	m_aaMarkedEDGE;
		std::queue<EdgeBase*>	m_qEdgeCandidates;

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
};


////////////////////////////////////////////////////////////////////////////////
template <class TAAPos>
bool MakeDelaunay(DelaunayInfo<TAAPos>& info)
{
	using namespace std;
	typedef typename TAAPos::ValueType vector_t;

	Grid& grid = info.grid();
	Face* nbrFaces[2];
	vector<EdgeBase*> edges;

//UG_LOG("MakeDelaunay...\n");
	TAAPos& aaPos = info.position_accessor();

	while(info.candidates_left()){
		EdgeBase* e = info.pop_candidate();

	//	we only perform swaps on regular manifolds.
		if(GetAssociatedFaces(nbrFaces, grid, e, 2) == 2){
		//	make sure that both neighbors are triangles
			if(nbrFaces[0]->num_vertices() != 3 || nbrFaces[1]->num_vertices() != 3)
				continue;

		//	This section is just temporary...
			VertexBase* conVrt0 = GetConnectedVertex(e, nbrFaces[0]);
			VertexBase* conVrt1 = GetConnectedVertex(e, nbrFaces[1]);

			vector_t& v0 = aaPos[e->vertex(0)];
			vector_t& v1 = aaPos[e->vertex(1)];
			vector_t& v2 = aaPos[conVrt0];
			vector_t& v3 = aaPos[conVrt1];

			vector_t cc1, cc2;

			bool cc1_ok = TriangleCircumcenter(cc1, v0, v1, v2);
			bool cc2_ok = TriangleCircumcenter(cc2, v1, v0, v3);

			number r1Sq = VecDistanceSq(cc1, v0);
			number r2Sq = VecDistanceSq(cc2, v0);

//			UG_LOG("  checking edge: " << CalculateCenter(e, aaPos) << "\n");
/*					<< "cc1_ok: " << cc1_ok << ", cc2_ok: " << cc2_ok
					<< ", cc1: " << cc1 << ", cc2: " << cc2
					<< ", r1: " << sqrt(r1Sq) << ", r2: " << sqrt(r2Sq) << endl);
*/

		//	for stability reasons, we're checking against the smaller circle
			if(cc1_ok){
				if(cc2_ok){
					if(r1Sq <= r2Sq){
						if(r1Sq <= VecDistanceSq(cc1, v3))
							continue; //the edge is fine
					}
					else{
						if(r2Sq <= VecDistanceSq(cc2, v2))
							continue; //the edge is fine
					}
				}
				else{
					if(r1Sq <= VecDistanceSq(cc1, v3))
						continue; //the edge is fine
				}
			}
			else if(cc2_ok){
				if(r2Sq <= VecDistanceSq(cc2, v2))
					continue; //the edge is fine
			}
			else{
				UG_LOG("TriangleCircumcenter failed! Excpect non-delaunay output!\n");
				UG_LOG("  This is most likely caused by two degenerated triangles which "
						"share an edge.\n");
				UG_LOG("edge center: " << CalculateCenter(e, aaPos) << endl);
				return false;
				//continue;
			}

//			UG_LOG("  swap-it!\n");

		//	before swapping, we have to make sure, that the generated triangle
		//	won't be degenerated.
		//	this is a costly operation... we check whether both circumcenters
		//	of the new triangles can be calculated...
/*
			if(!(TriangleCircumcenter(cc1, v0, v2, v3)
				 && TriangleCircumcenter(cc2, v2, v1, v3)))
			{
			//	we have to abort the swap
				continue;
			}
*/
		//	ok - everything is fine. Now swap the edge
			EdgeBase* eNew = SwapEdge(grid,  e);

			if(!eNew){
				UG_LOG("An edge-swap failed. Expect degenerated or flipped triangles "
						"and a non-delaunay output!\n");
				UG_LOG("edge center: " << CalculateCenter(e, aaPos) << endl);
				return false;
				//continue;
			}

			e = eNew;

			DelaunayDebugSave(grid, "Edge Swapped");

		//	all edges of associated triangles are candidates again (except e)
			GetAssociatedFaces(nbrFaces, grid, e, 2);
			for(size_t i = 0; i < 2; ++i){
				CollectAssociated(edges, grid, nbrFaces[i]);
				for(size_t j = 0; j < edges.size(); ++j){
					if(edges[j] != e && !(info.is_constrained(edges[j])
										  || info.is_candidate(edges[j])))
					{
						info.push_candidate(edges[j]);
					}
				}
			}
		}
	}
	return true;
}

inline
bool DelaunayLineLineIntersection(vector3& vOut,
								const vector3& lineFrom, const vector3& lineTo,
								const vector3& edgeVrt1, const vector3& edgeVrt2,
								vector3 areaNormal)
{
	number t;

	vector3 lineDir, edgeDir, planeNormal;
	VecSubtract(lineDir, lineTo, lineFrom);
	VecNormalize(lineDir, lineDir);
	VecSubtract(edgeDir, edgeVrt2, edgeVrt1);
	VecNormalize(edgeDir, edgeDir);
	VecCross(planeNormal, edgeDir, areaNormal);
	if(RayPlaneIntersection(vOut, t, lineFrom, lineDir, edgeVrt1, planeNormal)){
		if(t > 0 && VecDistanceSq(lineFrom, vOut) < VecDistanceSq(lineFrom, lineTo) + SMALL_SQ){
			return true;
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////
///	Transforms the given triangle-set into a delaunay set
/** Creates a delaunay triangulation. If a minAngle greater than 0 is specified,
 * then additional vertices are introduced, if required to generate fulfill
 * the min-angle-condition.
 */
template <class TriIter, class TAAPos>
bool QualityGridGeneration(Grid& grid, TriIter trisBegin, TriIter trisEnd,
						   TAAPos& aaPos, number minAngle = 0,
				  	  	   Grid::edge_traits::callback cbConstrainedEdge = Grid::edge_traits::cb_consider_none,
				  	  	   int maxSteps = -1/*remove this*/)
{
	using namespace std;
	typedef typename TAAPos::ValueType vector_t;

	int stepCount = 0;

//	helper to collect neighbors
	Face* nbrFaces[2];
	queue<VertexBase*> qvrts; // used during splits
	vector<VertexBase*> vrts;
	vector<EdgeBase*> edges;
	vector<EdgeBase*> closeEdges; // used during splits
	vector<Face*> faces;

//	set up a delaunay-info structure
	DelaunayInfo<TAAPos> info(grid, aaPos, cbConstrainedEdge);

//	first mark all triangles
//	mark all vertices of those tris
	for(TriIter iter = trisBegin; iter != trisEnd; ++iter){
		Face* f = *iter;
		if(f->num_vertices() != 3)
			continue;

		info.mark(f);
		for(size_t i = 0; i < 3; ++i){
			info.mark(f->vertex(i));
		}
	}

//	Collect all candidates for flips (only edges with two neighbors, both marked).
	for(TriIter triIter = trisBegin; triIter != trisEnd; ++triIter){
		Face* t = *triIter;
		CollectEdges(edges, grid, t);
		for(size_t i = 0; i < edges.size(); ++i){
			EdgeBase* e = edges[i];
		//	unmark associated vertices of constrained edges
		//...
			if(info.is_constrained(e)){
			//	unmark associated vertices
				info.mark(e->vertex(0), false);
				info.mark(e->vertex(1), false);
			}
			else{
				if(!info.is_candidate(e)){
					int numNbrs = GetAssociatedFaces(nbrFaces, grid, e, 2);
				//	two neighbors, both marked
					if(numNbrs == 2 && info.is_marked(nbrFaces[0])
						&& info.is_marked(nbrFaces[1]))
					{
					//	the edge is a flip candidate
						info.push_candidate(e);
					}
					else{
					//	unmark associated vertices
						info.mark(e->vertex(0), false);
						info.mark(e->vertex(1), false);
						info.mark_as_constrained(e);
					}
				}
			}
		}
	}

	MakeDelaunay(info);

	if(minAngle > 0 && maxSteps != 0){
		info.enable_face_classification(minAngle);

	//	while there are faces left which have to be improved
		while(info.classified_faces_left()){
			++stepCount;
			/*if(stepCount == 488){
				EnableDelaunayDebugSave();
			}*/

			Face* f = info.pop_classified_face();
			if(f->num_vertices() != 3)
				continue;

		//todo: the normal is only required for 3d-types. indeed this will crash
		//		for 2d. This should be moved to a Ray-Line_Intersection3d test.
			vector3 faceNormal;
			CalculateNormal(faceNormal, f, aaPos);
		//	we can't operate on degenerated faces. Let's hope, that this face
		//	will improve during improvement of some non-degenerated face.
			if(VecLengthSq(faceNormal) < SMALL)
				return false;

			vector_t faceCenter = CalculateCenter(f, aaPos);

		//	if two or more edges of this triangle are constrained, we'll ignore
		//	it to avoid infinite recursion.
			CollectAssociated(edges, grid, f);
			int numConstrained = 0;
			for(size_t i = 0; i < edges.size(); ++i){
				if(info.is_constrained(edges[i]))
					++numConstrained;
			}
			if(numConstrained > 1)
				continue;

		//	calculate triangle-circumcenter
			vector_t& v0 = aaPos[f->vertex(0)];
			vector_t& v1 = aaPos[f->vertex(1)];
			vector_t& v2 = aaPos[f->vertex(2)];
			vector_t cc;

			if(!TriangleCircumcenter(cc, v0, v1, v2)){
				UG_LOG("Couldn't calculate triangle-circumcenter. Expect unexpected results!\n");
				UG_LOG("triangle: " << faceCenter << "\n");
				//SaveGridToFile(grid, "delaunay_debug.ugx");
				return false;
				//continue;
			}

		//	locate the triangle which contains cc. Do this by traversing edges
		//	as required. Note that since the delaunay property holds, we're only
		//	traversing edges in a circle, which does not contain any vertices.
			EdgeBase* lastTraversedEdge = NULL;
			Face* curFace = f;
			//vector_t startPos = CalculateCenter(f, aaPos);
			vector_t rayDir;
			VecSubtract(rayDir, cc, faceCenter);
			VertexBase* pointInserted = NULL;

			while(pointInserted == NULL){
				//UG_LOG("curTri: " << CalculateCenter(curFace, aaPos) << "\n");
			//	to make things as robust as possible, we'll always intersect the
			//	same line with each edge
				EdgeBase* nextEdge = NULL;
				bool split = false;

				CollectAssociated(edges, grid, curFace);
				for(size_t i = 0; i < edges.size(); ++i){
					EdgeBase* e = edges[i];
					if(e == lastTraversedEdge)
						continue;

					vector_t a;
					if(DelaunayLineLineIntersection(a, faceCenter, cc,
							aaPos[e->vertex(0)], aaPos[e->vertex(1)], faceNormal))
					{
					//	this edge will be traversed next.
						nextEdge = e;
					//	check whether e has to be split, to avoid bad aspect ratios
						split = (VecDistanceSq(faceCenter, a) > VecDistanceSq(faceCenter, cc) - SMALL_SQ);
						break;
					}/*
					// This else-condition is not necessary - if they are parallel, that's not a problem...
					else{
						UG_LOG("Ray-Plane intersection failed in step " << stepCount << "... aborting.\n");
						UG_LOG("  curTri: " << CalculateCenter(curFace, aaPos)
								<< ", curEdge: " << CalculateCenter(e, aaPos) << endl);
						UG_LOG("  face-normal: " << faceNormal << ", plane-normal: " << planeNormal
								<< "  ray-dir: " << rayDir << endl);
						return false;
					}*/
				}

				if(nextEdge){
				//	check whether the edge has to be splitted
					bool isConstrained = info.is_constrained(nextEdge);
					split |= isConstrained;
					if(!split){
						int numNbrs = GetAssociatedFaces(nbrFaces, grid, nextEdge, 2);
						if(numNbrs != 2)
							split = true;
						else{
						//	get the next face
							if(nbrFaces[0] == curFace)
								curFace = nbrFaces[1];
							else
								curFace = nbrFaces[0];
						//	and set the last traversed edge
							lastTraversedEdge = nextEdge;
						}
					}

					if(split){
						//UG_LOG("SPLIT\n");
						vector_t center = CalculateCenter(nextEdge, aaPos);

						VertexBase* vrt0 = nextEdge->vertex(0);
						VertexBase* vrt1 = nextEdge->vertex(1);
						number radiusSq = VecDistanceSq(center, aaPos[vrt0]);
						number radius = sqrt(radiusSq);
						pointInserted = SplitEdge<Vertex>(grid, nextEdge, false);

						if(isConstrained){
							//UG_LOG("IS CONSTRAINED\n");
							aaPos[pointInserted] = center;

						//	we have to erase all vertices, which are marked and
						//	which are closer to the inserted point than the
						//	radius calculated above.
							assert(qvrts.empty());
							qvrts.push(pointInserted);

							vrts.clear();
							closeEdges.clear();

							grid.begin_marking();
							while(!qvrts.empty()){
								VertexBase* vrt = qvrts.front();
								qvrts.pop();

							//	collect all edges connected to vrt
								CollectAssociated(edges, grid, vrt);
							//	if the edge is not yet marked, we'll add it to the
							//	closeEdges list.
								for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
									EdgeBase* e = edges[i_edge];
									if(!grid.is_marked(e)){
									//	check whether the edge intersects the critical circle
										if(DistancePointToLine(center, aaPos[e->vertex(0)], aaPos[e->vertex(1)])
											< radiusSq)
										{
										//	if the edge is a constrained edge, we'll push it to
										//	closeEdges
											if(info.is_constrained(e) && !EdgeContains(e, pointInserted))
												closeEdges.push_back(e);

											grid.mark(e);

										//	the connected vertex has to be pushed to
										//	the queue, regardless whether it lies in the circle
										//	or not, since an edge from this vertex could
										//	reenter the circle
											VertexBase* vcon = GetConnectedVertex(e, vrt);
											if(grid.is_marked(vcon))
												continue;
											grid.mark(vcon);
											qvrts.push(vcon);

										//	vrt0 and vrt1 won't be removed anyways
											if(vcon == vrt0 || vcon == vrt1)
												continue;

										//	if the new vertrex is marked (created during remeshing)
										//	and lies in the circle, then we'll push it to vrts
											//UG_LOG("is_marked: " << info.is_marked(vcon) << "rad: " << radiusSq
											//	<< ", testvrt: " << VecDistanceSq(aaPos[vcon], center) << "\n");
											if(VecDistanceSq(aaPos[vcon], center) < radiusSq){
												if(info.is_marked(vcon))
													vrts.push_back(vcon);
											}
										}
									}
								}
							}
							grid.end_marking();

							//UG_LOG(vrts.size() << " vertices have to be deleted!\n");

						//	vrts now contains all vertices which lie in a circle of radius
						//	'radiusSq' around the center. We now have to check for each,
						//	whether it is visible from the center, by checking whether
						//	the line of sight intersects a constrained edge in
						//	closeEdges.
						//	if it is visible, then we'll delete it and locally remesh
						//	the grid.
							for(size_t i_vrts = 0; i_vrts < vrts.size(); ++i_vrts){
								VertexBase* vrt = vrts[i_vrts];
								vector_t& v = aaPos[vrt];
								//UG_LOG("ERASING VRT at " << v << "...\n");
								bool intersects = false;
								//UG_LOG("  checking for intersections\n");
								for(size_t i_edge = 0; i_edge < closeEdges.size(); ++i_edge){
									vector_t& ev0 = aaPos[closeEdges[i_edge]->vertex(0)];
									vector_t& ev1 = aaPos[closeEdges[i_edge]->vertex(1)];
									vector_t a;
									if(DelaunayLineLineIntersection(a, faceCenter, v, ev0, ev1, faceNormal)){
									//if(LineLineIntersection3d(a, b, center, v, ev0, ev1)){
										//UG_LOG("  intersection!\n");
										intersects = true;
										break;
									}
								}

								if(intersects)
									continue;

								//UG_LOG("  gathering surrounding edges\n");
							//	no intersection detected. erase the vertex and perform
							//	local retriangulation. Add surrounding edges to
							//	the 'edges' array
							//	Store one face, which will be the parent for all new faces
								Face* parent = NULL;
								edges.clear();
								CollectAssociated(faces, grid, vrt);
								for(size_t i_face = 0; i_face < faces.size(); ++i_face){
									Face* f = faces[i_face];
									parent = f;
									for(size_t i = 0; i < f->num_edges(); ++i){
										EdgeBase* e = grid.get_edge(f, i);
										if(!EdgeContains(e, vrt)){
											edges.push_back(e);
											//UG_LOG("surrounding edge: " << CalculateCenter(e, aaPos) << endl);
										}
									}
								}
								assert(parent);

								//UG_LOG("  retriangulation\n");
							//	perform retriangulation
							//	get the vertices of the poly-chain
								//UG_LOG("1");
								std::vector<VertexBase*> vrtPolyChain;
								CreatePolyChain(vrtPolyChain, grid, edges.begin(), edges.end());
								//UG_LOG("2");
								std::vector<vector_t> posPolyChain(vrtPolyChain.size());
								std::vector<int> edgeChain, faceIndices;
								edgeChain.reserve(vrtPolyChain.size() * 2);
								//UG_LOG("3");
								size_t numVrts = vrtPolyChain.size();
								for(size_t i = 0; i < numVrts; ++i){
									//UG_LOG(aaPos[vrtPolyChain[i]]);
									posPolyChain[i] = aaPos[vrtPolyChain[i]];
									edgeChain.push_back(i);
									edgeChain.push_back((i+1) % numVrts);
								}
/*
								("edges.size: " << edges.size() << ", posPolyChain.size: " << posPolyChain.size());
								for(size_t i = 0; i < edgeChain.size(); i+=2){
									UG_LOG(" (" << edgeChain[i] << ", " << edgeChain[i+1] << ")");
								}
								UG_LOG(endl);
*/

								//UG_LOG("4");

								if(!TriangleFill_SweepLine(faceIndices, posPolyChain, edgeChain)){
									UG_LOG("Sweepline failed in step " << stepCount << "... aborting\n");
									UG_LOG("  While examining face " << faceCenter << endl);
									return false;
								}
								//UG_LOG("5");
							//	create the faces
								for(size_t i = 0; i < faceIndices.size(); i+=3){
									grid.create<Triangle>(
											TriangleDescriptor(vrtPolyChain[faceIndices[i]],
															   vrtPolyChain[faceIndices[i+1]],
															   vrtPolyChain[faceIndices[i+2]]),
											parent);
								}
								//UG_LOG("6\n");

								//UG_LOG("  erase\n");
							//	now erase vrt
								grid.erase(vrt);

								DelaunayDebugSave(grid, "After TriangleFill");
							}

							//UG_LOG("adding new candidates\n");
						//	now add all edges which lie inside the circle as candidates
							qvrts.push(pointInserted);

							grid.begin_marking();
							while(!qvrts.empty()){
								VertexBase* vrt = qvrts.front();
								qvrts.pop();

							//	collect all edges connected to vrt
								CollectAssociated(edges, grid, vrt);
							//	if the edge is not yet marked, we'll make it a candidate
								for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
									EdgeBase* e = edges[i_edge];
									if(!grid.is_marked(e)){
									//	check whether the edge intersects the critical circle
										if(DistancePointToLine(center, aaPos[e->vertex(0)], aaPos[e->vertex(1)])
											< radius)
										{
										//	if the edge is not a constrained edge, we'll make it a candidate
											if(!info.is_constrained(e)){
											//	we'l only mark it, if both associated faces are marked.
												Face* nbrs[2];
												if(GetAssociatedFaces(nbrs, grid, e, 2) == 2){
													if(info.is_marked(nbrs[0])
													   && info.is_marked(nbrs[1]))
													{
														info.push_candidate(e);
													}
												}
											}

											grid.mark(e);

										//	the connected vertex has to be pushed to
										//	the queue, regardless whether it lies in the circle
										//	or not, since an edge from this vertex could
										//	reenter the circle
											VertexBase* vcon = GetConnectedVertex(e, vrt);
											qvrts.push(vcon);
										}
									}
								}
							}
							grid.end_marking();
						}
						else{
							aaPos[pointInserted] = cc;
						}
					}
				}
				else{
				//	we found the triangle, which contains cc. Insert the point
				//	and perform local delaunay (not necessarily local...).
				//	todo: make sure, that cc really lies in curTri
				//...
					VertexBase* vrt0 = curFace->vertex(0);
					VertexBase* vrt1 = curFace->vertex(1);
					VertexBase* vrt2 = curFace->vertex(2);

				//UG_LOG("creating new elements...\n");
					VertexBase* vrt = *grid.create<Vertex>(curFace);
					aaPos[vrt] = cc;
					grid.create<Triangle>(TriangleDescriptor(vrt0, vrt1, vrt), curFace);
					grid.create<Triangle>(TriangleDescriptor(vrt1, vrt2, vrt), curFace);
					grid.create<Triangle>(TriangleDescriptor(vrt2, vrt0, vrt), curFace);
				//UG_LOG("erasing old face\n");
					grid.erase(curFace);

					DelaunayDebugSave(grid, "After InsertPoint");

					pointInserted = vrt;
				}

				if(pointInserted){
					//UG_LOG("temp-save to delaunay_debug.ugx\n");
					//SaveGridToFile(grid, "delaunay_debug.ugx");

					//UG_LOG("adding candidates\n");
					//	if a vertex was inserted, we'll have to perform a delaunay step.
					//	mark all edges which are connected to a triangle, which is conencted
					//	to vrt.
						CollectAssociated(faces, grid, pointInserted);
						for(size_t i_face = 0; i_face < faces.size(); ++i_face){
							CollectAssociated(edges, grid, faces[i_face]);
							for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
								if(!info.is_constrained(edges[i_edge]))
									info.push_candidate(edges[i_edge]);
							}
						}

					//UG_LOG("redelaunaylizing\n");
						if(!MakeDelaunay(info)){
							UG_LOG("Make Delaunay failed in step " << stepCount << ".\n");
							UG_LOG("  While examining face " << faceCenter << endl);
							return false;
						}
					//UG_LOG("adaption step completed!\n");
				}
			}

/*
		//	check whether an illegal triangle has been inserted
			for(TriangleIterator iter = grid.begin<Triangle>();
				iter != grid.end<Triangle>(); ++iter)
			{
				vector3 n;
				CalculateNormal(n, *iter, aaPos);
				if(n.z <= 0){
					UG_LOG("ATTENTION: Illegal triangle created in step " << stepCount << "\n");
					return false;
				}
			}
*/
			if(stepCount == maxSteps)
				break;
/*
			if((stepCount % 1000) == 0){
				UG_LOG("temp-save to delaunay_debug.ugx in step " << stepCount << endl);
				//DelaunayDebugSave(grid, "");
				SaveGridToFile(grid, "delaunay_debug.ugx");
			}
*/
		}
	}
	return true;
}

}//	end of namespace

#endif
