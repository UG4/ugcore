//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__ATTACHMENT_UTIL__
#define __H__LIB_GRID__ATTACHMENT_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * Several methods that ease attachment-handling are grouped here.
 *
 * \defgroup lib_grid_algorithms_attachment_util attachment util
 * \ingroup lib_grid_algorithms
 * @{
 */

///	Accesses attachements in different element types at the same time.
/**	If the same attachment is attached to vertices, edges, faces and volumes of
 * a grid, then this accessor can be used to access them all at once.
 */
template <class TAttachment>
class MultiElementAttachmentAccessor
{
	public:
		typedef typename TAttachment::ValueType	ValueType;
		typedef typename attachment_value_traits<ValueType>::reference RefType;
		typedef typename attachment_value_traits<ValueType>::const_reference ConstRefType;

		MultiElementAttachmentAccessor()	{}
		MultiElementAttachmentAccessor(Grid& g, TAttachment& a, bool vrts = true,
							bool edges = true, bool faces = true, bool vols = true)
		{
			access(g, a, vrts, edges, faces, vols);
		}

		bool access(Grid& g, TAttachment& a, bool vrts = true, bool edges = true,
					bool faces = true, bool vols = true)
		{
			m_aaVrt.invalidate(); m_aaEdge.invalidate();
			m_aaFace.invalidate(); m_aaVol.invalidate();

			bool ok = true;
			if(vrts)
				ok &= m_aaVrt.access(g, a);
			if(edges)
				ok &= m_aaEdge.access(g, a);
			if(faces)
				ok &= m_aaFace.access(g, a);
			if(vols)
				ok &= m_aaVol.access(g, a);
			return ok;
		}

		bool is_valid_vertex_accessor() const		{return m_aaVrt.valid();}
		bool is_valid_edge_accessor() const			{return m_aaEdge.valid();}
		bool is_valid_face_accessor() const			{return m_aaFace.valid();}
		bool is_valid_volume_accessor() const		{return m_aaVol.valid();}

		RefType operator[](VertexBase* e)	{return m_aaVrt[e];}
		RefType operator[](EdgeBase* e)		{return m_aaEdge[e];}
		RefType operator[](Face* e)			{return m_aaFace[e];}
		RefType operator[](Volume* e)		{return m_aaVol[e];}
		RefType operator[](GeometricObject* e)
		{
			switch(e->base_object_id()){
				case VERTEX: return m_aaVrt[static_cast<VertexBase*>(e)];
				case EDGE: return m_aaEdge[static_cast<EdgeBase*>(e)];
				case FACE: return m_aaFace[static_cast<Face*>(e)];
				case VOLUME: return m_aaVol[static_cast<Volume*>(e)];
				default: UG_THROW("Unknown element type!"); break;
			}
		}

		ConstRefType operator[](VertexBase* e) const	{return m_aaVrt[e];}
		ConstRefType operator[](EdgeBase* e) const		{return m_aaEdge[e];}
		ConstRefType operator[](Face* e) const			{return m_aaFace[e];}
		ConstRefType operator[](Volume* e) const 		{return m_aaVol[e];}
		ConstRefType operator[](GeometricObject* e) const
		{
			switch(e->base_object_id()){
				case VERTEX: return m_aaVrt[static_cast<VertexBase*>(e)];
				case EDGE: return m_aaEdge[static_cast<EdgeBase*>(e)];
				case FACE: return m_aaFace[static_cast<Face*>(e)];
				case VOLUME: return m_aaVol[static_cast<Volume*>(e)];
				default: UG_THROW("Unknown element type!"); break;
			}
		}

	private:
		Grid::AttachmentAccessor<VertexBase, TAttachment>	m_aaVrt;
		Grid::AttachmentAccessor<EdgeBase, TAttachment>		m_aaEdge;
		Grid::AttachmentAccessor<Face, TAttachment>			m_aaFace;
		Grid::AttachmentAccessor<Volume, TAttachment>		m_aaVol;
};


////////////////////////////////////////////////////////////////////////
///	Instances can be used as compare operators, e.g., for std::sort.
/**	Make sure that the specified attachment is attached to the elements of
 * the given grid and that a < operator is defined for the attached values!
 */
template <class TElem, class TAttachment>
class CompareByAttachment
{
	public:
		CompareByAttachment(Grid& g, TAttachment& aGID)
		{
			UG_ASSERT(g.has_attachment<TElem>(aGID),
					  "Can't compare unattached attachments!");
			m_aaGID.access(g, aGID);
		}

		bool operator()(TElem* e1, TElem* e2)
		{
			return m_aaGID[e1] < m_aaGID[e2];
		}

	private:
		Grid::AttachmentAccessor<TElem, TAttachment>	m_aaGID;
};

////////////////////////////////////////////////////////////////////////
//	SetAttachmentValues
///	sets attachment-values for elements between elemsBegin and elemsEnd.
template <class TAttachmentAccessor, class TIter, class TVal>
void SetAttachmentValues(TAttachmentAccessor& aaVal,
						TIter elemsBegin, TIter elemsEnd,
						const TVal& val);

////////////////////////////////////////////////////////////////////////
//	ConvertMathVectorAttachmentValues
///	Fills the dest-attachment with values from the source-attachment.
/**
 * TSrcAttachment and TDestAttachment have to have ValueTypes that
 * are compatible with ug::MathVector.
 *
 * Copies values from srcAttachment to destAttachment.
 * if destAttachment is not already attached, it will be attached
 * automatically. The srcAttachment however has to be attached.
 *
 * Valid types for TElem are: VertexBase, EdgeBase, Face, Volume
 *
 * If the dimensions do not match, the algorithm behaves as follows:
 * dim(src) > dim(dest): Only dim(dest) values are copied per element.
 * dim(src) < dim(dest): Values in dimensions >= dim(src) are set to 0.
 */
template<class TElem, class TSrcAttachment, class TDestAttachment>
bool ConvertMathVectorAttachmentValues(Grid& grid,
							TSrcAttachment& srcAttachment,
							TDestAttachment& destAttachment);


////////////////////////////////////////////////////////////////////////
///	copies attachments from one grid to the other
/**
 * If aSrc is not attached to srcGrid, false is returned.
 * If aDest is not attached to destGrid, it is attached automatically.
 *
 * The method iterates through the elements specified by TElem
 * and copies the attachments.
 *
 * Call like this: CopyAttachments<VertexBase>(...);
 */
template <class TElem, class TAttachment>
bool CopyAttachments(Grid& srcGrid, TAttachment& aSrc,
					Grid& destGrid, TAttachment& aDest);

////////////////////////////////////////////////////////////////////////
///	copies attachments for the specified elements
/**
 * If aSrc is not attached to srcGrid, false is returned.
 * If aDest is not attached to destGrid, it is attached automatically.
 */
template <class TElemIter, class TAttachment>
bool CopyAttachments(Grid& grid, TElemIter elemsBegin, TElemIter elemsEnd,
					 TAttachment& aSrc, TAttachment& aDest);

////////////////////////////////////////////////////////////////////////
///	assigns indices to the elements between begin and end.
/**	Indices are stored in the given attachment. Make sure that the
 *	given attachment-accessor operates on an attachment-pipe at which
 *	those elements are registered.
 */
template <class TIterator, class TAAInt>
void AssignIndices(TIterator begin, TIterator end,
					TAAInt& aaInt, int baseIndex = 0);

////////////////////////////////////////////////////////////////////////
///	returns the iterator whose element has the specified attachment value.
/** If no element contains the given value, end is returned.
 *
 * Make sure that the specified attachment accessor operates on the
 * specified elements.
 */
template <class TIterator, class TAttAcc>
TIterator FindElementByValue(TIterator begin, TIterator end,
							 const typename TAttAcc::ValueType& val,
							 TAttAcc& aa);

/**@}*/ // end of doxygen defgroup command
}//	end of namespace

////////////////////////////////////////////////
// include implementations of template methods.
#include "attachment_util_impl.hpp"

#endif
