// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 16.03.2011 (m,d,y)

#ifndef __H__UG__grid_callbacks__
#define __H__UG__grid_callbacks__

#include <iostream>
#include "lib_grid/grid/grid.h"
#include "common/serialization.h"

namespace ug
{
///	Serialization callback for grid attachments
/**	template class usable as CB_Serialize{ELEMENT}Data, where
 * {ELEMENT} is specified by TGeomObj and must be one of the
 * following types: VertexBase, EdgeBase, Face, Volume.
 *
 * Note that the attachment is automatically attached, if not yet present.
 */
/*
template <class TGeomObj, class TAttachment>
class GridAttachmentSerializer
{
	public:
		GridAttachmentSerializer(Grid& g, TAttachment& a) :
			m_aa(g, a, true)	{}

		void operator() (std::ostream& out, TGeomObj* o) const
		{Serialize(out, m_aa[o]);}

	private:
		Grid::AttachmentAccessor<TGeomObj, TAttachment>	m_aa;
};
*/
///	Deserialization callback for grid attachments
/**	template class usable as CB_Deserialize{ELEMENT}Data, where
 * {ELEMENT} is specified by TGeomObj and must be one of the
 * following types: VertexBase, EdgeBase, Face, Volume.
 *
 * Note that the attachment is automatically attached, if not yet present.
 */
/*
template <class TGeomObj, class TAttachment>
class GridAttachmentDeserializer
{
	public:
		GridAttachmentDeserializer(Grid& g, TAttachment& a) :
			m_aa(g, a, true)	{}

		void operator() (std::istream& in, TGeomObj* o)
		{Deserialize(in, m_aa[o]);}

	private:
		Grid::AttachmentAccessor<TGeomObj, TAttachment>	m_aa;
};
*/

}//	end of namespace

#endif
