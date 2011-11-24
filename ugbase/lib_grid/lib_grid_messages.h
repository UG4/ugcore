// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)

#ifndef __H__UG__lib_grid_messages__
#define __H__UG__lib_grid_messages__

#include "common/util/message_hub.h"

namespace ug
{

///	constants which indicate the adaption type in a grid refinement message.
enum GridMessageAdaptionType{
	GMAT_UNKNOWN = 0,
	GMAT_GLOBAL_REFINEMENT_BEGINS,
	GMAT_HNODE_REFINEMENT_BEGINS,
	GMAT_GLOBAL_REFINEMENT_ENDS,
	GMAT_HNODE_REFINEMENT_ENDS,
	GMAT_GLOBAL_COARSENING_BEGINS,
	GMAT_HNODE_COARSENING_BEGINS,
	GMAT_GLOBAL_COARSENING_ENDS,
	GMAT_HNODE_COARSENING_ENDS
};

///	A message sent along with "GridRefinement" messages.
class GridMessage_Adaption : public MessageHub::IMessage
{
	public:
		GridMessage_Adaption(GridMessageAdaptionType adaptionType) :
			m_adaptionType(adaptionType)	{}

		GridMessageAdaptionType refinement_type()		{return m_adaptionType;}

	protected:
		GridMessageAdaptionType	m_adaptionType;
};


///	Returns the message id for grid adaption messages.
/**	The associated type is GridMessage_Adaption.
 * The associated messageIdName is "GridAdaption".
 *
 * Note: This method shouldn't be called on a regular basis. Instead the returned
 * message-id should be cached and used directly in calls to the message-hub.
 * The id will be determined on the first call with a given hub and then stay the
 * same in subsequent calls with the same hub.
 *
 * Note: Message-ids are not necessarily the same for different hubs.
 */
inline int GridMessageId_Adaption(SPMessageHub hub)
{
	return hub->get_message_id<GridMessage_Adaption>("GridAdaption");
}

}//	end of namespace

#endif
