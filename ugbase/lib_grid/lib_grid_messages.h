// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)

#ifndef __H__UG__lib_grid_messages__
#define __H__UG__lib_grid_messages__

#include "common/util/message_hub.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
enum GridMessageMultiGridChangedType{
	GMMGCT_LEVEL_ADDED = 1,
	GMMGCT_REMOVED = 2
};

///	A message sent by the MultiGrid, if something special happened.
class GridMessage_MultiGridChanged : public MessageHub::IMessage
{
	public:
		GridMessage_MultiGridChanged(GridMessageMultiGridChangedType msgType,
									 int numLevels) :
			m_msgType(msgType),
			m_numLevels(numLevels)
		{}

	///	returns the type of the message
		GridMessageMultiGridChangedType message_type() const	{return m_msgType;}

	///	returns the current number of levels in the multigrid, which did send the message.
		int num_levels_in_grid() const		{return m_numLevels;}

	private:
		GridMessageMultiGridChangedType	m_msgType;
		int		m_numLevels;
};

///	Returns the message id for grid adaption messages.
/**	The associated type is GridMessage_MultiGridChanged.
 * The associated messageIdName is "MultiGridChanged".
 *
 * Note: Message-ids are not necessarily the same for different hubs.
 */
inline int GridMessageId_MultiGridChanged(SPMessageHub hub)
{
	int id = hub->get_message_id<GridMessage_MultiGridChanged>("MultiGridChanged");
	return id;
}


////////////////////////////////////////////////////////////////////////////////
///	constants which indicate the adaption type in a grid refinement message.
enum GridMessageAdaptionType{
	GMAT_UNKNOWN = 0,
	GMAT_GLOBAL_ADAPTION_BEGINS,
	GMAT_HNODE_ADAPTION_BEGINS,
	GMAT_GLOBAL_ADAPTION_ENDS,
	GMAT_HNODE_ADAPTION_ENDS,
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

		GridMessageAdaptionType refinement_type() const	{return m_adaptionType;}

	///	tells whether grid adaption has just been started or has been finished.
	/**	Note that begins and ends may both be false. An adaption consists of
	 * multiple steps.
	 * \{ */
		bool adaption_begins() const;
		bool adaption_ends() const;
	/** \} */

	///	tells whether an adaption step has just been started or has been finished.
	/**	Note that both may be false. An adaption consists of multiple adaption steps.
	 * \{ */
		bool step_begins() const;
		bool step_ends() const;
	/** \} */

	///	tells whether refinement / coarsening is adaptive (adaptive() == !global())
		bool adaptive() const;

	///	tells whether refinement / coarsening is global (global() == !adaptive())
		bool global() const;

	///	tells whether a step is a refinement step.
	/**	This is always false if adaption_begins() or adaption_ends() returns true*/
		bool refinement() const;

	///	tells whether a step is a coarsen step
	/**	This is always false if adaption_begins() or adaption_ends() returns true*/
		bool coarsening() const;

	protected:
		GridMessageAdaptionType	m_adaptionType;
};


///	Returns the message id for grid adaption messages.
/**	The associated type is GridMessage_Adaption.
 * The associated messageIdName is "GridAdaption".
 *
 * Note: Message-ids are not necessarily the same for different hubs.
 */
inline int GridMessageId_Adaption(SPMessageHub hub)
{
	static int id = hub->get_message_id<GridMessage_Adaption>("GridAdaption");
	return id;
}


///	Instances of this class inform about distribution
enum GridMessageDistributionType{
	GMDT_NONE,
	GMDT_DISTRIBUTION_STARTS,
	GMDT_DISTRIBUTION_STOPS
};

class GridMessage_Distribution : public MessageHub::IMessage
{
	public:
		GridMessage_Distribution(GridMessageDistributionType msg = GMDT_NONE) :
			m_msg(msg)	{}

		GridMessageDistributionType msg() const	{return m_msg;}

	private:
		GridMessageDistributionType	m_msg;
};

///	Returns the message id for grid distribution messages.
/**	The associated type is GridMessage_Distribution.
 * The associated messageIdName is "GridDistribution".
 *
 * Note: Message-ids are not necessarily the same for different hubs.
 */
inline int GridMessageId_Distribution(SPMessageHub hub)
{
	static int id = hub->get_message_id<GridMessage_Distribution>("GridDistribution");
	return id;
}


}//	end of namespace

#endif
