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


///	Instances of this class inform about grid creation (e.g. during load or distribution)
enum GridMessageCreationType{
	GMCT_NONE,
	GMCT_CREATION_STARTS,
	GMCT_CREATION_STOPS
};

class GridMessage_Creation : public MessageHub::IMessage
{
	public:
		GridMessage_Creation(GridMessageCreationType msg = GMCT_NONE, int procId = -1) :
			m_msg(msg), m_procId(procId)	{}

		GridMessageCreationType msg() const	{return m_msg;}
		int proc_id() const					{return m_procId;}

	private:
		GridMessageCreationType	m_msg;
		int m_procId;
};

}//	end of namespace

#endif
