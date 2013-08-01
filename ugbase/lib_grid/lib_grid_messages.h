// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.11.2011 (m,d,y)

#ifndef __H__UG__lib_grid_messages__
#define __H__UG__lib_grid_messages__

#include "common/util/message_hub.h"
#include "lib_grid/grid/geometric_object_collection.h"
#include "lib_grid/algorithms/serialization.h"

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

		GridMessage_Adaption(GridMessageAdaptionType adaptionType,
							 const GeometricObjectCollection& affectedElements) :
			m_adaptionType(adaptionType),
			m_affectedElements(affectedElements){}

		GridMessageAdaptionType adaption_type() const	{return m_adaptionType;}

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

	///	returns the geometric-object-collection, which holds lists with affected elements.
	/**	- REFINEMENT_BEGINS: Elements on which new elements will occur or which are
	 * 						 transformed to a new element type (e.g. constrained/constraining)
	 *	- REFINEMENT_ENDS: Elements which were refined.
	 *	- COARSENING_BEGINS: Elements which will be removed or replaced during coarsening
	 *
	 *	For all other adaption types, the geometric-object-collection should be empty and
	 *	should be ignored.
	 *
	 *	\note	during refinement only elements on which a new child was actually created
	 *			are considered to be affected. Elements which were replaced by elements
	 *			of a different type are not considered in the list of affected elements.
	 *
	 *	\note	(THIS MAY CHANGE SOON!) Some of the elements marked for coarsening may
	 *			actually only be replaced by a normal/constraining/constrained element.
	 *			Note that elements may be marked during coarsening which actually have children
	 *			themselves, even if not all of those children are removed during coarsening
	 *			(e.g. edges which are replaced by constraining edges).*/
		const GeometricObjectCollection& affected_elements() const {return m_affectedElements;}

	protected:
		GridMessageAdaptionType		m_adaptionType;
		GeometricObjectCollection	m_affectedElements;
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
		GridMessage_Distribution(GridMessageDistributionType msg,
								 GridDataSerializationHandler& gdsh) :
			m_msg(msg),
			m_serializationHandler(gdsh)
		{}

		GridMessageDistributionType msg() const	{return m_msg;}

		GridDataSerializationHandler& serialization_handler() const
		{
			return m_serializationHandler;
		}

	private:
		GridMessageDistributionType		m_msg;
		GridDataSerializationHandler&	m_serializationHandler;
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
