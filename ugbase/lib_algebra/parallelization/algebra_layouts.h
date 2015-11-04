
#ifndef __H__UG4__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_LAYOUTS__
#define __H__UG4__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_LAYOUTS__

#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#endif

namespace ug{

#ifdef UG_PARALLEL

///	Holds Interfaces and communicators for horizontal communication.
/** The internal ProcessCommunicator is initialized as PCD_WORLD.*/
class HorizontalAlgebraLayouts
{
	public:
	///	clears the struct
		void clear()
		{
			masterLayout.clear();			slaveLayout.clear();
		}

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		const IndexLayout& master() const 				{return masterLayout;}
		const IndexLayout& slave() const 				{return slaveLayout;}
	/// \}

	///	returns process communicator
		const pcl::ProcessCommunicator& proc_comm() const 				{return processCommunicator;}

	///	returns (non-const !!!) communicator
	/**
	 * We return a non-const communicator, since it can only be used non-const.
	 * This is requirement of the InterfaceCommunicator design, since internally
	 * buffers are used during communication. Theoretically, a communicator could
	 * be created each time for communication, but for performance reasons we
	 * provide this one that can be shared and reused.
	 */
		pcl::InterfaceCommunicator<IndexLayout>& comm() const  	{return const_cast<HorizontalAlgebraLayouts*>(this)->communicator;}

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		IndexLayout& master()	{return masterLayout;}
		IndexLayout& slave()	{return slaveLayout;}
	/// \}

	///	returns communicator
	/// \{
		pcl::InterfaceCommunicator<IndexLayout>& comm() 	{return communicator;}
		pcl::ProcessCommunicator& proc_comm()			{return processCommunicator;}
	/// \}

	protected:
		///	(horizontal) master index layout
		IndexLayout masterLayout;

		///	(horizontal) slave index layout
		IndexLayout slaveLayout;

		///	process communicator
		pcl::ProcessCommunicator processCommunicator;

		///	communicator
		pcl::InterfaceCommunicator<IndexLayout> communicator;
};

///	Extends the HorizontalAlgebraLayouts by vertical layouts.
class AlgebraLayouts : public HorizontalAlgebraLayouts
{
	public:
	///	clears the struct
		void clear()
		{
			HorizontalAlgebraLayouts::clear();
			verticalMasterLayout.clear();	verticalSlaveLayout.clear();
		}

	public:
	/// returns the vertical slave/master index layout
	/// \{
		const IndexLayout& vertical_master() const 	{return verticalMasterLayout;}
		const IndexLayout& vertical_slave() const 	{return verticalSlaveLayout;}
	/// \}

	public:
	/// returns the vertical slave/master index layout
	/// \{
		IndexLayout& vertical_master() 		{return verticalMasterLayout;}
		IndexLayout& vertical_slave()  		{return verticalSlaveLayout;}
	/// \}

	protected:
		///	vertical master index layout
		IndexLayout verticalMasterLayout;

		///	vertical slave index layout
		IndexLayout verticalSlaveLayout;
};
#endif

std::ostream &operator << (std::ostream &out, const HorizontalAlgebraLayouts &layouts);
std::ostream &operator << (std::ostream &out, const AlgebraLayouts &layouts);


} // end namespace ug

#endif /* __H__UG4__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_LAYOUTS__ */
