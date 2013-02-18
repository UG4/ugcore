/*
 * algebra_layouts.h
 *
 *  Created on: 18.02.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG4__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_LAYOUTS__
#define __H__UG4__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_LAYOUTS__

#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#endif

namespace ug{

#ifdef UG_PARALLEL
struct AlgebraLayouts
{
	public:
	///	clears the struct
		void clear()
		{
			masterLayout.clear();			slaveLayout.clear();
			verticalMasterLayout.clear();	verticalSlaveLayout.clear();
		}

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		const IndexLayout& master() const 				{return masterLayout;}
		const IndexLayout& slave() const 				{return slaveLayout;}
	/// \}

	/// returns the vertical slave/master index layout
	/// \{
		const IndexLayout& vertical_master() const 	{return verticalMasterLayout;}
		const IndexLayout& vertical_slave() const 	{return verticalSlaveLayout;}
	/// \}

	///	returns communicator
	/// \{
		const pcl::InterfaceCommunicator<IndexLayout>& comm() const  	{return communicator;}
		const pcl::ProcessCommunicator& proc_comm() const 				{return processCommunicator;}
	/// \}

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		IndexLayout& master()	{return masterLayout;}
		IndexLayout& slave()	{return slaveLayout;}
	/// \}

	/// returns the vertical slave/master index layout
	/// \{
		IndexLayout& vertical_master() 		{return verticalMasterLayout;}
		IndexLayout& vertical_slave()  		{return verticalSlaveLayout;}
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

		///	vertical master index layout
		IndexLayout verticalMasterLayout;

		///	vertical slave index layout
		IndexLayout verticalSlaveLayout;

		///	process communicator
		pcl::ProcessCommunicator processCommunicator;

		///	communicator
		pcl::InterfaceCommunicator<IndexLayout> communicator;
};
#endif

} // end namespace ug

#endif /* __H__UG4__LIB_ALGEBRA__PARALLELIZATION__ALGEBRA_LAYOUTS__ */
