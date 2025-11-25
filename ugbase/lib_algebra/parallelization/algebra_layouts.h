/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
		HorizontalAlgebraLayouts() : m_overlapEnabled(false) 	{}

	///	clears the struct
		void clear()
		{
			masterLayout.clear();
			slaveLayout.clear();
		}

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		[[nodiscard]] const IndexLayout& master() const 			{return masterLayout;}
		[[nodiscard]] const IndexLayout& master_overlap() const 	{return masterOverlapLayout;}
		[[nodiscard]] const IndexLayout& slave() const 				{return slaveLayout;}
		[[nodiscard]] const IndexLayout& slave_overlap() const		{return slaveOverlapLayout;}
	/// \}

	///	returns process communicator
		[[nodiscard]] const pcl::ProcessCommunicator& proc_comm() const 	{return processCommunicator;}

	///	returns (non-const !!!) communicator
	/**
	 * We return a non-const communicator, since it can only be used non-const.
	 * This is requirement of the InterfaceCommunicator design, since internally
	 * buffers are used during communication. Theoretically, a communicator could
	 * be created each time for communication, but for performance reasons we
	 * provide this one that can be shared and reused.
	 */
		[[nodiscard]] pcl::InterfaceCommunicator<IndexLayout>& comm() const {
			return const_cast<HorizontalAlgebraLayouts*>(this)->communicator;
		}

	/**	It is important to enable or disable overlap on all involved processes
	 * at the same time. Otherwise communication issues may arise.*/
		void enable_overlap(bool enable)	{m_overlapEnabled = enable;}
	///	Tells whether overlap interfaces should be considered
		[[nodiscard]] bool overlap_enabled() const {return m_overlapEnabled;}

	public:
	/// returns the horizontal slave/master index layout
	/// \{
		IndexLayout& master()			{return masterLayout;}
		IndexLayout& master_overlap() 	{return masterOverlapLayout;}
		IndexLayout& slave()			{return slaveLayout;}
		IndexLayout& slave_overlap() 	{return slaveOverlapLayout;}
	/// \}

	///	returns communicator
	/// \{
		pcl::InterfaceCommunicator<IndexLayout>& comm() 	{return communicator;}
		pcl::ProcessCommunicator& proc_comm()				{return processCommunicator;}
	/// \}

	protected:
		///	(horizontal) master index layout
		IndexLayout masterLayout;
		IndexLayout masterOverlapLayout;

		///	(horizontal) slave index layout
		IndexLayout slaveLayout;
		IndexLayout slaveOverlapLayout;

		///	process communicator
		pcl::ProcessCommunicator processCommunicator;

		///	communicator
		pcl::InterfaceCommunicator<IndexLayout> communicator;

		bool m_overlapEnabled;
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
		[[nodiscard]] const IndexLayout& vertical_master() const 	{return verticalMasterLayout;}
		[[nodiscard]] const IndexLayout& vertical_slave() const 	{return verticalSlaveLayout;}
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

std::ostream &operator << (std::ostream &out, const HorizontalAlgebraLayouts &layouts);
std::ostream &operator << (std::ostream &out, const AlgebraLayouts &layouts);

#endif


} // end namespace ug

#endif