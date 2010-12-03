/*
 * dof_distribution.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__

#include <vector>

#include "./function_pattern.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_index_layout.h"
#endif

namespace ug{

// Base class for distributions
class DoFDistribution
{
	public:
		//DoFDistribution() : m_pFunctionPattern(NULL) {}

		DoFDistribution(FunctionPattern& dp) : m_pFunctionPattern(&dp) {}

		///////////////////////////
		// Infos
		///////////////////////////
		// TODO: The performance could be improved by caching the values instead of forwarding virtual functions

		/// number of discrete functions this dof distributor handles
		size_t num_fct() const {return m_pFunctionPattern->num_fct();}

		/// number of discrete functions this dof distributor handles an subset si
		size_t num_fct(int si) const {return m_pFunctionPattern->num_fct(si);}

		/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const  {return m_pFunctionPattern->local_shape_function_set_id(fct);}

		/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {return m_pFunctionPattern->name(fct);}

		/// returns the dimension in which solution lives
		int dim(size_t fct) const {return m_pFunctionPattern->dim(fct);}

		/// returns true if the discrete function nr_fct is defined on subset s
		bool is_def_in_subset(size_t fct, int si) const {return m_pFunctionPattern->is_def_in_subset(fct, si);}

		/// returns function id for local function on subset
		size_t fct_id(size_t loc_fct, int si) const {return m_pFunctionPattern->fct_id(loc_fct, si);}

	protected:
		// DoFDistributor for dofs
		FunctionPattern* m_pFunctionPattern;

#ifdef UG_PARALLEL
	public:
		IndexLayout& get_slave_layout()	{return m_slaveLayout;}
		IndexLayout& get_master_layout()	{return m_masterLayout;}
		IndexLayout& get_vertical_slave_layout()		{return m_verticalSlaveLayout;}
		IndexLayout& get_vertical_master_layout()	{return m_verticalMasterLayout;}

		pcl::ParallelCommunicator<IndexLayout>& get_communicator()	{return m_communicator;}
		pcl::ProcessCommunicator& get_process_communicator()	{return m_processCommunicator;}

		size_t num_master_dofs() {return num_dofs(m_masterLayout);}
		size_t num_slave_dofs() {return num_dofs(m_slaveLayout);}
		size_t num_vertical_master_dofs() {return num_dofs(m_verticalMasterLayout);}
		size_t num_vertical_slave_dofs() {return num_dofs(m_verticalSlaveLayout);}

	protected:
		size_t num_dofs(IndexLayout& Layout) const
		{
			size_t sum = 0;
			for(IndexLayout::iterator iter = Layout.begin();
					iter != Layout.end(); ++iter)
				sum += Layout.interface(iter).size();
			return sum;
		}

	protected:
		// index layout for each grid level
		IndexLayout m_slaveLayout;

		// index layout for each grid level
		IndexLayout m_masterLayout;

		// index layout for each grid level
		IndexLayout m_verticalMasterLayout;

		// index layout for each grid level
		IndexLayout m_verticalSlaveLayout;

		// process communicator
		pcl::ProcessCommunicator m_processCommunicator;

		// communicator
		pcl::ParallelCommunicator<IndexLayout> m_communicator;
#endif
};



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__ */
