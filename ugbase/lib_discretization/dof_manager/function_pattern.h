/*
 * function_pattern.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__FUNCTION_PATTERN__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__FUNCTION_PATTERN__

#include <vector>
#include <string>

#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"

namespace ug{

class FunctionPattern
{
	public:
		FunctionPattern() : m_bLocked(false) {}

		/// add a single solution of LocalShapeFunctionSetID to the entire domain (i.e. all elements of the (Multi-)grid)
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 */
		virtual bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
		{
			// if already fixed, return false
			if(m_bLocked) {UG_LOG("Already fixed. Cannot change Distributor.\n"); return false;}

			// add to function list, everywhere = true, empty SubsetGroup
			m_vFunction.push_back(Function(name, dim, id, true, SubsetGroup()));

			return true;
		}

		/// add a single solution of LocalShapeFunctionSetID to selected subsets
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 * \param[in] SubsetIndices	SubsetGroup of subset indices, where this solution lives
		 */
		virtual bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, const SubsetGroup& SubsetIndices, int dim)
		{
			// if already fixed, return false
			if(m_bLocked) {UG_LOG("Already fixed. Cannot change Distributor.\n"); return false;}

			// add to function list, everywhere = false, SubsetGroup as given
			m_vFunction.push_back(Function(name, dim, id, false, SubsetIndices));

			return true;
		}

		bool lock()
		{
			m_bLocked = true;
			return true;
		}

		bool is_locked() const {return m_bLocked;}

		/// number of discrete functions this dof distributor handles
		size_t num_fct() const {return m_vFunction.size();}

		/// number of discrete functions this dof distributor handles an subset si
		size_t num_fct(int si) const
		{
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
				if(m_vFunction[fct].is_def_in_subset(si)) num++;
			}
			return num;
		}

		/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const  {return m_vFunction[fct].id;}

		/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {return m_vFunction[fct].name;}

		/// returns fct_id of the loc_fct on subset si
		size_t fct_id(size_t loc_fct, int si) const
		{
			size_t fct = 0;
			for(size_t i = 0; i < loc_fct; ++i)
			{
				if(is_def_in_subset(fct, si)) fct++;
			}
			return fct;
		}

		/// returns the dimension in which solution lives
		int dim(size_t fct) const {return m_vFunction[fct].dim;}

		/// returns true if the discrete function nr_fct is defined on subset s
		bool is_def_in_subset(size_t fct, int si) const {return m_vFunction[fct].is_def_in_subset(si);}

		// virtual destructor
		virtual ~FunctionPattern() {}

	protected:
		struct Function
		{
			Function(std::string name_, int dim_, LocalShapeFunctionSetID id_, bool everywhere_, const SubsetGroup& subsetIndices_)
				: name(name_), dim(dim_), id(id_), everywhere(everywhere_), subsetIndices(subsetIndices_){};

			std::string name;
			int dim;
			LocalShapeFunctionSetID id;
			bool everywhere;
			SubsetGroup subsetIndices;

			bool is_def_in_subset(int si) const
			{
				if(everywhere) return true;
				if(subsetIndices.containes_subset(si)) return true;
				return false;
			}
		};

	protected:
		// bool if pattern is fixed
		bool m_bLocked;

		// informations about Functions
		std::vector<Function> m_vFunction;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__FUNCTION_PATTERN__ */
