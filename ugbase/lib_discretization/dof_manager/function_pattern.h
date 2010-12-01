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
#include "lib_discretization/local_shape_function_set/local_shape_function_set_id.h"
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	ERROR_FunctionPatternLocked
struct ERROR_FunctionPatternLocked{};


class FunctionPattern
{
	public:
		FunctionPattern() : m_bLocked(false), m_pSH(NULL) {clear();}

		/// set an underlying subset handler
		void set_subset_handler(const ISubsetHandler& sh) {m_pSH = &sh; clear();}

		/// get underlying subset handler
		const ISubsetHandler* get_subset_handler() const {return m_pSH;}

		/// add a single solution of LocalShapeFunctionSetID to the entire domain (i.e. all elements of the (Multi-)grid)
		/**
		 * \param[in] 	name		Name of this Single Solution
		 * \param[in] 	id			Shape Function set id
		 * \param[in]	dim			Dimension
		 */
		virtual bool add_discrete_function(const char* name,
		                                   LocalShapeFunctionSetID id, int dim)
		{
		// 	if already fixed, return false
			if(m_bLocked)
			{
				UG_LOG("Already fixed. Cannot change Distributor.\n");
				return false;
			}

		//	check that subset handler exists
			if(m_pSH == NULL)
			{
				UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
						"SubsetHandler not set.\n");
				return false;
			}

		//	create temporary subset group
			SubsetGroup tmpSubsetGroup;
			tmpSubsetGroup.set_subset_handler(*m_pSH);
			tmpSubsetGroup.add_all_subsets();

		// 	add to function list, everywhere = true, copy SubsetGroup
			m_vFunction.push_back(Function(name, dim, id, true, tmpSubsetGroup));

			return true;
		}

		/// add a single solution of LocalShapeFunctionSetID to selected subsets
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] id			Shape Function set id
		 * \param[in] SubsetIndices	SubsetGroup of subset indices, where this solution lives
		 * \param[in] dim			Dimension
		 */
		virtual bool add_discrete_function(const char* name,
		                                   LocalShapeFunctionSetID id,
		                                   const SubsetGroup& SubsetIndices,
		                                   int dim)
		{
		// 	if already fixed, return false
			if(m_bLocked) {UG_LOG("Already fixed. Cannot change Distributor.\n"); return false;}

		//	check that subset handler exists
			if(m_pSH == NULL)
			{
				UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
						"SubsetHandler not set.\n");
				return false;
			}

		//	check that subset handler are equal
			if(m_pSH != SubsetIndices.get_subset_handler())
			{
				UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
						"SubsetHandler of SubsetGroup does not match SubsetHandler of FunctionPattern.\n");
				return false;
			}

		// 	add to function list, everywhere = false, copy SubsetGroup as given
			m_vFunction.push_back(Function(name, dim, id, false, SubsetIndices));

			return true;
		}

		inline void lock()	{m_bLocked = true;}

		inline bool is_locked() const {return m_bLocked;}

		/// clear all functions
		void clear() {if(is_locked()) throw(ERROR_FunctionPatternLocked()); m_vFunction.clear();}

		/// number of discrete functions this dof distributor handles
		size_t num_fct() const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			return m_vFunction.size();
		}

		/// number of discrete functions this dof distributor handles an subset si
		size_t num_fct(int si) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
				if(m_vFunction[fct].is_def_in_subset(si)) num++;
			}
			return num;
		}

		/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].id;
		}

		/// returns the name of the discrete function nr_fct
		const char* name(size_t fct) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].name.c_str();
		}

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

		/// returns the function id if function with given name found in pattern, -1 else
		size_t fct_id_by_name(const char* name) const
		{
			for(size_t i = 0; i < m_vFunction.size(); ++i)
			{
				if(m_vFunction[i].name == name)
					return i;
			}

			return (size_t) -1;
		}

		/// returns the dimension in which solution lives
		int dim(size_t fct) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].dim;
		}

		/// returns true if the discrete function nr_fct is defined on subset s
		bool is_def_in_subset(size_t fct, int si) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].is_def_in_subset(si);
		}

		// virtual destructor
		virtual ~FunctionPattern() {}

	protected:
		struct Function
		{
			Function(const char* name_, int dim_, LocalShapeFunctionSetID id_, bool everywhere_, const SubsetGroup& subsetIndices_)
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

		// underlying subset handler
		const ISubsetHandler* m_pSH;

		// informations about Functions
		std::vector<Function> m_vFunction;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__FUNCTION_PATTERN__ */
