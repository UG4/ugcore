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

#include "common/common.h"
#include "lib_discretization/domain_util.h"
#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_id.h"
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	ERROR_FunctionPatternLocked
struct ERROR_FunctionPatternLocked{};

/// Describes the setup of discrete functions on a SubsetHandler
/**
 * A Function Pattern is used to describe the setup for discrete functions on
 * a domain. Therefore, it has the underlying SubsetHandler of the Domain.
 * The functions (sometimes also called grid functions) are defined with respect
 * to the subsets: A function can 'live' on parts of the subsets as well as on
 * the whole domain.
 */
class FunctionPattern
{
	public:
	///	Default Constructor
		FunctionPattern() :
			m_bLocked(false), m_pSH(NULL)
		{clear();}

	/// set an underlying subset handler
		void set_subset_handler(const ISubsetHandler& sh) {m_pSH = &sh; clear();}

	/// get underlying subset handler
		const ISubsetHandler* get_subset_handler() const {return m_pSH;}

	//	returns if ansatz space is supported
		virtual bool supports_trial_space(LocalShapeFunctionSetID& id) const = 0;

	/// add a single solution of LocalShapeFunctionSetID to the entire domain
	/**
	 * \param[in] 	name		Name of this Single Solution
	 * \param[in] 	id			Shape Function set id
	 * \param[in]	dim			Dimension (optional)
	 */
		virtual bool add_fct(const char* name,
		                     LocalShapeFunctionSetID id,
		                     int dim = -1);

	/// add a single solution of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] name			Name of this Single Solution
	 * \param[in] id			Shape Function set id
	 * \param[in] SubsetIndices	SubsetGroup, where solution lives
	 * \param[in] dim			Dimension
	 */
		virtual bool add_fct(const char* name,
		                     LocalShapeFunctionSetID id,
		                     const SubsetGroup& SubsetIndices,
		                     int dim = -1);

	/// add a single solution of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] name			Name of this Single Solution
	 * \param[in] id			Shape Function set id
	 * \param[in] subsets		Subsets separated by ','
	 * \param[in] dim			Dimension
	 */
		virtual bool add_fct(const char* name,
		                     LocalShapeFunctionSetID id,
		                     const char* subsets,
		                     int dim = -1);

		virtual bool add_fct(const char* name,
		                     const char* type,
		                     int order)
		{
			std::string strType(type);
			if(strType != "Lagrange" || order != 1) return false;

			return add_fct(name, LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE,1));
		}

		virtual bool add_fct_on_subset(const char* name,
		                               const char* type,
		                               int order,
		                               const char* subsets)
		{
			std::string strType(type);
			if(strType != "Lagrange" || order != 1) return false;

			return add_fct(name, LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE,1), subsets);
		}

	///	lock pattern (i.e. can not be changed then)
		inline void lock()	{m_bLocked = true;}

	///	returns if pattern is locked
		inline bool is_locked() const {return m_bLocked;}

	/// clear all functions
		void clear()
		{
			if(is_locked()) throw(ERROR_FunctionPatternLocked());
			m_vFunction.clear();
		}

	///	number of subsets
		int num_subsets() const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			return m_pSH->num_subsets();
		}

	///	dimension of subset
		int dim_subset(int si) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			return DimensionOfSubset(*m_pSH, si);
		}

	/// number of discrete functions this dof distributor handles
		size_t num_fct() const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			return m_vFunction.size();
		}

	/// number of discrete functions on a subset
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

	/// returns the trial space of a discrete function
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].id;
		}

	/// returns the name of a discrete function
		const char* name(size_t fct) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].name.c_str();
		}

	/// returns function id of the loc_fct's function on subset si
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

	/// returns true if the discrete function nr_fct is defined on subset si
		bool is_def_in_subset(size_t fct, int si) const
		{
			UG_ASSERT(m_pSH != NULL, "SubsetHandler not set.");
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].is_def_in_subset(si);
		}

	/// virtual destructor
		virtual ~FunctionPattern() {}

	protected:
	///	internal structure to hold all functions
		struct Function
		{
			Function(const char* name_, int dim_,
			         LocalShapeFunctionSetID id_, bool everywhere_,
			         const SubsetGroup& subsetIndices_)
				: name(name_), dim(dim_), id(id_),
				  everywhere(everywhere_), subsetIndices(subsetIndices_){};

			std::string name;
			int dim;
			LocalShapeFunctionSetID id;
			bool everywhere;
			SubsetGroup subsetIndices;

			bool is_def_in_subset(int si) const
			{
				if(everywhere) return true;
				if(subsetIndices.contains(si)) return true;
				return false;
			}
		};

	protected:
	// 	flag if pattern is fixed
		bool m_bLocked;

	// 	underlying subset handler
		const ISubsetHandler* m_pSH;

	// 	informations about Functions
		std::vector<Function> m_vFunction;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__FUNCTION_PATTERN__ */
