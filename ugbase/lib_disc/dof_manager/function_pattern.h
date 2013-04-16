/*
 * function_pattern.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__FUNCTION_PATTERN__
#define __H__UG__LIB_DISC__DOF_MANAGER__FUNCTION_PATTERN__

#include <vector>
#include <string>

#include "common/common.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_disc/common/subset_util.h"
#include "lib_disc/common/subset_group.h"
#include "common/util/string_util.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"

namespace ug{

// predeclaration
class FunctionGroup;

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
		FunctionPattern(ConstSmartPtr<ISubsetHandler> spSH) :
			m_bLocked(false), m_spSH(spSH)
		{clear();}

	///	sets new subsets handler
		void set_subset_handler(ConstSmartPtr<ISubsetHandler> spSH);

	/// get underlying subset handler
		ConstSmartPtr<ISubsetHandler> subset_handler() const {return m_spSH;}

	/// add single solutions of LocalShapeFunctionSetID to the entire domain
	/**
	 * \param[in] 	name		name(s) of single solution (comma separated)
	 * \param[in] 	id			Shape Function set id
	 * \param[in]	dim			Dimension (optional)
	 */
		void add(const std::vector<std::string>& vName, LFEID id, int dim = -1);

	/// add single solutions of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] name			name(s) of single solution (comma separated)
	 * \param[in] id			Shape Function set id
	 * \param[in] subsets		Subsets separated by ','
	 * \param[in] dim			Dimension
	 */
		void add(const std::vector<std::string>& vName, LFEID id,
				 const std::vector<std::string>& vSubset, int dim = -1);

	protected:
	/// add single solutions of LocalShapeFunctionSetID to selected subsets
	/**
	 * \param[in] name			name(s) of single solution (comma separated)
	 * \param[in] id			Shape Function set id
	 * \param[in] SubsetIndices	SubsetGroup, where solution lives
	 * \param[in] dim			Dimension
	 */
		void add(const std::vector<std::string>& vName, LFEID id, const SubsetGroup& ssGrp, int dim = -1);

	public:
	///	lock pattern (i.e. can not be changed then)
		void lock()	{m_bLocked = true;}

	///	returns if pattern is locked
		bool is_locked() const {return m_bLocked;}

	/// clear all functions
		void clear()
		{
			if(is_locked()) UG_THROW("Pattern locked.");
			m_vFunction.clear();
		}

	///	number of subsets
		int num_subsets() const {return m_spSH->num_subsets();}

		SubsetGroup subset_grp_by_name(const char* names) const
		{
			return SubsetGroup(subset_handler(), TokenizeString(names));
		}

	///	dimension of subset
		int dim_subset(int si) const {return DimensionOfSubset(*m_spSH, si);}

	///	returns the name of a subset
		const char* subset_name(int si) const {return m_spSH->subset_info(si).name.c_str();}

	/// number of discrete functions this dof distributor handles
		size_t num_fct() const {return m_vFunction.size();}

	/// number of discrete functions on a subset
		size_t num_fct(int si) const
		{
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
				if(m_vFunction[fct].is_def_in_subset(si)) num++;
			}
			return num;
		}

	/// returns the trial space of a discrete function
	/// \{
		const LFEID& local_finite_element_id(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].id;
		}
		const LFEID& lfeid(size_t fct) const {return local_finite_element_id(fct);}
	/// \}

	/// returns the name of a discrete function
		const char* name(size_t fct) const
		{
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

	/// returns the function id if function with given name found in pattern, exception else
		size_t fct_id_by_name(const char* name) const;

	/// returns the dimension in which solution lives
		int dim(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].dim;
		}

	/// returns true if the discrete function nr_fct is defined on subset si
		bool is_def_in_subset(size_t fct, int si) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].is_def_in_subset(si);
		}

	/// returns true if the discrete function nr_fct is defined on all subsets
		bool is_def_everywhere(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Invalid index.");
			return m_vFunction[fct].is_def_everywhere();
		}

	/// virtual destructor
		virtual ~FunctionPattern() {}

	protected:
	///	internal structure to hold all functions
		struct Function
		{
			Function(const char* name_, int dim_,
			         LFEID id_, bool everywhere_,
			         const SubsetGroup& subsetIndices_)
				: name(name_), dim(dim_), id(id_),
				  everywhere(everywhere_), subsetIndices(subsetIndices_){};

			std::string name;
			int dim;
			LFEID id;
			bool everywhere;
			SubsetGroup subsetIndices;

			bool is_def_in_subset(int si) const
			{
				if(everywhere) return true;
				if(subsetIndices.contains(si)) return true;
				return false;
			}

			bool is_def_everywhere() const {return everywhere;}
		};

	protected:
	// 	flag if pattern is fixed
		bool m_bLocked;

	// 	underlying subset handler
		ConstSmartPtr<ISubsetHandler> m_spSH;

	// 	informations about Functions
		std::vector<Function> m_vFunction;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__FUNCTION_PATTERN__ */
