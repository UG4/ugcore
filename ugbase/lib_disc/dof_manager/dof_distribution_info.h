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

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION_INFO__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION_INFO__

#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/dof_manager/function_pattern.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/local_finite_element/local_dof_set.h"

#include <vector>

namespace ug{

class DoFDistributionInfo : public FunctionPattern
{
	public:
		///	indication that function is not defined on a subset
		enum{NOT_DEF_ON_SUBSET = static_cast<size_t>(-1)};
		///	indication that function is not defined on a subset
		enum{NOT_YET_ASSIGNED = static_cast<size_t>(-2)};
		///	indicate not set value
		enum{NOT_SPECIFIED = CommonLocalDoFSet::NOT_SPECIFIED};

	public:
		/// constructor
		explicit DoFDistributionInfo(ConstSmartPtr<ISubsetHandler> spSH);

		///	initializes the DoFs
		void init();


		/// returns the maximum number of dofs on grid objects in a dimension
		[[nodiscard]] size_t max_dofs(const int dim) const {return m_vMaxDoFsInDim[dim];}

		/// returns the maximum number of dofs on a grid base object type
		[[nodiscard]] size_t max_dofs(const GridBaseObjectId gbo) const {return m_vMaxDoFsInDim[gbo];}

		/// returns the maximum number of dofs on reference object type
		[[nodiscard]] size_t max_dofs(const ReferenceObjectID roid) const {return m_vMaxDoFsOnROID[roid];}

		/// returns the maximum number of dofs in a dimension on a subset
		[[nodiscard]] size_t max_dofs(const int dim, const int si) const {return m_vvMaxDoFsInDimPerSubset[dim][si];}

		/// returns the maximum number of dofs on a grid base object on a subset
		[[nodiscard]] size_t max_dofs(const GridBaseObjectId gbo, const int si) const {return m_vvMaxDoFsInDimPerSubset[gbo][si];}

		///	returns the number of dofs on a Reference Object on a subset
		[[nodiscard]] size_t num_dofs(const ReferenceObjectID roid, const int si) const {return m_vvNumDoFsOnROIDPerSubset[roid][si];}


		/// returns the maximal number of dofs on a dimension for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, int dim) const {return m_vFctInfo[fct].vMaxDoFsInDim[dim];}

		/// returns the maximal number of dofs on a base object type for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const GridBaseObjectId gbo) const {return m_vFctInfo[fct].vMaxDoFsInDim[gbo];}

		///	returns the number of dofs on a reference object for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const ReferenceObjectID roid) const {return m_vFctInfo[fct].vMaxDoFsOnROID[roid];}

		/// returns the maximum number of dofs in a dimension on a subset for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const int dim, const int si) const {return m_vFctInfo[fct].vvMaxDoFsInDimPerSubset[dim][si];}

		/// returns the maximum number of dofs on a grid base object on a subset for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const GridBaseObjectId gbo, const int si) const {return m_vFctInfo[fct].vvMaxDoFsInDimPerSubset[gbo][si];}

		///	returns the number of dofs on a Reference Object on a subset for a function component
		[[nodiscard]] size_t num_fct_dofs(const size_t fct, const ReferenceObjectID roid, const int si) const {return m_vFctInfo[fct].vvNumDoFsOnROIDPerSubset[roid][si];}


		///	returns the offset for reference element, subset and function
		[[nodiscard]] size_t offset(const ReferenceObjectID roid, const int si, const size_t fct) const {return m_vFctInfo[fct].vvOffsets[roid][si];}


		///	prints informations
		void print_local_dof_statistic(int verboseLev) const;

		///	prints statistic on local dof distribution
		void print_local_dof_statistic() const {print_local_dof_statistic(1);}

	protected:
		/// creates offset arrays
		void create_offsets();

	protected:
		///	maximum number of DoFs on geometric objects in a dimension
		size_t m_vMaxDoFsInDim[NUM_GEOMETRIC_BASE_OBJECTS];

		///	maximum number of DoFs on a reference type
		size_t m_vMaxDoFsOnROID[NUM_REFERENCE_OBJECTS];

		///	maximum number of DoFs on geometric objects in a dimension per subset
		std::vector<size_t> m_vvMaxDoFsInDimPerSubset[NUM_GEOMETRIC_BASE_OBJECTS];

		///	number of DoFs on a reference element type on a subset
		std::vector<size_t> m_vvNumDoFsOnROIDPerSubset[NUM_REFERENCE_OBJECTS];


		struct FctInfo{
			///	number Dofs for local DoF set and subelement of element
			size_t vMaxDoFsInDim[NUM_GEOMETRIC_BASE_OBJECTS];

			///	number Dofs for local DoF set and subelement of element
			size_t vMaxDoFsOnROID[NUM_REFERENCE_OBJECTS];

			///	maximum number of DoFs on geometric objects in a dimension per subset
			std::vector<size_t> vvMaxDoFsInDimPerSubset[NUM_GEOMETRIC_BASE_OBJECTS];

			///	number of DoFs on a reference element type on a subset
			std::vector<size_t> vvNumDoFsOnROIDPerSubset[NUM_REFERENCE_OBJECTS];

			/// offset map
			std::vector<size_t> vvOffsets[NUM_REFERENCE_OBJECTS];
		};

		/// infos for a function component
		std::vector<FctInfo> m_vFctInfo;
};


class DoFDistributionInfoProvider{
	public:
		/// constructor
		explicit DoFDistributionInfoProvider(ConstSmartPtr<DoFDistributionInfo> spDDI)
			: m_spDDI(spDDI)
		{}

		/// constructor
		DoFDistributionInfoProvider() : m_spDDI(nullptr) {}

		/// sets the dd info
		void set_dof_distribution_info(ConstSmartPtr<DoFDistributionInfo> spDDI) {m_spDDI = spDDI;}

		///	returns underlying info
		[[nodiscard]] ConstSmartPtr<DoFDistributionInfo> dof_distribution_info() const {return m_spDDI;}

		///	returns the subset handler
		[[nodiscard]] ConstSmartPtr<ISubsetHandler> subset_handler() const {return m_spDDI->subset_handler();}

		///	returns the function pattern
		[[nodiscard]] ConstSmartPtr<FunctionPattern> function_pattern() const {return m_spDDI;}


		/// number of discrete functions on subset si
		[[nodiscard]] size_t num_fct() const {return m_spDDI->num_fct();}

		/// number of discrete functions on subset si
		[[nodiscard]] size_t num_fct(int si) const {return m_spDDI->num_fct(si);}

		/// returns the name of the discrete function nr_fct
		[[nodiscard]] std::string name(size_t fct) const {return m_spDDI->name(fct);}

		/// returns the names of the discrete functions
		[[nodiscard]] std::vector<std::string> names() const {return m_spDDI->names();}

		/// returns fct id by name
		size_t fct_id_by_name(const char* name) const{return m_spDDI->fct_id_by_name(name);}


		///	returns number of subsets
		[[nodiscard]] int num_subsets() const {return m_spDDI->num_subsets();}

		/// returns the dimension in which solution lives
		[[nodiscard]] int dim(size_t fct) const {return m_spDDI->dim(fct);}

		///	returns dimension of subset
		[[nodiscard]] int dim_subset(int si) const {return m_spDDI->dim_subset(si);}

		///	returns subset name
		[[nodiscard]] std::string subset_name(int si) const {return m_spDDI->subset_name(si);}

		/// returns the subset id
		int subset_id_by_name(const char* name) const {return m_spDDI->subset_id_by_name(name);}

		///	returns subset group by name
		SubsetGroup subset_grp_by_name(const char* names) const;

		///	returns subset group of all the subsets except for (by name)
		SubsetGroup all_subsets_grp_except_for(const char* names) const;

		///	returns a function group to a string of functions
		FunctionGroup fct_grp_by_name(const char* names) const;


		///	returns if a function is defined on a subset
		[[nodiscard]] bool is_def_in_subset(size_t fct, int si) const {return m_spDDI->is_def_in_subset(fct,si);}

		/// returns true if the discrete function nr_fct is defined everywhere
		[[nodiscard]] bool is_def_everywhere(size_t fct) const {return m_spDDI->is_def_everywhere(fct);}


		///	returns the local finite element id of a function
		/// \{
		[[nodiscard]] const LFEID& local_finite_element_id(size_t fct) const {return m_spDDI->local_finite_element_id(fct);}
		[[nodiscard]] const LFEID& lfeid(size_t fct) const {return m_spDDI->lfeid(fct);}
		/// \}


		/// returns the maximum number of dofs on grid objects in a dimension
		[[nodiscard]] size_t max_dofs(const int dim) const {return m_spDDI->max_dofs(dim);}

		/// returns the maximum number of dofs on a grid base object type
		[[nodiscard]] size_t max_dofs(const GridBaseObjectId gbo) const {return m_spDDI->max_dofs(gbo);}

		/// returns the maximum number of dofs on reference object type
		[[nodiscard]] size_t max_dofs(const ReferenceObjectID roid) const {return m_spDDI->max_dofs(roid);}

		/// returns the maximum number of dofs in a dimension on a subset
		[[nodiscard]] size_t max_dofs(const int dim, const int si) const {return m_spDDI->max_dofs(dim, si);}

		/// returns the maximum number of dofs on a grid base object on a subset
		[[nodiscard]] size_t max_dofs(const GridBaseObjectId gbo, const int si) const {return m_spDDI->max_dofs(gbo, si);}

		///	returns the number of dofs on a Reference Object on a subset
		[[nodiscard]] size_t num_dofs(const ReferenceObjectID roid, const int si) const {return m_spDDI->num_dofs(roid, si);}


		/// returns the maximal number of dofs on a dimension for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, int dim) const {return m_spDDI->max_fct_dofs(fct, dim);}

		/// returns the maximal number of dofs on a base object type for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const GridBaseObjectId gbo) const {return m_spDDI->max_fct_dofs(fct, gbo);}

		///	returns the number of dofs on a reference object for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const ReferenceObjectID roid) const {return m_spDDI->max_fct_dofs(fct, roid);}

		/// returns the maximum number of dofs in a dimension on a subset for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const int dim, const int si) const {return m_spDDI->max_fct_dofs(fct, dim, si);}

		/// returns the maximum number of dofs on a grid base object on a subset for a function component
		[[nodiscard]] size_t max_fct_dofs(const size_t fct, const GridBaseObjectId gbo, const int si) const {return m_spDDI->max_fct_dofs(fct, gbo, si);}

		///	returns the number of dofs on a Reference Object on a subset for a function component
		[[nodiscard]] size_t num_fct_dofs(const size_t fct, const ReferenceObjectID roid, const int si) const {return m_spDDI->num_fct_dofs(fct, roid, si);}


		///	returns the offset for reference element, subset and function
		[[nodiscard]] size_t offset(const ReferenceObjectID roid, const int si, const size_t fct) const {return m_spDDI->offset(roid,si,fct);}


		///	prints statistic on local dof distribution
		void print_local_dof_statistic() const {print_local_dof_statistic(1);}

		///	prints informations
		void print_local_dof_statistic(int verboseLev) const {m_spDDI->print_local_dof_statistic(verboseLev);}

	protected:
		///	Function Pattern
		ConstSmartPtr<DoFDistributionInfo> m_spDDI;
};


} // end namespace ug

#endif