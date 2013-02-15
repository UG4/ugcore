/*
 * dof_info.h
 *
 *  Created on: 15.02.2013
 *      Author: andreasvogel
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
	protected:
		///	indication that function is not defined on a subset
		enum{NOT_DEF_ON_SUBSET = (size_t) - 1};
		///	indication that function is not defined on a subset
		enum{NOT_YET_ASSIGNED = (size_t) - 2};

	public:
		/// constructor
		DoFDistributionInfo(ConstSmartPtr<ISubsetHandler> spSH);

		///	initializes the DoFs
		void init();

		///	returns the function pattern
		const FunctionPattern& function_pattern() const {return *this;}


		/// return the maximum number of dofs on grid objects in a dimension
		size_t max_dofs(const int dim) const {return m_vMaxDoFsInDim[dim];}

		///	returns if indices are defined on a geometric object
		size_t max_dofs(GeometricBaseObject gbo) const {return m_vMaxDoFsInDim[gbo];}

		///	returns if indices are defined on a reference object
		size_t max_dofs(ReferenceObjectID roid) const {return m_vMaxDoFsOnROID[roid];}

		/// returns the maximum number of dofs on grid objects in a dimension on a subset
		size_t max_dofs(const int dim, const int si) const {return m_vvMaxDoFsInDimPerSubset[dim][si];}


		///	returns the maximum number of dofs on a Reference Object
		inline size_t num_dofs(const ReferenceObjectID roid, const int si) const {
			UG_ASSERT((int)roid >= 0 && (int)roid < NUM_REFERENCE_OBJECTS, "Invalid subset.")
			UG_ASSERT(si >= 0 && si < (int)m_vvNumDoFsOnROID[roid].size(), "Invalid subset.")
			return m_vvNumDoFsOnROID[roid][si];
		}

		///	returns the number of dofs on a subelement of an element
		size_t num_dofs(size_t fct, const ReferenceObjectID roid, const ReferenceObjectID subRoid) const {return m_vNumDoFOnSubelem[fct](roid, subRoid);}


		/// returns maximal dimension where to dofs must be ordered
		int max_dim_to_order_dofs(size_t fct) const {return m_vMaxDimToOrderDoFs[fct];}


		///	prints informations
		void print_local_dof_statistic(int verboseLev) const;

		///	prints statistic on local dof distribution
		void print_local_dof_statistic() const {print_local_dof_statistic(1);}

		///	returns the offset for reference element, subset and function
		size_t offset(const ReferenceObjectID roid, const int si, const size_t fct) const {return m_vvvOffsets[roid][si][fct];}

	protected:
		/// creates offset arrays
		/// \{
		void create_offsets(ReferenceObjectID roid);
		void create_offsets();
		/// \}

	protected:
		/// offset map
		std::vector<std::vector<size_t> > m_vvvOffsets[NUM_REFERENCE_OBJECTS];

		///	number of DoFs on a reference element type on a subset
		std::vector<size_t> m_vvNumDoFsOnROID[NUM_REFERENCE_OBJECTS];

		///	maximum number of DoFs on a reference type
		size_t m_vMaxDoFsOnROID[NUM_REFERENCE_OBJECTS];

		///	maximum number of DoFs on geometric objects in a dimension
		size_t m_vMaxDoFsInDim[NUM_GEOMETRIC_BASE_OBJECTS];

		///	maximum number of DoFs on geometric objects in a dimension per subset
		std::vector<size_t> m_vvMaxDoFsInDimPerSubset[NUM_GEOMETRIC_BASE_OBJECTS];

		///	local dof sets
		std::vector<const ILocalDoFSet*> m_vLocalDoFSet[NUM_REFERENCE_OBJECTS];

		///	maximum dimensions where dofs must be ordered
		std::vector<int> m_vMaxDimToOrderDoFs;

		///	number Dofs for local DoF set and subelement of element
		std::vector<MathMatrix<NUM_REFERENCE_OBJECTS,NUM_REFERENCE_OBJECTS, int> > m_vNumDoFOnSubelem;
};


class DoFDistributionInfoProvider{
	public:
		/// constructor
		DoFDistributionInfoProvider(const DoFDistributionInfo& rDDI)
			: m_spDDI(&rDDI)
		{}

		///	returns underlying info
		const DoFDistributionInfo& dof_distribution_info() const {return *m_spDDI;}

		///	returns the subset handler
		ConstSmartPtr<ISubsetHandler> subset_handler() const {return m_spDDI->subset_handler();}

		///	returns the function pattern
		const FunctionPattern& function_pattern() const {return m_spDDI->function_pattern();}


		/// number of discrete functions on subset si
		size_t num_fct() const {return m_spDDI->num_fct();}

		/// number of discrete functions on subset si
		size_t num_fct(int si) const {return m_spDDI->num_fct(si);}

		/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {return m_spDDI->name(fct);}

		/// returns fct id by name
		size_t fct_id_by_name(const char* name) const{return m_spDDI->fct_id_by_name(name);}


		///	returns number of subsets
		int num_subsets() const {return m_spDDI->num_subsets();}

		/// returns the dimension in which solution lives
		int dim(size_t fct) const {return m_spDDI->dim(fct);}

		///	returns dimension of subset
		int dim_subset(int si) const {return m_spDDI->dim_subset(si);}

		///	returns subset name
		std::string subset_name(int si) const {return m_spDDI->subset_name(si);}

		///	returns subset group by name
		SubsetGroup subset_grp_by_name(const char* names) const;

		///	returns a function group to a string of functions
		FunctionGroup fct_grp_by_name(const char* names) const;


		///	retruns if a function is defined on a subset
		bool is_def_in_subset(size_t fct, int si) const {return m_spDDI->is_def_in_subset(fct,si);}

		/// returns true if the discrete function nr_fct is defined everywhere
		bool is_def_everywhere(size_t fct) const {return m_spDDI->is_def_everywhere(fct);}


		///	returns the local finite element id of a function
		const LFEID& local_finite_element_id(size_t fct) const {return m_spDDI->local_finite_element_id(fct);}


		/// return the maximum number of dofs on grid objects in a dimension
		size_t max_dofs(const int dim) const {return m_spDDI->max_dofs(dim);}

		///	returns if indices are defined on a geometric object
		size_t max_dofs(GeometricBaseObject gbo) const {return m_spDDI->max_dofs(gbo);}

		///	returns if indices are defined on a reference object
		size_t max_dofs(ReferenceObjectID roid) const {return m_spDDI->max_dofs(roid);}

		/// returns the maximum number of dofs on grid objects in a dimension on a subset
		size_t max_dofs(const int dim, const int si) const {return m_spDDI->max_dofs(dim, si);}


		///	returns the maximum number of dofs on a Reference Object
		size_t num_dofs(const ReferenceObjectID roid, const int si) const {return m_spDDI->num_dofs(roid, si);}

		///	returns the number of dofs on a subelement of an element
		size_t num_dofs(size_t fct, const ReferenceObjectID roid, const ReferenceObjectID subRoid) const {return m_spDDI->num_dofs(fct, roid, subRoid);}


		/// returns maximal dimension where to dofs must be ordered
		int max_dim_to_order_dofs(size_t fct) const {return m_spDDI->max_dim_to_order_dofs(fct);}


		///	prints informations
		void print_local_dof_statistic(int verboseLev) const {m_spDDI->print_local_dof_statistic(verboseLev);}

		///	returns the offset for reference element, subset and function
		size_t offset(const ReferenceObjectID roid, const int si, const size_t fct) const {return m_spDDI->offset(roid,si,fct);}

	protected:
		///	Function Pattern
		const DoFDistributionInfo* m_spDDI;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION_INFO__ */
