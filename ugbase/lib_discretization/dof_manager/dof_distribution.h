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
#include "lib_discretization/common/function_group.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_grid/grid/geometric_base_objects.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_index_layout.h"
#endif

namespace ug{

/// Base class for dof distributions
/**
 * A DoF Distribution handles the distribution of Degrees of Freedom (dofs) on
 * a Domain. From the point of view of the Domain, the degrees of freedom are
 * associated with grid elements. Thus, the DoFDistribution uses the attachment
 * system of the ug::Grid to attach indices to geometric elements. On the other
 * hand from the view of the algebra, dofs are just vector entries, that can
 * be accessed by indices. For the algebra exist no grid. The DoFDistribution
 * bridges this gap.
 *
 * This interface is implemented using the Barton-Nackman trick.
 *
 * \tparam 	TImpl		Special Implementation
 */
template <typename TImpl>
class IDoFDistribution
{
	public:
	/// type of multi index
		typedef MultiIndex<2> index_type;

	/// type of value container for element local indices
		typedef std::vector<index_type> multi_index_vector_type;

	/// type of algebra index vector
		typedef std::vector<size_t> algebra_index_vector_type;

	///	type of implementation
		typedef TImpl implementation_type;

	public:
	///	Constructor
		IDoFDistribution(GeometricObjectCollection goc, FunctionPattern& dp) :
			m_goc(goc), m_pFunctionPattern(&dp) {}

		///////////////////////////
		// Infos
		///////////////////////////

	/// number of discrete functions
		size_t num_fct() const {return m_pFunctionPattern->num_fct();}

	/// number of discrete functions on a subset
		size_t num_fct(int si) const {return m_pFunctionPattern->num_fct(si);}

	/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const
			{return m_pFunctionPattern->local_shape_function_set_id(fct);}

	/// name of discrete function
		std::string name(size_t fct) const
			{return m_pFunctionPattern->name(fct);}

	/// dimension where discrete function lives
		int dim(size_t fct) const {return m_pFunctionPattern->dim(fct);}

	/// returns true if the discrete function is defined on the subset
		bool is_def_in_subset(size_t fct, int si) const
			{return m_pFunctionPattern->is_def_in_subset(fct, si);}

	/// returns function id for local function on subset
		size_t fct_id(size_t loc_fct, int si) const
			{return m_pFunctionPattern->fct_id(loc_fct, si);}

		///////////////////////////////////////
		// Elements where dofs are distributed
		///////////////////////////////////////

	/// number of subsets where dofs are defined
		int num_subsets() const {return m_goc.num_levels();}

	/// number of elements where dofs are defined
		template<typename TElem>
		size_t num() const {return m_goc.num<TElem>();}

	/// number of elements in a subset where dofs are defined
		template<typename TElem>
		size_t num(int si) const {return m_goc.num<TElem>(si);}

	///	iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::iterator begin()
			{return m_goc.begin<TElem>();}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::iterator end()
			{return m_goc.end<TElem>();}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator begin() const
			{return m_goc.begin<TElem>();}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator end() const
			{return m_goc.end<TElem>();}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::iterator begin(int si)
			{return m_goc.begin<TElem>(si);}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::iterator end(int si)
			{return m_goc.end<TElem>(si);}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator begin(int si) const
			{return m_goc.begin<TElem>(si);}

	/// iterator for elements where dofs are defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator end(int si) const
			{return m_goc.end<TElem>(si);}

		///////////////////////////
		// general informations
		///////////////////////////

	/// return the number of dofs distributed
		size_t num_dofs() const {return getImpl().num_dofs();}

	/// return the number of dofs distributed on subset si
		size_t num_dofs(int si) const {return getImpl().num_dofs(si);}

	///	Size algebra block on all subset
		int blocksize() const {return getImpl().blocksize();}

	///	Size algebra block on subset
		int blocksize(int si) const {return getImpl().blocksize(si);}

		///////////////////////////////////////
		// LocalIndex update
		///////////////////////////////////////

	/// number of algebra indices (Element + Closure of Element)
		size_t num_indices(ReferenceObjectID refID, int si,
		                   const FunctionGroup& funcGroup) const
			{return getImpl().num_indices(refID, si, funcGroup);}

	/// number of algebra indices (only inner part of Element)
		size_t num_inner_indices(ReferenceObjectID refID, int si,
		                         const FunctionGroup& funcGroup) const
			{return getImpl().num_inner_indices(refID, si, funcGroup);}

	/// fill local informations in LocalIndex (Element + Closure of Element)
		bool prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind,
		                     bool withHanging = false) const
			{return getImpl().prepare_indices(refID, si, ind, withHanging);}

	/// fill local informations in LocalIndex (only inner part of Element)
		bool prepare_inner_indices(ReferenceObjectID refID, int si,
		                           LocalIndices& ind) const
			{return getImpl().prepare_inner_indices(refID, si, ind);}

	/// fill the global algebra indices in LocalIndex (Element + Closure of Element)
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind,
		                    bool withHanging = false) const
			{getImpl().update_indices(elem, ind, withHanging);}

	/// fill the global algebra indices in LocalIndex (only inner part of Element)
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const
			{getImpl().update_inner_indices(elem, ind);}

		///////////////////////////////////////
		// Multi index access
		///////////////////////////////////////

	/// number of multi indices (Element + Closure of Element)
		template<typename TElem>
		size_t num_multi_indices(TElem* elem, size_t fct) const
			{return getImpl().num_multi_indices(elem, fct);}

	/// number of multi indices (only inner part of Element)
		template<typename TElem>
		size_t num_inner_multi_indices(TElem* elem, size_t fct) const
			{return getImpl().num_inner_multi_indices(elem, fct);}

	/// get multi indices (Element + Closure of Element)
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct,
		                         multi_index_vector_type& ind) const
			{return getImpl().get_multi_indices(elem, fct, ind);}

	/// get multi indices (only inner part of Element)
		template<typename TElem>
		size_t get_inner_multi_indices(TElem* elem, size_t fct,
		                               multi_index_vector_type& ind) const
			{return getImpl().get_inner_multi_indices(elem, fct, ind);}

		///////////////////////////////////////
		// Algebra index access
		///////////////////////////////////////

	/// number of algebra indices (Element + Closure of Element)
		template<typename TElem>
		size_t num_algebra_indices(TElem* elem, size_t fct) const
			{return getImpl().num_algebra_indices(elem, fct);}

	/// number of algebras indices (only inner part of Element)
		template<typename TElem>
		size_t num_inner_algebra_indices(TElem* elem, size_t fct) const
			{return getImpl().num_inner_algebra_indices(elem, fct);}

	/// get algebra indices (Element + Closure of Element)
		template<typename TElem>
		void get_algebra_indices(TElem* elem,
		                         algebra_index_vector_type& ind) const
			{getImpl().get_algebra_indices(elem, ind);}

	/// get algebra indices (only inner part of Element)
		template<typename TElem>
		void get_inner_algebra_indices(TElem* elem,
		                               algebra_index_vector_type& ind) const
			{getImpl().get_inner_algebra_indices(elem,ind);}

		///////////////////////////
		// Creation
		///////////////////////////

	/// distribute dofs
		bool distribute_dofs() {return getImpl().distribute_dofs();}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

	protected:
	/// geometric object collection for this Distributor
		GeometricObjectCollection m_goc;

	/// Function Pattern
		FunctionPattern* m_pFunctionPattern;

#ifdef UG_PARALLEL
	public:
	/// returns the slave index layout for domain decompostion level
		IndexLayout& get_slave_layout(size_t ddlev = 0)	{require_ddlev(ddlev); return m_vSlaveLayout.at(ddlev);}

	/// returns all slave indey layouts in a std::vector (all dd level)
		std::vector<IndexLayout>& get_slave_layouts()	{return m_vSlaveLayout;}

	/// returns the master index layout for domain decompostion level
		IndexLayout& get_master_layout(size_t ddlev = 0){require_ddlev(ddlev); return m_vMasterLayout.at(ddlev);}

	/// returns all master indey layouts in a std::vector (all dd level)
		std::vector<IndexLayout>& get_master_layouts()	{return m_vMasterLayout;}

		IndexLayout& get_vertical_slave_layout()		{return m_verticalSlaveLayout;}
		IndexLayout& get_vertical_master_layout()	{return m_verticalMasterLayout;}

		pcl::ParallelCommunicator<IndexLayout>& get_communicator(size_t ddlev = 0)	{return m_vCommunicator.at(ddlev);}
		pcl::ProcessCommunicator& get_process_communicator(size_t ddlev = 0)	{return m_vProcessCommunicator.at(ddlev);}

		std::vector<pcl::ParallelCommunicator<IndexLayout> >& get_communicators()	{return m_vCommunicator;}
		std::vector<pcl::ProcessCommunicator>& get_process_communicators()	{return m_vProcessCommunicator;}

		size_t num_master_dofs(size_t ddlev = 0) {return num_dofs(m_vMasterLayout.at(ddlev));}
		size_t num_slave_dofs(size_t ddlev = 0) {return num_dofs(m_vSlaveLayout.at(ddlev));}

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

		void require_ddlev(size_t ddlev)
		{
			if(m_vSlaveLayout.size() <= ddlev) m_vSlaveLayout.resize(ddlev+1);
			if(m_vMasterLayout.size() <= ddlev) m_vMasterLayout.resize(ddlev+1);
			if(m_vProcessCommunicator.size() <= ddlev) m_vProcessCommunicator.resize(ddlev+1);
			if(m_vCommunicator.size() <= ddlev) m_vCommunicator.resize(ddlev+1);
		}

	protected:
		// index layout for this grid level (an each domain decomposition)
		std::vector<IndexLayout> m_vSlaveLayout;

		// index layout for this grid level (an each domain decomposition)
		std::vector<IndexLayout> m_vMasterLayout;

		// index layout for each grid level
		IndexLayout m_verticalMasterLayout;

		// index layout for each grid level
		IndexLayout m_verticalSlaveLayout;

		// process communicator
		std::vector<pcl::ProcessCommunicator> m_vProcessCommunicator;

		// communicator
		std::vector<pcl::ParallelCommunicator<IndexLayout> > m_vCommunicator;
#endif
};



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__ */
