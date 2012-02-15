/*
 * dof_distribution.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__

#include <vector>
#include <algorithm>

#include "./function_pattern.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "lib_disc/function_spaces/grid_function.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_index_layout.h"
	#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{


// predeclaration
template <typename TDoFDistribution>
class IGridFunction;

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
			m_goc(goc), m_pFuncPattern(&dp), m_pSurfaceView(NULL) {}

	///	Constructor
		IDoFDistribution(GeometricObjectCollection goc, FunctionPattern& dp,
		                 const SurfaceView& surfView) :
			m_goc(goc), m_pFuncPattern(&dp), m_pSurfaceView(&surfView) {}

		///////////////////////////
		// Infos
		///////////////////////////

	///	set grouped
		void set_grouping(bool bGrouped){getImpl().set_grouping(bGrouped);}

	///	returns if a trial space is supported
		static bool supports_trial_space(LFEID& id)
		{
			return implementation_type::supports_trial_space(id);
		}

	/// dimension of subset
		int dim_subset(int si) const {return m_pFuncPattern->dim_subset(si);}

	/// name of subset
		const char* subset_name(int si) const {return m_pFuncPattern->subset_name(si);}

	/// number of discrete functions
		size_t num_fct() const {return m_pFuncPattern->num_fct();}

	/// number of discrete functions on a subset
		size_t num_fct(int si) const {return m_pFuncPattern->num_fct(si);}

	/// returns the trial space of the discrete function fct
		LFEID local_finite_element_id(size_t fct) const
			{return m_pFuncPattern->local_finite_element_id(fct);}

	/// name of discrete function
		std::string name(size_t fct) const
			{return m_pFuncPattern->name(fct);}

	/// dimension where discrete function lives
		int dim(size_t fct) const {return m_pFuncPattern->dim(fct);}

	/// returns true if the discrete function is defined on the subset
		bool is_def_in_subset(size_t fct, int si) const
			{return m_pFuncPattern->is_def_in_subset(fct, si);}

	/// returns true if the discrete function is defined everywhere
		bool is_def_everywhere(size_t fct) const
			{return m_pFuncPattern->is_def_everywhere(fct);}

	/// returns function id for local function on subset
		size_t fct_id(size_t loc_fct, int si) const
			{return m_pFuncPattern->fct_id(loc_fct, si);}

	///	returns the function pattern
		const FunctionPattern& get_function_pattern() const
			{return *m_pFuncPattern;}

		///////////////////////////////////////
		// Elements where dofs are distributed
		///////////////////////////////////////

	/// number of subsets where dofs are defined
		int num_subsets() const {return m_goc.num_levels();}

	/// number of elements in a subset where dofs are defined
		template<typename TElem>
		size_t num(int si) const {return m_goc.num<TElem>(si);}

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

	///	returns if the dof distribution distributes dofs on a given element type
		bool has_indices_on(ReferenceObjectID roid) const {return getImpl().has_indices_on(roid);}

	///	returns if the dof distribution distributes dofs on a base element type
		bool has_indices_on(GeometricBaseObject gbo) const {return getImpl().has_indices_on(gbo);}

	/// return the number of distributed indices
		size_t num_indices() const {return getImpl().num_indices();}

	/// return the number of distributed indices on subset si
		size_t num_indices(int si) const {return getImpl().num_indices(si);}

	///	Size algebra block on all subset
		int blocksize() const {return getImpl().blocksize();}

	///	Size algebra block on subset
		int blocksize(int si) const {return getImpl().blocksize(si);}

	///	prints information about the local dof pattern
		void print_local_dof_statistic(int verboseLev = 1) const
			{return getImpl().print_local_dof_statistic(verboseLev);}

		///////////////////////////////////////
		// Index Access
		///////////////////////////////////////

	///	returns all indices of the element
		template<typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const
			{getImpl().indices(elem, ind, bHang);}

	//	overloads for base types, forwarding call to concrete type
		void indices(GeometricObject* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(VertexBase* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(EdgeBase* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(Face* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(Volume* elem, LocalIndices& ind, bool bHang = false) const;

	/// get multi indices (Element + Closure of Element)
		template<typename TElem>
		size_t multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
			{return getImpl().multi_indices(elem, fct, ind);}

	//	overloads for base types, forwarding call to concrete type
		size_t multi_indices(GeometricObject* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t multi_indices(VertexBase* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t multi_indices(EdgeBase* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t multi_indices(Face* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t multi_indices(Volume* elem, size_t fct, multi_index_vector_type& ind) const;

	/// get multi indices (only inner part of Element)
		template<typename TElem>
		size_t inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
			{return getImpl().inner_multi_indices(elem, fct, ind);}

	//	overloads for base types, forwarding call to concrete type
		size_t inner_multi_indices(GeometricObject* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t inner_multi_indices(VertexBase* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t inner_multi_indices(EdgeBase* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t inner_multi_indices(Face* elem, size_t fct, multi_index_vector_type& ind) const;
		size_t inner_multi_indices(Volume* elem, size_t fct, multi_index_vector_type& ind) const;

	/// get algebra indices (Element + Closure of Element)
		template<typename TElem>
		size_t algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
			{return getImpl().algebra_indices(elem, ind);}

	//	overloads for base types, forwarding call to concrete type
		size_t algebra_indices(GeometricObject* elem, algebra_index_vector_type& ind) const;
		size_t algebra_indices(VertexBase* elem, algebra_index_vector_type& ind) const;
		size_t algebra_indices(EdgeBase* elem, algebra_index_vector_type& ind) const;
		size_t algebra_indices(Face* elem, algebra_index_vector_type& ind) const;
		size_t algebra_indices(Volume* elem, algebra_index_vector_type& ind) const;

	/// get algebra indices (only inner part of Element)
		template<typename TElem>
		size_t inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
			{return getImpl().inner_algebra_indices(elem,ind);}

	//	overloads for base types, forwarding call to concrete type
		size_t inner_algebra_indices(GeometricObject* elem, algebra_index_vector_type& ind) const;
		size_t inner_algebra_indices(VertexBase* elem, algebra_index_vector_type& ind) const;
		size_t inner_algebra_indices(EdgeBase* elem, algebra_index_vector_type& ind) const;
		size_t inner_algebra_indices(Face* elem, algebra_index_vector_type& ind) const;
		size_t inner_algebra_indices(Volume* elem, algebra_index_vector_type& ind) const;

		///////////////////////////
		// Creation
		///////////////////////////

	///	schedule grid function for adaption
		void manage_grid_function(IGridFunction<TImpl>& gf);

	///	unschedule grid function for adaption
		void unmanage_grid_function(IGridFunction<TImpl>& gf);

	///	add indices to elements
		void grid_obj_added(GeometricObject* obj);
		void grid_obj_added(VertexBase* vrt){getImpl().grid_obj_added(vrt);}
		void grid_obj_added(EdgeBase* edge) {getImpl().grid_obj_added(edge);}
		void grid_obj_added(Face* face) 	{getImpl().grid_obj_added(face);}
		void grid_obj_added(Volume* vol) 	{getImpl().grid_obj_added(vol);}

	///	remove indices from elements
		void grid_obj_to_be_removed(GeometricObject* obj);
		void grid_obj_to_be_removed(VertexBase* vrt){getImpl().grid_obj_to_be_removed(vrt);}
		void grid_obj_to_be_removed(EdgeBase* edge) {getImpl().grid_obj_to_be_removed(edge);}
		void grid_obj_to_be_removed(Face* face) 	{getImpl().grid_obj_to_be_removed(face);}
		void grid_obj_to_be_removed(Volume* vol) 	{getImpl().grid_obj_to_be_removed(vol);}

	//	copies the indices on replacement of objects
		void grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld){getImpl().grid_obj_replaced(vrtNew, vrtOld);}
		void grid_obj_replaced(EdgeBase* edgeNew, EdgeBase* edgeOld) {getImpl().grid_obj_replaced(edgeNew, edgeOld);}
		void grid_obj_replaced(Face* faceNew, Face* faceOld) 	{getImpl().grid_obj_replaced(faceNew, faceOld);}
		void grid_obj_replaced(Volume* volNew, Volume* volOld) 	{getImpl().grid_obj_replaced(volNew, volOld);}

	/// distribute dofs
		bool distribute_indices(){return getImpl().distribute_indices();}

	///	permutes all indices
	/**
	 * This method swaps the indices according to the passed mapping vector, i.e.
	 * it performs a permutation of the whole index set. The vector vIndNew must
	 * have the size of the number of indices and for each index it must return
	 * the new index, i.e. newIndex = vIndNew[oldIndex].
	 *
	 * \param[in]	vIndNew		mapping for each index
	 * \returns 	success flag
	 */
		bool permute_indices(std::vector<size_t>& vIndNew);

	///	compress indices
		bool defragment() {return getImpl().defragment();}

	///	returns the connectivity for the indices
	/**
	 * This method returns a vector, that contains for each index a list of
	 * all indices, that are connected to the index. This gives a representation
	 * of the connectivity graph of the dof distribution.
	 *
	 * \param[out]	vvConnection	vector with connected indices for each index
	 * \returns 	success flag
	 */
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection)
			{return  getImpl().get_connections(vvConnection);}

	protected:
	///	copy values in managed grid functions for pairs of indices
	/**
	 * This method copies the values for the pairs of indices in all
	 * managed grid function. The method is typically used to fill "holes" in
	 * vectors after adaptive grid refinement.
	 *
	 * \param[in]	vIndexMap		vector of index mappings (indexOld, indexNew)
	 * \param[in]	bDisjunct		flag, if permutation disjunct
	 * \returns 	success flag
	 */
		bool indices_swaped(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
							bool bDisjunct = false);

	///	resizes managed grid functions when number of indices changed
	/**
	 * This method resizes all managed grid function. The method is typically
	 * used to cut vectors after adaptive grid refinement.
	 *
	 * \param[in]	newSize			new vector size
	 */
		void num_indices_changed(size_t newSize);

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

	protected:
	///	vector of all handled gridfunctions
		std::vector<IGridFunction<TImpl>*> m_vManagedGridFunc;

	/// geometric object collection for this Distributor
		GeometricObjectCollection m_goc;

	/// Function Pattern
		FunctionPattern* m_pFuncPattern;

	///	Surface View for handling of shadows
		const SurfaceView* m_pSurfaceView;

#ifdef UG_PARALLEL
	public:
	/// returns the slave index layout for domain decompostion level
		IndexLayout& get_slave_layout()	{return m_slaveLayout;}

	/// returns the master index layout for domain decompostion level
		IndexLayout& get_master_layout(){return m_masterLayout;}

		IndexLayout& get_vertical_slave_layout()  {return m_verticalSlaveLayout;}
		IndexLayout& get_vertical_master_layout() {return m_verticalMasterLayout;}

		pcl::ParallelCommunicator<IndexLayout>& get_communicator()	{return m_communicator;}
		pcl::ProcessCommunicator& get_process_communicator()	{return m_processCommunicator;}
		const pcl::ProcessCommunicator& get_process_communicator() const	{return m_processCommunicator;}

		size_t num_master_indices() const {return num_indices(*const_cast<IndexLayout*>(&m_masterLayout));}
		size_t num_slave_indices() const {return num_indices(*const_cast<IndexLayout*>(&m_slaveLayout));}

		size_t num_vertical_master_indices() const {return num_indices(*const_cast<IndexLayout*>(&m_verticalMasterLayout));}
		size_t num_vertical_slave_indices() const {return num_indices(*const_cast<IndexLayout*>(&m_verticalSlaveLayout));}

	protected:
		size_t num_indices(IndexLayout& Layout) const
		{
			size_t sum = 0;
			for(IndexLayout::iterator iter = Layout.begin();
					iter != Layout.end(); ++iter)
				sum += Layout.interface(iter).size();
			return sum;
		}

	protected:
	//	(horizontal) master index layout
		IndexLayout m_masterLayout;

	//	(horizontal) slave index layout
		IndexLayout m_slaveLayout;

	///	vertical master index layout
		IndexLayout m_verticalMasterLayout;

	///	vertical slave index layout
		IndexLayout m_verticalSlaveLayout;

	// 	process communicator
		pcl::ProcessCommunicator m_processCommunicator;

	// 	communicator
		pcl::ParallelCommunicator<IndexLayout> m_communicator;
#endif
};

} // end namespace ug

#include "dof_distribution_impl.h"

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_DISTRIBUTION__ */
