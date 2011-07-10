/*
 * dof_distribution.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__

#include <vector>
#include <algorithm>

#include "./function_pattern.h"
#include "lib_discretization/common/function_group.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "lib_discretization/function_spaces/grid_function.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_index_layout.h"
	#include "lib_discretization/parallelization/parallelization_util.h"
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

	///	returns if a trial space is supported
		static bool supports_trial_space(LSFSID& id)
		{
			return implementation_type::supports_trial_space(id);
		}

	/// dimension of subset
		int dim_subset(int si) const {return m_pFuncPattern->dim_subset(si);}

	/// number of discrete functions
		size_t num_fct() const {return m_pFuncPattern->num_fct();}

	/// number of discrete functions on a subset
		size_t num_fct(int si) const {return m_pFuncPattern->num_fct(si);}

	/// returns the trial space of the discrete function fct
		LSFSID local_shape_function_set_id(size_t fct) const
			{return m_pFuncPattern->local_shape_function_set_id(fct);}

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
		template <typename TElem>
		bool has_dofs_on() const {return getImpl().has_dofs_on<TElem>();}

	/// return the number of dofs distributed
		size_t num_dofs() const {return getImpl().num_dofs();}

	/// return the number of dofs distributed on subset si
		size_t num_dofs(int si) const {return getImpl().num_dofs(si);}

	///	Size algebra block on all subset
		int blocksize() const {return getImpl().blocksize();}

	///	Size algebra block on subset
		int blocksize(int si) const {return getImpl().blocksize(si);}

		///////////////////////////////////////
		// Index Access
		///////////////////////////////////////

	///	returns all indices of the element
		template<typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const
			{getImpl().indices(elem, ind, bHang);}

	/// get multi indices (Element + Closure of Element)
		template<typename TElem>
		size_t multi_indices(TElem* elem, size_t fct,
		                         multi_index_vector_type& ind) const
			{return getImpl().multi_indices(elem, fct, ind);}

	/// get multi indices (only inner part of Element)
		template<typename TElem>
		size_t inner_multi_indices(TElem* elem, size_t fct,
		                               multi_index_vector_type& ind) const
			{return getImpl().inner_multi_indices(elem, fct, ind);}

	/// get algebra indices (Element + Closure of Element)
		template<typename TElem>
		size_t algebra_indices(TElem* elem,
		                         algebra_index_vector_type& ind) const
			{return getImpl().algebra_indices(elem, ind);}

	/// get algebra indices (only inner part of Element)
		template<typename TElem>
		size_t inner_algebra_indices(TElem* elem,
		                               algebra_index_vector_type& ind) const
			{return getImpl().inner_algebra_indices(elem,ind);}

		///////////////////////////
		// Creation
		///////////////////////////

	///	schedule grid function for adaption
		void manage_grid_function(IGridFunction<TImpl>& gf)
		{
		//	search for grid function in list
			typename std::vector<IGridFunction<TImpl>*>::iterator it;
			it = find(m_vManagedGridFunc.begin(), m_vManagedGridFunc.end(), &gf);

		//	add if not found
			if(it == m_vManagedGridFunc.end())
				m_vManagedGridFunc.push_back(&gf);
		}

	///	unschedule grid function for adaption
		void unmanage_grid_function(IGridFunction<TImpl>& gf)
		{
		//	search for grid function in list
			typename std::vector<IGridFunction<TImpl>*>::iterator it;
			it = find(m_vManagedGridFunc.begin(), m_vManagedGridFunc.end(), &gf);

		//	remove if found
			if(it != m_vManagedGridFunc.end())
				m_vManagedGridFunc.erase(it);
		}

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
		bool permute_indices(std::vector<size_t>& vIndNew)
		{
		//	check size of permuting array. Must have same size as index set
			if(vIndNew.size() != num_dofs())
			{
				UG_LOG("ERROR in 'IDoFDistribution::permute_indices': Passed "
						" permutation does not have the size of the index set "
						<<num_dofs()<<", but has size "<<vIndNew.size()<<"\n");
				return false;
			}

		//	swap indices as implemented
			if(!getImpl().permute_indices(vIndNew)) return false;

		//	in parallel adjust also the layouts
#ifdef UG_PARALLEL
			PermuteIndicesInIndexLayout(m_slaveLayout, vIndNew);
			PermuteIndicesInIndexLayout(m_masterLayout, vIndNew);
			PermuteIndicesInIndexLayout(m_verticalSlaveLayout, vIndNew);
			PermuteIndicesInIndexLayout(m_verticalMasterLayout, vIndNew);
#endif

		//	swap values of handled grid functions
			for(size_t i = 0; i < m_vManagedGridFunc.size(); ++i)
				m_vManagedGridFunc[i]->permute_values(vIndNew);

		//	we're done
			return true;
		}

	///	compress indices
		bool defragment() {return getImpl().defragment();}

	///	add indices to elements
		void grid_obj_added(GeometricObject* obj)
		{
			uint type = obj->base_object_type_id();
			switch(type)
			{
				case VERTEX:grid_obj_added(reinterpret_cast<VertexBase*>(obj)); return;
				case EDGE: 	grid_obj_added(reinterpret_cast<EdgeBase*>(obj)); return;
				case FACE:	grid_obj_added(reinterpret_cast<Face*>(obj)); return;
				case VOLUME:grid_obj_added(reinterpret_cast<Volume*>(obj)); return;
			}
			throw(UGFatalError("GeomObject type not known."));
		}
		void grid_obj_added(VertexBase* vrt){getImpl().grid_obj_added(vrt);}
		void grid_obj_added(EdgeBase* edge) {getImpl().grid_obj_added(edge);}
		void grid_obj_added(Face* face) 	{getImpl().grid_obj_added(face);}
		void grid_obj_added(Volume* vol) 	{getImpl().grid_obj_added(vol);}

	///	remove indices from elements
		void grid_obj_to_be_removed(GeometricObject* obj)
		{
			uint type = obj->base_object_type_id();
			switch(type)
			{
				case VERTEX:grid_obj_to_be_removed(reinterpret_cast<VertexBase*>(obj)); return;
				case EDGE: 	grid_obj_to_be_removed(reinterpret_cast<EdgeBase*>(obj)); return;
				case FACE:	grid_obj_to_be_removed(reinterpret_cast<Face*>(obj)); return;
				case VOLUME:grid_obj_to_be_removed(reinterpret_cast<Volume*>(obj)); return;
			}
			throw(UGFatalError("GeomObject type not known."));
		}
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
		bool distribute_dofs(){return getImpl().distribute_dofs();}

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
							bool bDisjunct = false)
		{
		//	swap values of handled grid functions
			for(size_t i = 0; i < m_vManagedGridFunc.size(); ++i)
				m_vManagedGridFunc[i]->copy_values(vIndexMap, bDisjunct);

		//	we're done
			return true;
		}

	///	resizes managed grid functions when number of indices changed
	/**
	 * This method resizes all managed grid function. The method is typically
	 * used to cut vectors after adaptive grid refinement.
	 *
	 * \param[in]	newSize			new vector size
	 */
		void num_indices_changed(size_t newSize)
		{
		//	swap values of handled grid functions
			for(size_t i = 0; i < m_vManagedGridFunc.size(); ++i)
				m_vManagedGridFunc[i]->resize_values(newSize);
		}

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

		size_t num_master_dofs() const {return num_dofs(*const_cast<IndexLayout*>(&m_masterLayout));}
		size_t num_slave_dofs() const {return num_dofs(*const_cast<IndexLayout*>(&m_slaveLayout));}

		size_t num_vertical_master_dofs() const {return num_dofs(*const_cast<IndexLayout*>(&m_verticalMasterLayout));}
		size_t num_vertical_slave_dofs() const {return num_dofs(*const_cast<IndexLayout*>(&m_verticalSlaveLayout));}

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

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_DISTRIBUTION__ */
