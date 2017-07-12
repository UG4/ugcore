/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION__

#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/dof_manager/function_pattern.h"
#include "lib_disc/common/local_algebra.h"
#include "approximation_space.h"
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug{

/// Base class for all Grid Functions
/**
 * This class is the base class for all grid functions. It basically only
 * stores the Dof distribution and registers itself at the DoFDistribution on
 * creation, such that the Grid function is adapted when the Distribution is
 * changed.
 */
class IGridFunction
{
	public:
		virtual ~IGridFunction()	{}

	///	permutes all values
	/**
	 * This method permutes the values according to the passed mapping vector, i.e.
	 * it performs a permutation of the whole index set. The vector vIndNew must
	 * have the size of the number of indices and for each index it must return
	 * the new index, i.e. newIndex = vIndNew[oldIndex].
	 *
	 * \param[in]	vIndNew		mapping for each index
	 * \returns 	success flag
	 */
		virtual void permute_values(const std::vector<size_t>& vIndNew) = 0;

	///	copy values
	/**
	 * This method copies values between indices according to the passed mapping.
	 * The copy of the values is are performed as:
	 *
	 * 	for all i:	newIndex = vIndexMap[i].second
	 * 				oldIndex = vIndexMap[i].first
	 * 				value[newIndex] <- value[oldIndex]
	 *
	 * If the copy operation is known to be a disjunct composition (i.e. each index
	 * appears only in one operation), this can be specified by a flag. In
	 * this case the order in which the copying is performed is arbitrary and
	 * this will save a copy operation of the whole vector.
	 *
	 * \param[in]	vIndexMap		vector of index mappings (indexOld, indexNew)
	 * \param[in]	bDisjunct		flag, if permutation disjunct
	 * \returns 	success flag
	 */
		virtual void copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                         bool bDisjunct = false) = 0;

	///	resize
	/**
	 * This method resizes the length of the vector.
	 *
	 * \param[in]	s				new size
	 * \param[in]	defaultValue	default value for new entries
	 */
		virtual void resize_values(size_t s, number defaultValue = 0.0) = 0;
};

template <typename TDomain> class AdaptionSurfaceGridFunction;

/// represents numerical solutions on a grid using an algebraic vector
/**
 * A grid function brings approximation space and algebra together. For a given
 * DoF Distribution (e.g. a level dof distribution or a surface dof distribution)
 * the grid function stores the values of the DoFs in an algebraic vector.
 * In addition access to the grid elements is provided and the mapping between
 * grid elements and DoFs is provided.
 *
 * \tparam 	TDomain				domain type
 * \tparam	TAlgebra			algebra type
 */
template <typename TDomain, typename TAlgebra>
class GridFunction
	: 	public TAlgebra::vector_type,
	  	public IGridFunction,
	  	public DoFDistributionInfoProvider
{
	public:
	///	This type
		typedef GridFunction<TDomain, TAlgebra> this_type;

	///	Type of Approximation space
		typedef ApproximationSpace<TDomain> approximation_space_type;

	///	Domain
		typedef TDomain domain_type;

	///	World Dimension
		static const int dim = domain_type::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type used to store dof values
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	iterator traits
		template <typename TElem>
		struct traits
		{
			typedef typename DoFDistribution::traits<TElem>::grid_object grid_object;
			typedef typename DoFDistribution::traits<TElem>::iterator iterator;
			typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
		};

		template <int dim>
		struct dim_traits
		{
			typedef typename DoFDistribution::dim_traits<dim>::grid_base_object grid_base_object;
			typedef typename DoFDistribution::dim_traits<dim>::iterator iterator;
			typedef typename DoFDistribution::dim_traits<dim>::const_iterator const_iterator;
		};

		typedef typename dim_traits<dim>::grid_base_object element_type;
		typedef typename dim_traits<dim>::iterator element_iterator;
		typedef typename dim_traits<dim>::const_iterator const_element_iterator;

		typedef typename dim_traits<dim-1>::grid_base_object side_type;
		typedef typename dim_traits<dim-1>::iterator side_iterator;
		typedef typename dim_traits<dim-1>::const_iterator const_side_iterator;

	protected:
	/// virtual clone using covariant return type
		virtual this_type* virtual_clone() const {return new this_type(*this);}

	/// virtual clone using covariant return type excluding values
		virtual this_type* virtual_clone_without_values() const;

	public:
	/// Initializing Constructor
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		             SmartPtr<DoFDistribution> spDoFDistr, bool bManage = true);

	/// Initializing Constructor using surface dof distribution
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, bool bManage = true);

	/// Initializing Constructor using surface dof distribution on a level
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, int level, bool bManage = true);

	/// Initializing Constructor using a grid level
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const GridLevel& gl, bool bManage = true);

	protected:
	///	checks the algebra
		void check_algebra();

	///	inits the grid function
		void init(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		          SmartPtr<DoFDistribution> spDoFDistr, bool bManage);

	public:
	/// Copy constructor
		GridFunction(const this_type& v) : IGridFunction(v) {assign(v);}

	///	assigns another grid function
		this_type& operator=(const this_type& v) 
		{ 
		  if (this!= &v) assign(v); 
		  return *this;
		}

        ///	assigns constant value		
		this_type& operator=(number d) 
		{
		  vector_type::operator=(d); 
		  return *this;
		}

	/// clone including values
		SmartPtr<this_type> clone() const {return SmartPtr<this_type>(this->virtual_clone());}


	/// clone excluding values
		SmartPtr<this_type> clone_without_values() const {return SmartPtr<this_type>(this->virtual_clone_without_values());}

	/// copies the GridFunction v, except that the values are copied.
		virtual void clone_pattern(const this_type& v);

	///	assigns another GridFunction
		void assign(const this_type& v);

	///	assigns the values of a vector
		void assign(const vector_type& v);

	/// Destructor
		virtual ~GridFunction() {m_spDD->unmanage_grid_function(*this);}

	public:
	///	returns dof distribution
	/// \{
		SmartPtr<DoFDistribution> dof_distribution() {return m_spDD;}
		SmartPtr<DoFDistribution> dd() {return dof_distribution();}
	/// \}

	///	returns dof distribution
	/// \{
		ConstSmartPtr<DoFDistribution> dof_distribution() const {return m_spDD;}
		ConstSmartPtr<DoFDistribution> dd() const {return dof_distribution();}
	/// \}

	///	returns the grid level
		const GridLevel& grid_level() const {return m_spDD->grid_level();}

	/// iterator for elements where this grid function is defined
	/// \{
		template <typename TElem>
		typename traits<TElem>::const_iterator begin() const
			{return m_spDD->template begin<TElem>();}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		begin(SurfaceView::SurfaceConstants validStates) const
			{return m_spDD->template begin<TElem>(validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator end() const
			{return m_spDD->template end<TElem>();}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		end(SurfaceView::SurfaceConstants validStates) const
			{return m_spDD->template end<TElem>(validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator begin(int si) const
			{return m_spDD->template begin<TElem>(si);}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		begin(int si, SurfaceView::SurfaceConstants validStates) const
			{return m_spDD->template begin<TElem>(si, validStates);}

		template <typename TElem>
		typename traits<TElem>::const_iterator end(int si) const
			{return m_spDD->template end<TElem>(si);}

		template <typename TElem>
		typename traits<TElem>::const_iterator
		end(int si, SurfaceView::SurfaceConstants validStates) const
			{return m_spDD->template end<TElem>(si, validStates);}
	/// \}

	///	returns the adjacend elements
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
								TElem* elem, bool clearContainer = true) const{
				m_spDD->collect_associated(vAssElem, elem, clearContainer);
		}

	/// returns if the grid object is part of this grid function
		template <class TGeomObj>
		bool is_contained(TGeomObj* obj) const{
			return m_spDD->is_contained(obj);
		}

	public:
	/// return the number of indices distributed (proc local)
		size_t num_indices() const {return m_spDD->num_indices();}

	/// return the number of indices distributed on subset si (proc local)
		size_t num_indices(int si) const {return m_spDD->num_indices(si);}

	/// return the number of dofs distributed (global)
	/// \{
		size_t num_dofs() const {return num_dofs(DoFCount::ALL_FCT);}
		size_t num_dofs(int fct) const {return num_dofs(fct, DoFCount::ALL_SUBSET);}
		size_t num_dofs(int fct, int si) const;
	/// \}

	/// get all indices of the element
		template <typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const
			{m_spDD->indices(elem, ind, bHang);}

	/// get multi indices on an finite element in canonical order
		template <typename TElem>
		size_t dof_indices(TElem* elem, size_t fct, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const
			{return m_spDD->dof_indices(elem, fct, ind, bHang, bClear);}

	/// get multi indices on an geometric object in canonical order
		template <typename TElem>
		size_t inner_dof_indices(TElem* elem, size_t fct,	std::vector<DoFIndex>& ind, bool bClear = true) const
			{return m_spDD->inner_dof_indices(elem, fct, ind, bClear);}

	/// get algebra indices on an geometric object in canonical order
		template <typename TElem>
		size_t algebra_indices(TElem* elem, std::vector<size_t>& ind, bool bClear = true) const
			{return m_spDD->algebra_indices(elem, ind, bClear);}

	/// get algebra indices on an geometric object in canonical order
		template <typename TElem>
		size_t inner_algebra_indices(TElem* elem, std::vector<size_t>& ind, bool bClear = true) const
			{return m_spDD->inner_algebra_indices(elem, ind, bClear);}

	public:
	///	returns domain
		SmartPtr<TDomain> domain() {return m_spApproxSpace->domain();}

	///	returns const domain
		ConstSmartPtr<TDomain> domain() const {return m_spApproxSpace->domain();}

	///	returns approx space
		SmartPtr<ApproximationSpace<TDomain> > approx_space() {return m_spApproxSpace;}

	///	returns const domain
		ConstSmartPtr<ApproximationSpace<TDomain> > approx_space() const {return m_spApproxSpace;}

	public:
	/// return m_bManaged
		bool managed() const {return m_bManaged;}
		
	///	retruns true if the grid-function is redistributed together with the grid in parallel applications
	/**	\note	if redistribution is disabled, values corresponding to elements
	 *			which were received during redistribution will not be initialized
	 *			and my be completely random.*/
		bool redistribution_enabled() const	{return m_bRedistribute;}

	///	enables or disables redistribution for this grid function
	/**	\note	if redistribution is disabled, values corresponding to elements
	 *			which were received during redistribution will not be initialized
	 *			and my be completely random.*/
		void enable_redistribution(bool enable)	{m_bRedistribute = enable;}

	///	\copydoc IGridFunction::resize_values
		virtual void resize_values(size_t s, number defaultValue = 0.0);

	///	\copydoc IGridFunction::permute_values
		virtual	void permute_values(const std::vector<size_t>& vIndNew);

	///	\copydoc IGridFunction::copy_values
		virtual void copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                         bool bDisjunct = false);

	protected:
	///	message hub id
		MessageHub::SPCallbackId m_spGridAdaptionCallbackID;
		MessageHub::SPCallbackId m_spGridDistributionCallbackID;

	///	registers at message hub for grid adaption
		void register_at_adaption_msg_hub();

	/**	this callback is called by the message hub, when a grid change has been
	 * performed. It will call all necessary actions in order to keep the grid
	 * correct for computations.*/
		void grid_changed_callback(const GridMessage_Adaption& msg);

	///	called during parallel redistribution
		void grid_distribution_callback(const GridMessage_Distribution& msg);

	/// adaption grid function for temporary storage of values in grid
		SmartPtr<AdaptionSurfaceGridFunction<TDomain> > m_spAdaptGridFct;

	protected:
	///	DoF Distribution this GridFunction relies on
		SmartPtr<DoFDistribution> m_spDD;

	/// Approximation Space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	/// boolean for DoF Distribution management of grid function
		bool m_bManaged;

	///	specifies whether the gridfunction should be redistributed together with the grid.
	/**	enabled by default */
		bool m_bRedistribute;

	///	stores the storage-type from before redistribution
		uint m_preDistStorageType;
};

template <typename TDomain, typename TAlgebra>
const typename TAlgebra::vector_type &getVector(const GridFunction<TDomain, TAlgebra> &t)
{
	return *dynamic_cast<const GridFunction<TDomain, TAlgebra>*>(&t);
}


template <typename TDomain, typename TAlgebra>
inline std::ostream& operator<< (std::ostream& outStream, const GridFunction<TDomain, TAlgebra>& v)
{
	outStream << *dynamic_cast<const GridFunction<TDomain, TAlgebra>*>(&v);
	return outStream;
}

} // end namespace ug

// include implementation
#include "grid_function_impl.h"

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION__ */
