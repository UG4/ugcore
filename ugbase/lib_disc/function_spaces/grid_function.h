/*
 * grid_function.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION__

#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/dof_manager/function_pattern.h"
#include "lib_disc/common/local_algebra.h"
#include "local_transfer_interface.h"
#include "approximation_space.h"

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
	  	public IGridFunction
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

	///	type of multi indices
		typedef DoFDistribution::multi_index_type multi_index_type;

	public:
	///	iterator traits
		template <typename TElem>
		struct traits
		{
			typedef typename DoFDistribution::traits<TElem>::geometric_object geometric_object;
			typedef typename DoFDistribution::traits<TElem>::iterator iterator;
			typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
		};

		template <int dim>
		struct dim_traits
		{
			typedef typename DoFDistribution::dim_traits<dim>::geometric_base_object geometric_base_object;
			typedef typename DoFDistribution::dim_traits<dim>::iterator iterator;
			typedef typename DoFDistribution::dim_traits<dim>::const_iterator const_iterator;
		};

		typedef typename dim_traits<dim>::geometric_base_object element_type;
		typedef typename dim_traits<dim>::iterator element_iterator;
		typedef typename dim_traits<dim>::const_iterator const_element_iterator;

		typedef typename dim_traits<dim-1>::geometric_base_object side_type;
		typedef typename dim_traits<dim-1>::iterator side_iterator;
		typedef typename dim_traits<dim-1>::const_iterator const_side_iterator;

	protected:
	/// virtual clone using covariant return type
		virtual this_type* virtual_clone() const {return new this_type(*this);}

	/// virtual clone using covariant return type excluding values
		virtual this_type* virtual_clone_without_values() const;

	public:
	/// Initializing Constructor
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > approxSpace,
		             SmartPtr<DoFDistribution> spDoFDistr, bool bManage = true);

	/// Initializing Constructor using surface dof distribution
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > approxSpace, bool bManage = true);

	/// Initializing Constructor using surface dof distribution on a level
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > approxSpace, int level, bool bManage = true);

	/// Initializing Constructor using a grid level
		GridFunction(SmartPtr<ApproximationSpace<TDomain> > approxSpace, const GridLevel& gl, bool bManage = true);

	///	checks the algebra
		void check_algebra();

	/// Copy constructor
		GridFunction(const this_type& v) : IGridFunction(v) {assign(v);}

	///	assigns another grid function
		this_type& operator=(const this_type& v) {assign(v); return *this;}

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

	///	\copydoc IGridFunction::resize_values
		virtual void resize_values(size_t s, number defaultValue = 0.0);

	///	\copydoc IGridFunction::permute_values
		virtual	void permute_values(const std::vector<size_t>& vIndNew);

	///	\copydoc IGridFunction::copy_values
		virtual void copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                         bool bDisjunct = false);

	/// Destructor
		virtual ~GridFunction() {m_spDD->unmanage_grid_function(*this);}

	public:
	///	returns dof distribution
		SmartPtr<DoFDistribution> dof_distribution() {return m_spDD;}

	///	returns dof distribution
		ConstSmartPtr<DoFDistribution> dof_distribution() const {return m_spDD;}

	///	returns the grid level
		const GridLevel& grid_level() {return m_spDD->grid_level();}

	/// number of discrete functions
		size_t num_fct() const {return m_spDD->num_fct();}

	/// number of discrete functions on subset si
		size_t num_fct(int si) const {return m_spDD->num_fct(si);}

	/// returns the trial space of the discrete function fct
		LFEID local_finite_element_id(size_t fct) const
			{return m_spDD->local_finite_element_id(fct);}

	/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {return m_spDD->name(fct);}

	///	returns subset id by name
		size_t subset_id_by_name(const char* name) const {return this->subset_grp_by_name(name)[0];}

	///	returns subset group by name
		SubsetGroup subset_grp_by_name(const char* names) const {return m_spDD->subset_grp_by_name(names);}

	/// returns fct id by name
		size_t fct_id_by_name(const char* name) const{return m_spDD->fct_id_by_name(name);}

	///	returns a function group to a string of functions
		FunctionGroup fct_grp_by_name(const char* names) const {return m_spDD->fct_grp_by_name(names);}

	/// returns the dimension in which solution lives
		int get_dim(size_t fct) const {return m_spDD->dim(fct);}

	///	returns function pattern
		const FunctionPattern& function_pattern() const {return m_spDD->function_pattern();}

	/// returns true if the discrete function nr_fct is defined on subset s
		bool is_def_in_subset(size_t fct, int si) const {return m_spDD->is_def_in_subset(fct, si);}

	/// returns true if the discrete function nr_fct is defined everywhere
		bool is_def_everywhere(size_t fct) const {return m_spDD->is_def_everywhere(fct);}

	/// number of subsets
		int num_subsets() const {return m_spDD->num_subsets();}

	/// iterator for elements where this grid function is defined
	/// \{
		template <typename TElem>
		typename traits<TElem>::const_iterator begin() const
			{return m_spDD->template begin<TElem>();}

		template <typename TElem>
		typename traits<TElem>::const_iterator end() const
			{return m_spDD->template end<TElem>();}

		template <typename TElem>
		typename traits<TElem>::const_iterator begin(int si) const
			{return m_spDD->template begin<TElem>(si);}

		template <typename TElem>
		typename traits<TElem>::const_iterator end(int si) const
			{return m_spDD->template end<TElem>(si);}
	/// \}

	///	returns the adjacend elements
		template <typename TElem, typename TBaseElem>
		void collect_associated(std::vector<TBaseElem*>& vAssElem,
								TElem* elem, bool clearContainer = true) const{
				m_spDD->collect_associated(vAssElem, elem, clearContainer);
		}

	public:
	/// return the number of dofs distributed
		size_t num_indices() const {return m_spDD->num_indices();}

	/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return m_spDD->num_indices(si);}

	/// returns the maximum number of dofs on grid objects in a dimension on a subset
		size_t max_dofs(const int dim, const int si) const {return m_spDD->max_dofs(dim,si);}

	/// return the maximum number of dofs on grid objects in a dimension
		size_t max_dofs(const int dim) const {return m_spDD->max_dofs(dim);}

	/// get all indices of the element
		template <typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const
			{m_spDD->indices(elem, ind, bHang);}

	/// get multi indices on an finite element in canonical order
		template <typename TElem>
		size_t multi_indices(TElem* elem, size_t fct, std::vector<multi_index_type>& ind, bool bHang = false, bool bClear = true) const
			{return m_spDD->multi_indices(elem, fct, ind, bHang, bClear);}

	/// get multi indices on an geometric object in canonical order
		template <typename TElem>
		size_t inner_multi_indices(TElem* elem, size_t fct,	std::vector<multi_index_type>& ind, bool bClear = true) const
			{return m_spDD->inner_multi_indices(elem, fct, ind, bClear);}

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
	///	returns the position of the dofs of a function if available
		template <typename TElem>
		bool dof_positions(TElem* elem, size_t fct,
		                   std::vector<MathVector<dim> >& vPos) const;

		template <typename TElem>
		bool inner_dof_positions(TElem* elem, size_t fct,
		                         std::vector<MathVector<dim> >& vPos) const;

	public:
	///	add a transfer callback
		void add_transfer(SmartPtr<ILocalTransferAlgebra<TAlgebra> > transfer);

	///	add a transfer callback
		void add_transfer(SmartPtr<ILocalTransfer> transfer);

	///	add a transfer callback
		void remove_transfer(SmartPtr<ILocalTransfer> transfer);

	///	add a transfer callback
		void clear_transfers();

	protected:
	///	registered transfers
		std::vector<SmartPtr<ILocalTransfer> > m_vTransfer;

	protected:
	///	DoF Distribution this GridFunction relies on
		SmartPtr<DoFDistribution> m_spDD;

	/// Approximation Space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;
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
