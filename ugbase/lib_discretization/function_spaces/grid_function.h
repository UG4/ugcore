/*
 * grid_function_space.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__

#include "lib_algebra/operator/operator_base_interface.h"
#include "lib_discretization/dof_manager/dof_distribution.h"

namespace ug{

// predeclaration
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
class ApproximationSpace;

template <typename TImpl>
class IDoFDistribution;

/// Base class for all Grid Functions
/**
 * This class is the base class for all grid functions. It basically only
 * stores the Dof distribution and registers itself at the DoFDistribution on
 * creation, such that the Grid function is adapted when the Distribution is
 * changed.
 */
template <typename TDoFDistribution>
class IGridFunction
{
	public:
	///	DoF Distribution used
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Multi index
		typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;

	/// Algebra index
		typedef typename dof_distribution_type::algebra_index_vector_type algebra_index_vector_type;

	public:
	///	default Constructor
		IGridFunction() : m_pDD(NULL) {};

	///	assigns the dof distribution
		void assign_dof_distribution(dof_distribution_type& DoFDistr, bool adapt = true);

	/// get assigned dof distribution
		const dof_distribution_type& get_dof_distribution() const
			{check(); return *m_pDD;}

	/// get assigned dof distribution
		dof_distribution_type& get_dof_distribution() {check(); return *m_pDD;}

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
		virtual bool permute_values(const std::vector<size_t>& vIndNew) = 0;

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
		virtual bool copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                         bool bDisjunct = false) = 0;

	///	resize
	/**
	 * This method resizes the length of the vector.
	 *
	 * \param[in]	s		new size
	 */
		virtual void resize_values(size_t s, number defaultValue = 0.0) = 0;

	///	Destructor
		virtual ~IGridFunction() {check(); m_pDD->unmanage_grid_function(*this);}

		/////////////////////////////////
		// DoF Distribution
		/////////////////////////////////

	/// number of discrete functions
		size_t num_fct() const {check(); return m_pDD->num_fct();}

	/// number of discrete functions on subset si
		size_t num_fct(int si) const {check(); return m_pDD->num_fct(si);}

	/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const
			{check(); return m_pDD->local_shape_function_set_id(fct);}

	/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {check(); return m_pDD->name(fct);}

	/// returns the dimension in which solution lives
		int dim(size_t fct) const {check(); return m_pDD->dim(fct);}

	/// returns true if the discrete function nr_fct is defined on subset s
		bool is_def_in_subset(size_t fct, int si) const {check(); return m_pDD->is_def_in_subset(fct, si);}

	/// returns true if the discrete function nr_fct is defined everywhere
		bool is_def_everywhere(size_t fct) const {check(); return m_pDD->is_def_everywhere(fct);}

	/// number of subsets
		int num_subsets() const {check(); return m_pDD->num_subsets();}

	/// return the number of dofs distributed
		size_t num_dofs() const {check(); return m_pDD->num_dofs();}

	/// return the number of dofs distributed on subset si
		size_t num_dofs(int si) const {check(); return m_pDD->num_dofs(si);}

	/// number of elements for a subset
		template <typename TElem>
		size_t num() const
			{check(); return m_pDD->template num<TElem>();}

	/// iterator for elements where this grid function is defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator begin() const
			{check(); return const_cast<const dof_distribution_type*>(m_pDD)->template begin<TElem>();}

	/// iterator for elements where this grid function is defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator end() const
			{check(); return const_cast<const dof_distribution_type*>(m_pDD)->template end<TElem>();}

	/// number of elements of this type for a subset
		template <typename TElem>
		size_t num(int si) const
			{check(); return m_pDD->template num<TElem>(si);}

	/// iterator for elements where this grid function is defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator begin(int si) const
			{check(); return const_cast<const dof_distribution_type*>(m_pDD)->template begin<TElem>(si);}

	/// iterator for elements where this grid function is defined
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator end(int si) const
			{check(); return const_cast<const dof_distribution_type*>(m_pDD)->template end<TElem>(si);}

		////////// Local Algebra ////////////

	/// number of algebra indices on an element
		template<typename TElem>
		size_t num_indices(int si, const FunctionGroup& funcGroup) const
			{check(); return m_pDD->num_indices<TElem>(si, funcGroup);}

	/// number of algebra indices on an element
		template<typename TElem>
		size_t num_inner_indices(int si, const FunctionGroup& funcGroup) const
			{check(); return m_pDD->num_inner_indices<TElem>(si, funcGroup);}

	/// fill local informations in LocalIndex
		template<typename TElem>
		bool prepare_indices(int si, LocalIndices& ind, bool useHanging = false) const
			{check(); return m_pDD->prepare_indices<TElem>(si, ind, useHanging);}

	/// fill local informations in LocalIndex
		template<typename TElem>
		bool prepare_inner_indices(int si, LocalIndices& ind) const
			{check(); return m_pDD->prepare_inner_indices<TElem>(si, ind);}

	/// fill the global algebra indices in LocalIndex
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind, bool useHanging = false) const
			{check(); m_pDD->update_indices(elem, ind, useHanging);}

	/// fill the global algebra indices in LocalIndex
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const
			{check(); m_pDD->update_inner_indices(elem, ind);}

		////////// Multi indices ////////////

	/// number of multi indices on an finite element in canonical order
		template <typename TElem>
		size_t num_multi_indices(TElem* elem, size_t fct) const
			{check(); return m_pDD->num_multi_indices(elem, fct);}

	/// number of multi indices on an geometric object in canonical order
		template <typename TGeomObj>
		size_t num_inner_multi_indices(TGeomObj* elem, size_t fct) const
			{check(); return m_pDD->num_inner_multi_indices(elem, fct);}

	/// get multi indices on an finite element in canonical order
		template <typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
			{check(); return m_pDD->get_multi_indices(elem, fct, ind);}

	/// get multi indices on an geometric object in canonical order
		template <typename TGeomObj>
		size_t get_inner_multi_indices(TGeomObj* elem, size_t fct,	multi_index_vector_type& ind) const
			{check(); return m_pDD->get_inner_multi_indices(elem, fct, ind);}

		////////// Algebra indices ////////////

	/// number of algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		size_t num_algebra_indices(TGeomObj* elem, size_t fct) const
			{check(); return m_pDD->num_algebra_indices(elem, fct);}

	/// number of algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		size_t num_inner_algebra_indices(TGeomObj* elem, size_t fct) const
			{check(); return m_pDD->num_inner_algebra_indices(elem, fct);}

	/// get algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		void get_algebra_indices(TGeomObj* elem, algebra_index_vector_type& ind) const
			{check(); m_pDD->get_algebra_indices(elem, ind);}

	/// get algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		void get_inner_algebra_indices(TGeomObj* elem, algebra_index_vector_type& ind) const
			{check(); m_pDD->get_inner_algebra_indices(elem, ind);}

	protected:
	//	check that object can be used
		void check() const {UG_ASSERT(m_pDD != NULL, "DoF Distribution not set.\n");}

	protected:
	//	DoF Distribution this GridFunction relies on
		dof_distribution_type* m_pDD;
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
 * \tparam	TDoFDistribution	dof distribution type
 * \tparam	TAlgebra			algebra type
 */
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
class GridFunction
	: 	public TAlgebra::vector_type,
	  	public IGridFunction<TDoFDistribution>
{
	public:
	///	This type
		typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> this_type;

	///	Type of Approximation space
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

	///	Domain
		typedef TDomain domain_type;

	///	World Dimension
		static const int dim = domain_type::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type used to store dof values
		typedef typename algebra_type::vector_type vector_type;

	///	Local vector type
		typedef LocalVector<typename vector_type::value_type> local_vector_type;

	///	Local index type
		typedef LocalIndices local_index_type;

	/// DoF Distribution type
		typedef typename IGridFunction<TDoFDistribution>::dof_distribution_type
				dof_distribution_type;

	public:
	/// Default constructor
		GridFunction() : m_pApproxSpace(NULL){}

	/// Initializing Constructor
		GridFunction(	approximation_space_type& approxSpace,
						dof_distribution_type& DoFDistr)
			: m_pApproxSpace(&approxSpace)
		{
			IGridFunction<TDoFDistribution>::assign_dof_distribution(DoFDistr);
		};

	/// Copy constructor
		GridFunction(const this_type& v) {assign(v);}

	///	assigns another grid function
		this_type& operator=(const this_type& v) {assign(v); return *this;}

	/// clone
		this_type& clone(){return *(new this_type(*this));}

	/// copies the GridFunction v, except that the values are copied.
		virtual void clone_pattern(const this_type& v);

	///	assigns another GridFunction
		bool assign(const this_type& v);

	///	assigns the values of a vector
		bool assign(const vector_type& v);

	///	\copydoc IGridFunction::resize_values
		virtual void resize_values(size_t s, number defaultValue = 0.0);

	///	\copydoc IGridFunction::permute_values
		virtual	bool permute_values(const std::vector<size_t>& vIndNew);

	///	\copydoc IGridFunction::copy_values
		virtual bool copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                         bool bDisjunct = false);

	/// Destructor
		virtual ~GridFunction() {}

	/// returns world dimension
		int get_dim() const {return domain_type::dim;}

	///	assign approximation space
		void assign_approximation_space(approximation_space_type& ApproxSpace)
			{m_pApproxSpace = &ApproxSpace;}

	/// returns approximation space
		const approximation_space_type& get_approximation_space() const
			{return *m_pApproxSpace;}

	///	returns approximation space
		approximation_space_type& get_approximation_space()
			{return *m_pApproxSpace;}

	///	returns domain
		domain_type& get_domain() {return m_pApproxSpace->get_domain();}

	///	returns const domain
		const domain_type& get_domain() const {return m_pApproxSpace->get_domain();}

	/// access dof values
		inline number get_dof_value(size_t i, size_t comp) const
			{return BlockRef( (vector_type::operator[](i)), comp);}

	protected:
	/// Approximation Space
		approximation_space_type* m_pApproxSpace;
};

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
const typename TAlgebra::vector_type &getVector(const GridFunction<TDomain, TDoFDistribution, TAlgebra> &t)
{
	return *dynamic_cast<const GridFunction<TDomain, TDoFDistribution, TAlgebra>*>(&t);
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
inline std::ostream& operator<< (std::ostream& outStream, const GridFunction<TDomain, TDoFDistribution, TAlgebra>& v)
{
	outStream << *dynamic_cast<const GridFunction<TDomain, TDoFDistribution, TAlgebra>*>(&v);
	return outStream;
}

} // end namespace ug

// include implementation
#include "grid_function_impl.h"

#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__ */
