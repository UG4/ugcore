/*
 * ass_adapter.h
 *
 *  Created on: 04.02.2013
 *      Author: raphaelprohl, Andreas Vogel
 */

#ifndef ASS_ADAPTER_H_
#define ASS_ADAPTER_H_

#include "lib_grid/tools/bool_marker.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class IDomainConstraint;

/// Types of constraint
/**
 * This types control the order in with the constraints are performed.
 * constraints with a lower number will be performed first. The order of
 * constraints with equal number is undefined.
 */
enum ConstraintType
{
	CT_NONE = 0,
	CT_CONSTRAINTS = 1 << 0,
	CT_DIRICHLET = 1 << 1,
	CT_ALL = CT_NONE | CT_CONSTRAINTS | CT_DIRICHLET
};

template <typename TAlgebra>
class LocalToGlobalMapper : public ILocalToGlobalMapper<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default constructor
		LocalToGlobalMapper() {}

	///	adds a local vector to the global one
		void AddLocalVec(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
			{ AddLocalVector(vec,lvec);}

	///	adds a local matrix to the global one
		void AddLocalMatToGlobal(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
			{ AddLocalMatrixToGlobal(mat,lmat);}

	///	destructor
		~LocalToGlobalMapper() {};
};

/// The AssAdapter class combines tools to adapt the assemble routine
template <typename TAlgebra>
class AssAdapter
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of algebra value
		typedef typename vector_type::value_type value_type;

	public:
	/// constructor
		AssAdapter(): m_pBoolMarker(NULL), m_pSelector(NULL),
		m_bForceRegGrid(false), m_ConstraintTypesEnabled(CT_ALL),
		m_ElemTypesEnabled(EDT_ALL)
		{
			m_assIndex.index_set = false;
			m_pMapper = &m_pMapperCommon;
		}

	/// set local to global mapping
		void set_mapping(ILocalToGlobalMapper<TAlgebra>* pMapper = NULL)
		{
			if(pMapper){ m_pMapper = pMapper;}
			else{ m_pMapper = &m_pMapperCommon;}
		}

	/// LocalToGlobalMapper-function calls
		void AddLocalVec(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
		{ m_pMapper->AddLocalVec(vec, lvec, dd);}
		void AddLocalMatToGlobal(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
		{ m_pMapper->AddLocalMatToGlobal(mat, lmat, dd);}


	///	sets a marker to exclude elements from assembling
	/**
	 * This methods sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
		void set_marker(BoolMarker* mark = NULL){ m_pBoolMarker = mark; }

	///	sets a selector of elements for assembling
	/**
	 * This methods sets an element list. Only elements of this list will be
	 * assembled during assembling process. Especially the list defines the begin
	 * and end of the element-iterator in the element assembling-loop.
	 * If no element list is set, this corresponds to a assembling where the loop is
	 * carried out over all elements of a subset.
	 *
	 * \param[in]	sel		Selector
	 */
		void set_selector(Selector* sel = NULL){ m_pSelector = sel; }

	///	sets an index for which the assembling should be carried out
	/**
	 * This methods sets a boolean if an index-wise assemble routine should be used.
	 * This proceeding is e.g. useful for a nonlinear Gauss-Seidel or nonlinear
	 * Jacobi solver. The specific index is passed.
	 *
	 * \param[in]	ind			size_t
	 * \param[in]	index_set	bool
	 */
		void set_ass_index(){ set_ass_index(0, false);}
		void set_ass_index(size_t ind, bool index_set = true)
		{
			m_assIndex.index = ind; m_assIndex.index_set = index_set;
		}
	///	checks whether the assemble index is set or not
		size_t is_ass_index_set(){ return m_assIndex.index_set;}


	/// forces the assembling to consider the grid as regular
		void force_regular_grid(bool bForce) {m_bForceRegGrid = bForce;}

	///	returns if constraints enabled
		int constraints_enabled() const {return m_ConstraintTypesEnabled;}

	///	enables constraints
		void enable_constraints(int bEnableTypes) {m_ConstraintTypesEnabled = bEnableTypes;}

	///	returns type of boundary elem discs enabled
		int elem_discs_enabled() const {return m_ElemTypesEnabled;}

	///	enables boundary elem discs
		void enable_elem_discs(int bEnableTypes) {m_ElemTypesEnabled = bEnableTypes;}


	///	resize functions used in assemble funcs
		void resize(ConstSmartPtr<DoFDistribution> dd, vector_type& vec);
		void resize(ConstSmartPtr<DoFDistribution> dd, matrix_type& mat);

	///	gets the element iterator from the Selector
		template <typename TElem>
		void elemIter_fromSel(ConstSmartPtr<DoFDistribution> dd, int si,
				std::vector<TElem*>& elems);

	///	adapts the constraints in case of index-wise assembling
		template <typename TDomain>
		void adaptConstraint(SmartPtr<IDomainConstraint<TDomain, TAlgebra> >& constraint);

	///	only one index will be set to Dirichlet in case of index-wise assembling
	///	instead of setting a complete matrix row to Dirichlet
		void adjust_matrix(matrix_type& mat, const size_t index, const size_t alpha);
		void adjust_vector(vector_type& vec, const size_t index, const size_t alpha, double val);

	public:

	///	default LocalToGlobalMapper
		LocalToGlobalMapper<TAlgebra> 	m_pMapperCommon;
	///	LocalToGlobalMapper
		ILocalToGlobalMapper<TAlgebra>* m_pMapper;

	///	marker used to skip elements
		BoolMarker* m_pBoolMarker;

	///	selector used to set a list of elements for the assembling
		Selector* 	m_pSelector;

	///	index-wise assembling
		struct AssIndex{
			///	current index
			size_t index;
			bool index_set;
		};

	///	object for index-wise assemble routine
		AssIndex 	m_assIndex;

	/// forces the assembling to regard the grid as regular
		bool m_bForceRegGrid;

	///	enables the constraints
		int m_ConstraintTypesEnabled;

	///	enables the constraints
		int m_ElemTypesEnabled;
};

} // end namespace ug

#include "ass_adapter_impl.h"

#endif /* ASS_ADAPTER_H_ */
