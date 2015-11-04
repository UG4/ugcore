/*
 * ass_tuner.h
 *
 *  Created on: 04.02.2013
 *      Author: raphaelprohl, Andreas Vogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER__

#include "lib_grid/tools/bool_marker.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class IDomainConstraint;

/// Types of constraint
/**
 * These types control the order the constraints are applied in.
 * Constraints with a lower number will be applied first. The order of
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
		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
			{ AddLocalVector(vec,lvec);}

	///	adds a local matrix to the global one
		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
			{ AddLocalMatrixToGlobal(mat,lmat);}

	///	modifies local solution vector for adapted defect computation
		void modify_LocalSol(LocalVector& vecMod, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd){};

	///	destructor
		~LocalToGlobalMapper() {};
};

/// The AssemblingTuner class combines tools to adapt the assembling routine.
template <typename TAlgebra>
class AssemblingTuner
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
		AssemblingTuner(): m_pBoolMarker(NULL), m_pSelector(NULL),
		m_bSingleAssIndex(false), m_SingleAssIndex(0),
		m_bForceRegGrid(false), m_bModifySolutionImplemented(false),
		m_ConstraintTypesEnabled(CT_ALL), m_ElemTypesEnabled(EDT_ALL),
		m_bMatrixIsConst(false)
		{
			m_pMapper = &m_pMapperCommon;
		}

	/// set local to global mapping
		void set_mapping(ILocalToGlobalMapper<TAlgebra>* pMapper = NULL)
		{
			if(pMapper)
				m_pMapper = pMapper;
			else
				m_pMapper = &m_pMapperCommon;
		}

	/// LocalToGlobalMapper-function calls
		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
		                 ConstSmartPtr<DoFDistribution> dd) const
		{ m_pMapper->add_local_vec_to_global(vec, lvec, dd);}

		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
		                         ConstSmartPtr<DoFDistribution> dd) const
		{ m_pMapper->add_local_mat_to_global(mat, lmat, dd);}

		void modify_LocalSol(LocalVector& vecMod, const LocalVector& lvec,
		                         ConstSmartPtr<DoFDistribution> dd) const
		{ m_pMapper->modify_LocalSol(vecMod, lvec, dd);}
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
	 * assembled during assembling process. The list especially defines the begin
	 * and end of the element-iterator in the element assembling-loop.
	 * If no element list is set, this corresponds to an assembling where the loop is
	 * carried out over all elements of a subset.
	 *
	 * \param[in]	sel		Selector
	 */
		void set_selector(Selector* sel = NULL){ m_pSelector = sel; }

	///	sets an index for which the assembling should be carried out
	/**
	 * This methods sets a boolean if an DoFindex-wise assemble routine should be used.
	 * This proceeding is e.g. useful for a nonlinear Gauss-Seidel or nonlinear
	 * Jacobi solver. The specific index is passed.
	 *
	 * \param[in]	ind			DoFIndex
	 */
		void disable_single_index_assembling() {m_bSingleAssIndex = false;}
		void set_single_index_assembling(const size_t index)
		{
			m_SingleAssIndex = index; m_bSingleAssIndex = true;
		}

	///	checks whether the assemble DoFindex is set or not
		bool single_index_assembling_enabled() const {return m_bSingleAssIndex;}


	///	enables the usage of modify solution
		void enable_modify_solution(bool bEnable) {m_bModifySolutionImplemented = bEnable;}

	///	checks whether the assemble index is set or not
		bool modify_solution_enabled() const {return m_bModifySolutionImplemented;}


	/// forces the assembling to consider the grid as regular
		void set_force_regular_grid(bool bForce) {m_bForceRegGrid = bForce;}

	/// returns if assembling is to considered as regular grid
		bool regular_grid_forced() const {return m_bForceRegGrid;}


	///	enables constraints
		void enable_constraints(int bEnableTypes) {m_ConstraintTypesEnabled = bEnableTypes;}

	///	returns flags of enabled constraints
		int enabled_constraints() const {return m_ConstraintTypesEnabled;}

	///	returns if constraint type enabled
		bool constraint_type_enabled(int type) const {return (type & m_ConstraintTypesEnabled);}


	///	enables elem discs
		void enable_elem_discs(int bEnableTypes) {m_ElemTypesEnabled = bEnableTypes;}

	///	returns flags of enabled elem discs
		int enabled_elem_discs() const {return m_ElemTypesEnabled;}

	///	returns if elem disc type enabled
		bool elem_disc_type_enabled(int type) const {return (type & m_ElemTypesEnabled);}


	///	resize functions used in assemble funcs
		void resize(ConstSmartPtr<DoFDistribution> dd, vector_type& vec) const;
		void resize(ConstSmartPtr<DoFDistribution> dd, matrix_type& mat) const;

	///	gets the element iterator from the Selector
		template <typename TElem>
		void collect_selected_elements(std::vector<TElem*>& vElem, ConstSmartPtr<DoFDistribution> dd, int si) const;

	///	returns if only selected elements used for assembling
		bool selected_elements_used() const {return (m_pSelector != NULL);}

	///	returns if element is to be used in assembling
		template <typename TElem>
		bool element_used(TElem* elem) const;

	///	only one index will be set to Dirichlet in case of index-wise assembling
	///	instead of setting a complete matrix row to Dirichlet
		void set_dirichlet_row(matrix_type& mat, const DoFIndex& ind) const;
		void set_dirichlet_val(vector_type& vec, const DoFIndex& ind, const double val) const;

	/**
	 * specify whether matrix will be modified by assembling
	 * disables matrix assembling if set to true
	 *
	 * @param bCh set true if matrix is not to be changed during assembling
	 */
		void set_matrix_is_const(bool bCh) {m_bMatrixIsConst = bCh;}

	/**
	 * whether matrix is to be modified by assembling
	 *
	 * @return true iff matrix is not to be modified
	 */
		bool matrix_is_const() const {return m_bMatrixIsConst;}

	protected:
	///	default LocalToGlobalMapper
		LocalToGlobalMapper<TAlgebra> m_pMapperCommon;

	///	LocalToGlobalMapper
		ILocalToGlobalMapper<TAlgebra>* m_pMapper;

	///	marker used to skip elements
		BoolMarker* m_pBoolMarker;

	///	selector used to set a list of elements for the assembling
		Selector* m_pSelector;

	///	object for DoFindex-wise assemble routine
		bool m_bSingleAssIndex;
		size_t m_SingleAssIndex;

	/// forces the assembling to regard the grid as regular
		bool m_bForceRegGrid;

	/// calls the 'modify_solution()' method of constraints;
	///	gives the modified solution to the assembling methods
		bool m_bModifySolutionImplemented;

	///	enables the constraints
		int m_ConstraintTypesEnabled;

	///	enables the constraints
		int m_ElemTypesEnabled;

	/// disables matrix assembling if set to false
		bool m_bMatrixIsConst;
};

} // end namespace ug

#include "ass_tuner_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER__ */
