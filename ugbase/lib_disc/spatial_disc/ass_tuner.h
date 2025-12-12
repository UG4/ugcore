/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Raphael Prohl, Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER__

#include "lib_grid/tools/bool_marker.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug {

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
	CT_ASSEMBLED = 1,
	CT_MAY_DEPEND_ON_HANGING = 1 << 1,	// constraints which may depend on hanging DoFs; but NOT vice-versa
	CT_HANGING = 1 << 2,			// constraint defined for hanging DoFs
	CT_CONSTRAINTS = 1 << 3,		// any other constraint which MUST NOT depend on a hanging DoF
	CT_DIRICHLET = 1 << 4,			// Dirichlet constraints
	CT_ALL = CT_NONE | CT_ASSEMBLED | CT_MAY_DEPEND_ON_HANGING | CT_HANGING | CT_CONSTRAINTS | CT_DIRICHLET
};

template <typename TAlgebra>
class LocalToGlobalMapper
{
	public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Type of algebra matrix
	using matrix_type = typename algebra_type::matrix_type;

	///	Type of algebra vector
	using vector_type = typename algebra_type::vector_type;

	public:

	///	adds a local vector to the global one
		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec) const
			{ AddLocalVector(vec, lvec);}

	///	adds a local matrix to the global one
		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat) const
			{ AddLocalMatrixToGlobal(mat, lmat);}
};

/// The AssemblingTuner class combines tools to adapt the assembling routine.
template <typename TAlgebra>
class AssemblingTuner
{
	public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Type of algebra matrix
		using matrix_type = typename algebra_type::matrix_type;

	///	Type of algebra vector
		using vector_type = typename algebra_type::vector_type;

	///	Type of algebra value
		using value_type = typename vector_type::value_type;

	public:
	/// constructor
		AssemblingTuner(): m_pMapper(nullptr), m_pBoolMarker(nullptr), m_pSelector(nullptr),
		m_bSingleAssIndex(false), m_SingleAssIndex(0),
		m_bForceRegGrid(false), m_bModifySolutionImplemented(false),
		m_ConstraintTypesEnabled(CT_ALL), m_ElemTypesEnabled(EDT_ALL),
		m_bMatrixIsConst(false), m_bMatrixStructureIsConst(false), m_bClearOnResize(true) {}

	/// destructor
		virtual ~AssemblingTuner() = default;

	/// set local to global mapping
		void set_mapping(ILocalToGlobalMapper<TAlgebra>* pMapper = nullptr)
		{
			m_pMapper = pMapper;
		}

	/// LocalToGlobalMapper-function calls
		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
		                 ConstSmartPtr<DoFDistribution> dd) const
		{
			if (m_pMapper)
				m_pMapper->add_local_vec_to_global(vec, lvec, dd);
			else
				m_defaultMapper.add_local_vec_to_global(vec, lvec);
		}

		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
		                         ConstSmartPtr<DoFDistribution> dd) const
		{
			if (m_pMapper)
				m_pMapper->add_local_mat_to_global(mat, lmat, dd);
			else
				m_defaultMapper.add_local_mat_to_global(mat, lmat);
		}

		void modify_LocalSol(LocalVector& vecMod, const LocalVector& lvec,
		                         ConstSmartPtr<DoFDistribution> dd) const
		{
			if (m_pMapper)
				m_pMapper->modify_LocalSol(vecMod, lvec, dd);
		}
	///	sets a marker to exclude elements from assembling
	/**
	 * This methods sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
		void set_marker(BoolMarker* mark = nullptr){ m_pBoolMarker = mark; }

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
		void set_selector(Selector* sel = nullptr){ m_pSelector = sel; }

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
		virtual void set_force_regular_grid(bool bForce) {m_bForceRegGrid = bForce;}

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
		bool selected_elements_used() const {return (m_pSelector != nullptr);}

	///	returns if element is to be used in assembling
		template <typename TElem>
		bool element_used(TElem* elem) const;

	///	only one index will be set to Dirichlet in case of index-wise assembling
	///	instead of setting a complete matrix row to Dirichlet
		void set_dirichlet_row(matrix_type& mat, const DoFIndex& ind) const;
		void set_dirichlet_val(vector_type& vec, const DoFIndex& ind, const double val) const;

	/// Disable clearing of matrix/vector when resizing.
	/// This is useful when an IAssemble object consists of more than one
	/// domain disc, e.g., CompositeTimeDisc.
		void disable_clear_on_resize() {m_bClearOnResize = false;}

	/**
	 * specify whether matrix will be modified by assembling
	 * disables matrix assembling if set to true
	 *
	 * @param bCh set true if matrix is not to be changed during assembling
	 */
		void set_matrix_is_const(bool bCh) {m_bMatrixIsConst = bCh;}

		void set_matrix_structure_is_const(bool b) {m_bMatrixStructureIsConst = b;}

	/**
	 * whether matrix is to be modified by assembling
	 *
	 * @return true iff matrix is not to be modified
	 */
		bool matrix_is_const() const {return m_bMatrixIsConst;}

	protected:
	///	default LocalToGlobalMapper
		LocalToGlobalMapper<TAlgebra> m_defaultMapper;

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

	/// keeps matrix structure from last call if set to true
		bool m_bMatrixStructureIsConst;

	/// disables clearing of vector/matrix on resize
		bool m_bClearOnResize;
};

} // end namespace ug

#include "ass_tuner_impl.h"

#endif