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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

#include "lib_disc/assemble_interface.h"

namespace ug{

///	matrix operator based on the assembling of a problem
/**
 * This operator implements the MatrixOperator interface, thus is basically a
 * matrix that can be applied to vectors. In addition the class allows to set
 * an IAssemble object and an the GridLevel. Invoking the init method the
 * matrix is created using the IAssemble routine on the given GridLevel.
 *
 * \tparam	TAlgebra			algebra type
 */
template <typename TAlgebra>
class AssembledLinearOperator :
	public virtual MatrixOperator<	typename TAlgebra::matrix_type,
									typename TAlgebra::vector_type>
{
	public:
	///	Type of Algebra
		using algebra_type = TAlgebra;

	///	Type of Vector
		using vector_type = typename TAlgebra::vector_type;

	///	Type of Matrix
		using matrix_type = typename TAlgebra::matrix_type;

	///	Type of base class
		using base_type = MatrixOperator<matrix_type,vector_type>;

	public:
	///	Default Constructor
		AssembledLinearOperator() :	m_spAss(nullptr) {};

	///	Constructor
		AssembledLinearOperator(SmartPtr<IAssemble<TAlgebra> > ass) : m_spAss(ass) {};

	///	Constructor
		AssembledLinearOperator(SmartPtr<IAssemble<TAlgebra> > ass, const GridLevel& gl)
			: m_spAss(ass), m_gridLevel(gl) {};

	///	sets the discretization to be used
		void set_discretization(SmartPtr<IAssemble<TAlgebra> > ass) {m_spAss = ass;}

	///	returns the discretization to be used
		SmartPtr<IAssemble<TAlgebra> > discretization() {return m_spAss;}

	///	sets the level used for assembling
		void set_level(const GridLevel& gl) {m_gridLevel = gl;}

	///	returns the level
		const GridLevel& level() const {return m_gridLevel;}

	///	initializes the operator that may depend on the current solution
		void init(const vector_type& u) override;

	///	initialize the operator
		void init() override;

	///	initializes the operator and assembles the passed rhs vector
		void init_op_and_rhs(vector_type& b);

	///	compute d = J(u)*c (here, J(u) is a Matrix)
		void apply(vector_type& d, const vector_type& c) override;

	///	Compute d := d - J(u)*c
		void apply_sub(vector_type& d, const vector_type& c) override;

	///	Set Dirichlet values
		void set_dirichlet_values(vector_type& u);

	///	Destructor
		~AssembledLinearOperator() override = default;

	protected:
	// 	assembling procedure
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

	// 	DoF Distribution used
		GridLevel m_gridLevel;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// help function to assemble a linear operator
/**
 * This function initializes the operator, sets the vector b to the computed rhs
 * and sets the dirichlet post processes for the vector u.
 *
 * \param[out]	op		Operator
 * \param[out]	u		Solution
 * \param[out]	b		Rigth-Hand side vector
 *
 * \tparam	TAlgebra			algebra type
 */
template <typename TAlgebra>
void AssembleLinearOperatorRhsAndSolution
		(AssembledLinearOperator<TAlgebra>& op,
		 typename TAlgebra::vector_type& u,
		 typename TAlgebra::vector_type& b);

} // namespace ug

// include implementation
#include "assembled_linear_operator_impl.h"

#endif