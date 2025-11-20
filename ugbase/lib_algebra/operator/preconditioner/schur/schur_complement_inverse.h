/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef SCHUR_COMPLEMENT_INVERSE_H_
#define SCHUR_COMPLEMENT_INVERSE_H_

#ifdef UG_PARALLEL



#include "lib_algebra/operator/linear_solver/auto_linear_solver.h"


#include "schur.h"
#include "schur_complement_inverse_interface.h"


namespace ug{



template<typename TAlgebra>
class SchurInverseWithOperator : public ISchurComplementInverse<TAlgebra>
{
public:
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;

	SchurInverseWithOperator(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op) override {
		op->set_skeleton_debug(m_linOpInv);
		return m_linOpInv->init(op);
	}

	bool apply(vector_type& u, const vector_type& f) override {
		return m_linOpInv->apply(u, f);
	}

	bool apply_return_defect(vector_type& u, vector_type& f) override {
		return m_linOpInv->apply_return_defect(u, f);
	}

	std::string config_string() const override {
		std::stringstream ss; ss << "SchurInverseWithOperator\n";
		ss << " Solver: " << ConfigShift(m_linOpInv->config_string()) << "\n";
		return ss.str();
	}

	bool supports_parallel() const override {
		return m_linOpInv->supports_parallel();
	}

protected:
	SmartPtr<ILinearOperatorInverse<vector_type> > m_linOpInv;

};


template<typename TAlgebra>
class SchurInverseWithFullMatrix : public ISchurComplementInverse<TAlgebra>
{
public:
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;

	SchurInverseWithFullMatrix(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op) override {
		PROFILE_BEGIN(SchurInverseWithFullMatrix_init)

		m_exactSchurOp = make_sp(new MatrixOperator<matrix_type, vector_type>);

		PROFILE_BEGIN(SchurInverseWithFullMatrix_compute_matrix)
			op->compute_matrix(m_exactSchurOp->get_matrix());
		PROFILE_END();

		op->set_skeleton_debug(m_linOpInv);

		PROFILE_BEGIN(SchurInverseWithFullMatrix_init_solver)
		return m_linOpInv->init(m_exactSchurOp);
	}

	bool apply(vector_type& u, const vector_type& f) override {
		PROFILE_BEGIN(SchurInverseWithFullMatrix_apply)
		return m_linOpInv->apply(u, f);
	}

	bool apply_return_defect(vector_type& u, vector_type& f) override {
		PROFILE_BEGIN(SchurInverseWithFullMatrix_apply_return_defect)
		return m_linOpInv->apply_return_defect(u, f);
	}

	std::string config_string() const override {
		std::stringstream ss; ss << "SchurInverseWithFullMatrix\n";
		ss << " Solver: " << ConfigShift(m_linOpInv->config_string()) << "\n";
		return ss.str();
	}

	bool supports_parallel() const override {
		return m_linOpInv->supports_parallel();
	}

protected:
	SmartPtr<MatrixOperator<matrix_type, vector_type> > m_exactSchurOp;
	SmartPtr<ILinearOperatorInverse<vector_type> > m_linOpInv;

};


template<typename TAlgebra>
class SchurInverseWithAGammaGamma : public ISchurComplementInverse<TAlgebra>
{
public:
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;

	SchurInverseWithAGammaGamma(SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > linSolver)
	{
		UG_COND_THROW(!linSolver.valid(), "?");
		m_linSolver = linSolver;
	}

	bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op) override {
		PROFILE_BEGIN(SchurInverseWithAGammaGamma_init)

		op->set_skeleton_debug(m_linSolver);
		SmartPtr<IPreconditioner<TAlgebra> > precond =
				m_linSolver->preconditioner().template cast_dynamic< IPreconditioner<TAlgebra> > ();
		UG_COND_THROW(!precond.valid(), "?");
		precond->set_approximation(op->sub_operator(SD_SKELETON, SD_SKELETON));
//		precond->set_damp(0.5);
//		m_linSolver->set_preconditioner(precond);
		return m_linSolver->init(op);
		// do init of linsolver before or after set approximation?
	}

	bool apply(vector_type& u, const vector_type& f) override {
		PROFILE_BEGIN(SchurInverseWithAGammaGamma_apply)
		return m_linSolver->apply(u, f);
	}

	bool apply_return_defect(vector_type& u, vector_type& f) override {
		PROFILE_BEGIN(SchurInverseWithAGammaGamma_apply_return_defect)
		return m_linSolver->apply_return_defect(u, f);
	}

	std::string config_string() const override {
		std::stringstream ss; ss << "SchurInverseWithAGammaGamma\n";
		ss << " Solver: " << ConfigShift(m_linSolver->config_string()) << "\n";
		return ss.str();
	}

	bool supports_parallel() const override {
		return m_linSolver->supports_parallel();
	}
protected:
	SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > m_linSolver;
};



template <typename TAlgebra, typename M, typename X, typename Y = X>
class SchurComplementMatrixOperator : public MatrixOperator<M, X, Y>, public UpdateableMatrixOperator
{
	using matrix_type = M;

	SmartPtr<SchurComplementOperator<TAlgebra> > m_op;
public:
	bool invalid;
	SchurComplementMatrixOperator(SmartPtr<SchurComplementOperator<TAlgebra> > op)
	{
		set_op(op);
	}

	void set_op(SmartPtr<SchurComplementOperator<TAlgebra> > op)
	{
		m_op = op;
		invalid = true;
	}

// 	Init Operator J(u)
	void init(const X& u) override { init(); }

// 	Init Operator L
	void init() override {
		if(invalid)
		{
			m_op->compute_matrix(get_matrix());
			invalid = false;
		}
	}

	void calculate_matrix() override {
		init();
	}

// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
	void apply(Y& f, const X& u) override {m_op->apply(f,u);}

// 	Apply Operator, i.e. f = f - L*u;
	void apply_sub(Y& f, const X& u) override {m_op->apply_sub(f,u);}

// 	Access to matrix
	M& get_matrix() override {return *this;};
};

// not completely working at the moment
template<typename TAlgebra>
class SchurInverseWithAutoFullMatrix : public ISchurComplementInverse<TAlgebra>
{
public:
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;

	SchurInverseWithAutoFullMatrix(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op) override {
		if(m_exactSchurOp.valid() == false)
			m_exactSchurOp = make_sp(new SchurComplementMatrixOperator<TAlgebra, matrix_type, vector_type>(op));
		else
			m_exactSchurOp->set_op(op);
		return m_linOpInv->init(m_exactSchurOp);
	}

	bool apply(vector_type& u, const vector_type& f) override {
		return m_linOpInv->apply(u, f);
	}

	bool apply_return_defect(vector_type& u, vector_type& f) override {
		return m_linOpInv->apply_return_defect(u, f);
	}

	std::string config_string() const override {
		std::stringstream ss; ss << "SchurInverseWithAutoFullMatrix\n";
		ss << " Solver: " << ConfigShift(m_linOpInv->config_string()) << "\n";
		return ss.str();
	}

	bool supports_parallel() const override {
		return m_linOpInv->supports_parallel();
	}

protected:
	SmartPtr<SchurComplementMatrixOperator<TAlgebra, matrix_type, vector_type> > m_exactSchurOp;
	SmartPtr<ILinearOperatorInverse<vector_type> > m_linOpInv;

};

}

#endif
#endif