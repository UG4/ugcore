/*
 * schur_complement_inverse.h
 *
 *  Created on: 08.01.2014
 *      Author: mrupp
 */

#ifndef SCHUR_COMPLEMENT_INVERSE_H_
#define SCHUR_COMPLEMENT_INVERSE_H_

#ifdef UG_PARALLEL


#include "schur.h"
#include "slicing.h"
#include "lib_algebra/operator/linear_solver/auto_linear_solver.h"

namespace ug{



template<typename TAlgebra>
class SchurInverseWithOperator : public ISchurComplementInverse<TAlgebra>
{
public:
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	SchurInverseWithOperator(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	virtual bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op)
	{
		op->set_skeleton_debug(m_linOpInv);
		return m_linOpInv->init(op);
	}

	virtual bool apply(vector_type& u, const vector_type& f)
	{
		return m_linOpInv->apply(u, f);
	}

	virtual bool apply_return_defect(vector_type& u, vector_type& f)
	{
		return m_linOpInv->apply_return_defect(u, f);
	}

	virtual std::string config_string() const
	{
		std::stringstream ss; ss << "SchurInverseWithOperator\n";
		ss << " Solver: " << ConfigShift(m_linOpInv->config_string()) << "\n";
		return ss.str();
	}
	virtual bool supports_parallel() const
	{
		return m_linOpInv->supports_parallel();
	}

protected:
	SmartPtr<ILinearOperatorInverse<vector_type> > m_linOpInv;

};


template<typename TAlgebra>
class SchurInverseWithFullMatrix : public ISchurComplementInverse<TAlgebra>
{
public:
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	SchurInverseWithFullMatrix(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	virtual bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op)
	{
		m_exactSchurOp = make_sp(new MatrixOperator<matrix_type, vector_type>);
		op->compute_matrix(m_exactSchurOp->get_matrix());
		op->set_skeleton_debug(m_linOpInv);
		return m_linOpInv->init(m_exactSchurOp);
	}

	virtual bool apply(vector_type& u, const vector_type& f)
	{
		return m_linOpInv->apply(u, f);
	}

	virtual bool apply_return_defect(vector_type& u, vector_type& f)
	{
		return m_linOpInv->apply_return_defect(u, f);
	}

	virtual std::string config_string() const
	{
		std::stringstream ss; ss << "SchurInverseWithFullMatrix\n";
		ss << " Solver: " << ConfigShift(m_linOpInv->config_string()) << "\n";
		return ss.str();
	}
	virtual bool supports_parallel() const
	{
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
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	SchurInverseWithAGammaGamma(SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > linSolver)
	{
		UG_COND_THROW(!linSolver.valid(), "?");
		m_linSolver = linSolver;
	}

	virtual bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op)
	{
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
	virtual bool apply(vector_type& u, const vector_type& f)
	{
		return m_linSolver->apply(u, f);
	}

	virtual bool apply_return_defect(vector_type& u, vector_type& f)
	{
		return m_linSolver->apply_return_defect(u, f);
	}
	virtual std::string config_string() const
	{
		std::stringstream ss; ss << "SchurInverseWithAGammaGamma\n";
		ss << " Solver: " << ConfigShift(m_linSolver->config_string()) << "\n";
		return ss.str();
	}
	virtual bool supports_parallel() const
	{
		return m_linSolver->supports_parallel();
	}
protected:
	SmartPtr<IPreconditionedLinearOperatorInverse<vector_type> > m_linSolver;
};



template <typename TAlgebra, typename M, typename X, typename Y = X>
class SchurComplementMatrixOperator : public MatrixOperator<M, X, Y>, public UpdateableMatrixOperator
{
	typedef M matrix_type;

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
	virtual void init(const X& u) { init(); }

// 	Init Operator L
	virtual void init()
	{
		if(invalid)
		{
			m_op->compute_matrix(get_matrix());
			invalid = false;
		}
	}

	virtual void calculate_matrix()
	{
		init();
	}

// 	Apply Operator f = L*u (e.g. d = J(u)*c in iterative scheme)
	virtual void apply(Y& f, const X& u) {m_op->apply(f,u);}

// 	Apply Operator, i.e. f = f - L*u;
	virtual void apply_sub(Y& f, const X& u) {m_op->apply_sub(f,u);}

// 	Access to matrix
	virtual M& get_matrix() {return *this;};
};

// not completely working at the moment
template<typename TAlgebra>
class SchurInverseWithAutoFullMatrix : public ISchurComplementInverse<TAlgebra>
{
public:
	typedef TAlgebra algebra_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	SchurInverseWithAutoFullMatrix(SmartPtr<ILinearOperatorInverse<vector_type> > linOpInv )
	{
		m_linOpInv = linOpInv;
	}

	virtual bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op)
	{
		if(m_exactSchurOp.valid() == false)
			m_exactSchurOp = make_sp(new SchurComplementMatrixOperator<TAlgebra, matrix_type, vector_type>(op));
		else
			m_exactSchurOp->set_op(op);
		return m_linOpInv->init(m_exactSchurOp);
	}

	virtual bool apply(vector_type& u, const vector_type& f)
	{
		return m_linOpInv->apply(u, f);
	}

	virtual bool apply_return_defect(vector_type& u, vector_type& f)
	{
		return m_linOpInv->apply_return_defect(u, f);
	}

	virtual std::string config_string() const
	{
		std::stringstream ss; ss << "SchurInverseWithAutoFullMatrix\n";
		ss << " Solver: " << ConfigShift(m_linOpInv->config_string()) << "\n";
		return ss.str();
	}
	virtual bool supports_parallel() const
	{
		return m_linOpInv->supports_parallel();
	}

protected:
	SmartPtr<SchurComplementMatrixOperator<TAlgebra, matrix_type, vector_type> > m_exactSchurOp;
	SmartPtr<ILinearOperatorInverse<vector_type> > m_linOpInv;

};

}

#endif /* UG_PARALLEL */
#endif /* SCHUR_COMPLEMENT_INVERSE_H_ */
