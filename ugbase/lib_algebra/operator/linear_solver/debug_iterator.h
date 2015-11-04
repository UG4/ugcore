
#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__DEBUG_ITERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__DEBUG_ITERATOR__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/preconditioned_linear_operator_inverse.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"
#include "lib_algebra/operator/debug_writer.h"


#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{

/// Debugging iterator
/**
 * This class implements an iterator that can be used for debugging:
 *
 * It wraps for an ILinearIterator<typename TAlgebra::vector_type> object, which can be used as usual,
 * i.e., all call are forwarded. Moreover, this class can determine the algebarically smooth error. In order to doso, one must specify a solver.
 * This solver is called once during initialization.
 *
 * This class is a VectorDebugWritingObject <typename TAlgebra::vector_type>. If a DebugWriter is supplied using use the set_debug routines, the smooth error is written to file.
 *
 */
template <typename TAlgebra>
class DebugIterator :
public virtual ILinearIterator<typename TAlgebra::vector_type>,
public VectorDebugWritingObject <typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	// 	Base type
		typedef ILinearIterator<typename TAlgebra::vector_type> base_type;

		typedef IPreconditionedLinearOperatorInverse<vector_type> solver_type;
	private:
		typedef VectorDebugWritingObject <typename TAlgebra::vector_type> vdwo_type;

	public:
	///	Constructor
		DebugIterator() : from(0.0), to (1.0)
		{
			MatrixOperator<matrix_type, vector_type> *ptr =new MatrixOperator<matrix_type, vector_type>();
			m_pOperator  = SmartPtr<MatrixOperator<matrix_type, vector_type> >(ptr);
		};

		DebugIterator(SmartPtr<base_type> pprecond, SmartPtr<solver_type> psolver) : from(0.0), to (1.0)
		{
			set_preconditioner(pprecond);
			set_solver(psolver);
			MatrixOperator<matrix_type, vector_type> *ptr =new MatrixOperator<matrix_type, vector_type>();
			m_pOperator  = SmartPtr<MatrixOperator<matrix_type, vector_type> >(ptr);
		};

	/// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<DebugIterator<algebra_type> > newInst(new DebugIterator<algebra_type> ());
			newInst->set_damp(this->damping());
			newInst->set_preconditioner(this->get_preconditioner()->clone());  // clone preconditioner
			newInst->set_solver(this->get_solver());  						   // keep solver (initialized below)
			newInst->set_debug(this->vdwo_type::vector_debug_writer());
			newInst->from = this->from;
			newInst->to = this->to;
			return newInst;
		}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			if(m_pprecond.valid())
				return m_pprecond->supports_parallel();
			return true;
		}

	/// specify the real preconditioner (used for iterating)
		void set_preconditioner(SmartPtr<base_type> pprecond) {m_pprecond=pprecond;}
	/// specify the solver that will be used for debugging (optional)
		void set_solver(SmartPtr<solver_type> psolver) {m_solver=psolver;}
	/// specify bounds for random initial guess (optional)
		void set_random_bounds(double a, double b){from =a; to=b;}

	protected:
	/// get reference to 'real' preconditioner
		SmartPtr<base_type> get_preconditioner() {return m_pprecond;}
	/// get reference to aux. solver
		SmartPtr<solver_type> get_solver() {return m_solver;}


	///	name of preconditioner
		virtual const char* name() const
		{return "DebugIterator";}


	/// init (expensive)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J,
				                  const vector_type& u)
				{
				//	cast to matrix based operator
					SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
							J.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

				//	Check that matrix if of correct type
					if(pOp.invalid())
						UG_THROW(name() << "::init': Passed Operator is "
								"not based on matrix. This Preconditioner can only "
								"handle matrix-based operators.");

				//	forward request to matrix based implementation
					return init(pOp);
				}
	/// init (expensive)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L)
				{
				//	cast to matrix based operator
					SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
							L.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

				//	Check that matrix if of correct type
					if(pOp.invalid())
						UG_THROW(name() << "::init': Passed Operator is "
								"not based on matrix. This Preconditioner can only "
								"handle matrix-based operators.");

				//	forward request to matrix based implementation
					return init(pOp);
				}

		/// init (expensive, since it calls \sa find_smooth_error)
		bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			m_pOperator->get_matrix() = pOp->get_matrix();
			//m_pOperator->get_matrix().set_as_copy_of(pOp->get_matrix());
#ifdef UG_PARALLEL
			//m_pOperator->get_matrix().set_storage_Type(pOp->get_matrix()->get_storage_type());
#endif
			m_pprecond->init(pOp);

			//write_debug(pOp->get_matrix(), "DebugMatrix");
			find_smooth_error();

			return true;
		}

	/// Determines algebraically smooth error
	/** Solves Ae=0 starting from a random initial guess.
	 * \returns true, iff error was determined*/
		bool find_smooth_error()
		{
			if (m_solver.invalid())
				{
					UG_LOG("WARNING: cannot find smooth error; no solver supplied!");
					return false;
				}

			vector_type myRhs(m_pOperator->get_matrix().num_rows());
			vector_type myError(m_pOperator->get_matrix().num_rows());

			myError.set(0.0);
			myError.set_random(from, to);
			myRhs.set(0.0);

			this->write_debug(myError, "DebugIterError0");
			m_solver->init(m_pOperator);
			m_solver->apply(myError, myRhs);

			this->write_debug(myError, "DebugIterErrorS");

			return true;
		}

	/// forwarding call to original preconditioner
		virtual bool apply(vector_type& c, const vector_type& d)
		{
			return m_pprecond->apply(c, d);
		}

	/// forwarding call to original preconditioner
		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{
			return (m_pprecond->apply_update_defect(c, d));
		}

	protected:

		SmartPtr<base_type> m_pprecond;

		SmartPtr<solver_type> m_solver;
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_pOperator;

		double from, to;

};



} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFORMING__
