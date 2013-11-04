/*
 * sor.h
 *
 *  Created on: 9.09.2012
 *      Author: Arne Naegel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__OVERRELAXATION__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__OVERRELAXATION__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_vector_impl.h"
#endif

namespace ug{

// It would be nicer to implement this object as an
// IPreconditioner<TAlgebra>,  but IPreconditioner is not an interface...

template <typename TAlgebra>
class SOR : public GaussSeidel<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef GaussSeidel<TAlgebra> base_type;

	protected:
		//using base_type::set_debug;
		using base_type::debug_writer;
		//using base_type::write_debug;

	public:
	//	Constructor
		SOR() {};

		// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<SOR<algebra_type> > newInst(new SOR<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

		// init
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J,
					                  const vector_type& u)
		{
			m_u = &u;

			return true;
		}

	protected:
		//	Name of preconditioner
		virtual const char* name() const {return "SOR<GaussSeidel>";}


		//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			if (base_type::step(pOp, c, d))
			{
				VecScaleAdd(c, 1.0, c, -1.0, *m_u);
				return true;
			}

			return false;
		}

	protected:
// base_type m_basePrecond;
		const vector_type *m_u;  // WARNING: This could be a reference to a smart ptr!
};



} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SOR__
