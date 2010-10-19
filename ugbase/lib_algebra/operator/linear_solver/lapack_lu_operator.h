/*
 * amg_solver.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__

#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/martin_algebra/lapack_lu.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TAlgebra>
class LapackLUSolver : public IMatrixOperatorInverse<	typename TAlgebra::vector_type,
														typename TAlgebra::vector_type,
														typename TAlgebra::matrix_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
		LapackLUSolver() :
			m_pOperator(NULL), m_lapacklu()
		{};

	//	set operator L, that will be inverted
		virtual bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& Op)
		{
		// 	remember operator
			m_pOperator = &Op;

		//	get matrix of Operator
			m_pMatrix = &m_pOperator->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
				{UG_LOG("ERROR in LapackLUOperator::init: No Matrix given,\n"); return false;}

		//	init lapack operator
			if(!m_lapacklu.init(*m_pMatrix))
				{UG_LOG("ERROR in LapackLUOperator::init: Cannot init LapackLU.\n"); return false;}

		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f
		virtual bool apply(vector_type& u, const vector_type& f)
		{
#ifdef UG_PARALLEL
			if(!f.has_storage_type(PST_ADDITIVE))
			{
				UG_LOG("ERROR: In 'LaplackLUSolver::apply':Inadequate storage format of Vector f.\n");
				return false;
			}
			if(!u.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR: In 'LaplackLUSolver::apply':Inadequate storage format of Vector u.\n");
				return false;
			}
#endif
			UG_ASSERT(f.size() == m_pMatrix->num_rows(),	"Vector and Row sizes have to match!");
			UG_ASSERT(u.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(f.size() == u.size(), "Vector sizes have to match!");

			// TODO: This must be inverted
			if(!m_lapacklu.apply(f, u))
				{UG_LOG("ERROR in LapackLUOperator::apply: Cannot init LapackLU.\n"); return false;}

		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f AND return defect f := f - L*u
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
		//	solve u
			if(!apply(u, f)) return false;

		//	update defect
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LapackLUSolver::apply_return_defect': Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

	// 	Destructor
		virtual ~LapackLUSolver() {};

	protected:
		// Operator to invert
		IMatrixOperator<vector_type, vector_type, matrix_type>* m_pOperator;

		// matrix to invert
		matrix_type* m_pMatrix;

		// lapack inverse
		LapackLU m_lapacklu;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
