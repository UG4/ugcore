/*
 * \file external_solver.h
 *
 *
 * \date 16.01.2014
 * \author Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__EXTERNAL_SOLVER_
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__EXTERNAL_SOLVER_

#include "common/common.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"


#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "common/progress.h"
#include "common/util/ostream_util.h"

#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/cpu_algebra_types.h"


namespace ug{

class IExternalSolverImplementation
{
public:
	virtual void init(const CPUAlgebra::matrix_type &A) = 0;
	virtual void apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d) = 0;
	virtual const char* name() const = 0;
	virtual ~IExternalSolverImplementation() {}
};

template <typename TAlgebra>
class IExternalSolver
		: public IMatrixOperatorInverse<typename TAlgebra::matrix_type,
			  	  	  	  	  	  	  	    typename TAlgebra::vector_type>
{
	protected:
		IExternalSolverImplementation *impl;


	public:
		virtual const char* name() const
		{
			return impl->name();
		}

	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		IExternalSolver()
		{
		};

	// 	Clone

		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			UG_THROW("");
			return NULL;
		}


	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return false;}


	public:



		void mat_preprocess(const matrix_type &A)
		{
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

			CPUAlgebra::matrix_type mat;

			#ifdef UG_PARALLEL
				mat.set_storage_type(PST_ADDITIVE);
				mat.set_layouts(CreateLocalAlgebraLayouts());
			#endif

			m_size = GetDoubleSparseFromBlockSparse(mat, A);

			impl->init(mat);
		}

		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spOperator;
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
		// 	remember operator
			m_spOperator = Op;

		//	init LU operator
			mat_preprocess(m_spOperator->get_matrix());
		//	we're done
			return true;
		}


	protected:

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type &A = *pOp;
			mat_preprocess(A);

			return true;
		}

		void get_vector(CPUAlgebra::vector_type &v_scalar, const vector_type& v)
		{
			for(size_t i=0, k=0; i<v.size(); i++)
			{
				for(size_t j=0; j<GetSize(v[i]); j++)
					v_scalar[k++] = BlockRef(v[i],j);
			}
		}

		void set_vector(CPUAlgebra::vector_type &v_scalar, vector_type& v)
		{
			for(size_t i=0, k=0; i<v.size(); i++)
			{
				for(size_t j=0; j<GetSize(v[i]); j++)
					BlockRef(v[i],j) = v_scalar[k++];
			}
		}

	//	Stepping routine
		virtual bool apply(vector_type& c, const vector_type& d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);
#endif
			get_vector(m_d, d);
			m_c.set(0.0);


			impl->apply(m_c, m_d);

			set_vector(m_c, c);

#ifdef 	UG_PARALLEL
		//	Correction is always consistent
		//	todo: We set here correction to consistent, but it is not. Think about how to use ilu in parallel.
			c.set_storage_type(PST_CONSISTENT);
#endif
			return true;
		}

		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
		//	solve u
			if(!apply(u, f)) return false;

		//	update defect
			if(!m_spOperator->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LU::apply_return_defect': "
						"Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}



	protected:
	//	Postprocess routine
		virtual bool postprocess() {return true;}

	CPUAlgebra::vector_type m_c, m_d;
	size_t m_size;
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__EXTERNAL_SOLVER_ */
