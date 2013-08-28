/*
 * \file ilut.h
 *
 * \date 20.07.2013
 * \author Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT_SCALAR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT_SCALAR___

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "common/progress.h"
#include "common/util/ostream_util.h"

#include "lib_algebra/algebra_common/vector_util.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "ilut.h"

namespace ug{

template <typename TAlgebra>
class ILUTScalarPreconditioner : public IPreconditioner<TAlgebra>
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

	private:
		typedef typename matrix_type::value_type block_type;
		using IPreconditioner<TAlgebra>::debug_writer;
		using IPreconditioner<TAlgebra>::set_debug;

	public:
	//	Constructor
		ILUTScalarPreconditioner(double eps=1e-6) :
			m_eps(eps), m_info(false)
		{};

	// 	Clone

		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ILUTScalarPreconditioner<algebra_type> > newInst(new ILUTScalarPreconditioner<algebra_type>(m_eps));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_info(m_info);
			return newInst;
		}

	// 	Destructor
		virtual ~ILUTScalarPreconditioner()
		{
		};

	///	sets threshold for incomplete LU factorisation (added 01122010ih)
		void set_threshold(number thresh)
		{
			m_eps = thresh;
		}
		
	///	sets storage information output to true or false
		void set_info(bool info)
		{
			m_info = info;
		}
		
	public:
		void preprocess(const matrix_type &A)
		{
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

			const size_t nrOfRows = block_traits<typename matrix_type::value_type>::static_num_rows;
			UG_ASSERT(nrOfRows == block_traits<typename matrix_type::value_type>::static_num_cols, "only square matrices supported");
			m_size = A.num_rows() * nrOfRows;

			ilut = new ILUTPreconditioner<CPUAlgebra>(m_eps);
			ilut->set_info(m_info);

			mo = new MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type>;
			CPUAlgebra::matrix_type &mat = mo->get_matrix();
			mat.resize_and_clear(m_size, m_size);
			#ifdef UG_PARALLEL
				mat.set_storage_type(PST_ADDITIVE);
			#endif

			for(size_t r=0; r<A.num_rows(); r++)
				for(typename matrix_type::const_row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
				{
					size_t rr = r*nrOfRows;
					size_t cc = it.index()*nrOfRows;
					for(size_t r2=0; r2<nrOfRows; r2++)
						for(size_t c2=0; c2<nrOfRows; c2++)
							mat(rr + r2, cc + c2) = BlockRef(it.value(), r2, c2);
				}
			mat.defragment();

			SmartPtr<StdConvCheck<CPUAlgebra::vector_type> > convCheck = new StdConvCheck<CPUAlgebra::vector_type>;
			convCheck->set_maximum_steps(100);
			convCheck->set_minimum_defect(1e-50);
			convCheck->set_reduction(1e-20);
			convCheck->set_verbose(false);

			linearSolver.set_preconditioner(ilut);
			linearSolver.set_convergence_check(convCheck);
			linearSolver.init(mo);
		}
	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUTScalar";}


	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type &A = *pOp;
			preprocess(A);

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
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);
#endif
			get_vector(m_d, d);
			m_c.set(0.0);

			ilut->apply(m_c, m_d);

			set_vector(m_c, c);

#ifdef 	UG_PARALLEL
		//	Correction is always consistent
		//	todo: We set here correction to consistent, but it is not. Think about how to use ilu in parallel.
			c.set_storage_type(PST_CONSISTENT);
#endif
			return true;
		}

	public:
		bool solve(vector_type &c, const vector_type &d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);

			/*AlgebraLayouts *p = new AlgebraLayouts;
			p->proc_comm() = pcl::ProcessCommunicator(pcl::PCD_LOCAL);
			SmartPtr<AlgebraLayouts> layouts(p);
			m_c.set_layouts(c.layouts());
			m_d.set_layouts(d.layouts());*/
#endif
			get_vector(m_d, d);
			m_c.set(0.0);

			//ApplyLinearSolver(mo, m_u, m_b, linearSolver);
			//linearSolver.apply_return_defect(m_c, m_d);
			ilut->apply(m_c, m_d);


			set_vector(m_c, c);
			return true;
		}

	protected:
	//	Postprocess routine
		virtual bool postprocess() {return true;}

protected:
	SmartPtr<ILUTPreconditioner<CPUAlgebra> > ilut;
	SmartPtr<MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type> > mo;
	LinearSolver<CPUAlgebra::vector_type> linearSolver;
	CPUAlgebra::vector_type m_c, m_d;
	double m_eps;
	bool m_info;
	size_t m_size;
};

} // end namespace ug

#endif
