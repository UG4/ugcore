/*
 * \file ilut.h
 *
 * \date 20.07.2013
 * \author Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT_SCALAR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT_SCALAR__

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
			m_eps(eps), m_info(false), m_sort(true)
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

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

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
		
		void set_sort(bool b)
		{
			m_sort = b;
		}

	public:
		void preprocess(const matrix_type &M)
		{
#ifdef 	UG_PARALLEL
			matrix_type A;
			A = M;

			MatAddSlaveRowsToMasterRowOverlap0(A);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex, M.layouts()->slave());
			SetDirichletRow(A, vIndex);
#else
			const matrix_type &A = M;
#endif

			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

			ilut = make_sp(new ILUTPreconditioner<CPUAlgebra>(m_eps));
			ilut->set_info(m_info);
			ilut->set_sort(m_sort);

			mo = make_sp(new MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type>);
			CPUAlgebra::matrix_type &mat = mo->get_matrix();

			#ifdef UG_PARALLEL
				mat.set_storage_type(PST_ADDITIVE);
				mat.set_layouts(CreateLocalAlgebraLayouts());
			#endif

			m_size = GetDoubleSparseFromBlockSparse(mat, A);

			SmartPtr<StdConvCheck<CPUAlgebra::vector_type> > convCheck(new StdConvCheck<CPUAlgebra::vector_type>);
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
#ifdef UG_PARALLEL
			SmartPtr<vector_type> spDtmp = d.clone();
			spDtmp->change_storage_type(PST_UNIQUE);
			bool b = apply_double(c, *spDtmp);

			c.set_storage_type(PST_ADDITIVE);
			c.change_storage_type(PST_CONSISTENT);
			return b;
#else
			return apply_double(c, d);
#endif
			return true;
		}

		bool apply_double(vector_type &c, const vector_type &d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);
			m_c.set_layouts(CreateLocalAlgebraLayouts());
			m_d.set_layouts(CreateLocalAlgebraLayouts());
#endif
			get_vector(m_d, d);
			m_c.set(0.0);

			//ApplyLinearSolver(mo, m_u, m_b, linearSolver);
			ilut->apply(m_c, m_d);
			//ilut->apply(m_c, m_d);

			set_vector(m_c, c);
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
			m_c.set_layouts(CreateLocalAlgebraLayouts());
			m_d.set_layouts(CreateLocalAlgebraLayouts());
#endif
			get_vector(m_d, d);
			m_c.set(0.0);

			//ApplyLinearSolver(mo, m_u, m_b, linearSolver);
			linearSolver.apply_return_defect(m_c, m_d);
			//ilut->apply(m_c, m_d);

			set_vector(m_c, c);
			return true;
		}

	public:
		virtual std::string config_string() const
		{
			std::stringstream ss ; ss << "ILUTScalar(threshold = " << m_eps << ", sort = " << (m_sort?"true":"false") << ")";
			if(m_eps == 0.0) ss << " = Sparse LU";
			return ss.str();
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
	bool m_info, m_sort;
	size_t m_size;
};

} // end namespace ug

#endif
