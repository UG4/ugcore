/*
 * GPUjacobi.h
 *
 *  Created on: 27.05.2013
 *      Author: mrupp
 */

#ifdef CUDA_AVAILABLE
#ifdef UG_GPU
#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GPUJACOBI__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GPUJACOBI__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/gpu_algebra/gpusparsematrix.h"
#include "lib_algebra/gpu_algebra/gpuvector.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

#include "lib_algebra/small_algebra/blocks.h"
#include "lib_algebra/gpu_algebra/gpuvector.h"
#include "lib_algebra/cpu_algebra_types.h"


extern "C" bool
CUDA_JacobiApply(const double *diagInv, double *corr, const double *defect, const int N);

namespace ug{


///	Jacobi Preconditioner
class GPUJacobi : public IPreconditioner<GPUAlgebra>
{
	public:
	///	Algebra type
		typedef GPUAlgebra TAlgebra;
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;
		using base_type::m_spOperator;
		using base_type::damping;

	public:
	///	default constructor
		GPUJacobi() {this->set_damp(1.0);};

	///	constructor setting the damping parameter
		GPUJacobi(number damp) {this->set_damp(damp);};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<GPUJacobi> newInst(new GPUJacobi());
			newInst->set_debug(debug_writer());
			newInst->set_damp(damping());
			return newInst;
		}

	///	Destructor
		virtual ~GPUJacobi()
		{};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "GPUJacobi";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			PROFILE_BEGIN_GROUP(GPUJacobi_preprocess, "algebra GPUJacobi");

			matrix_type &mat = *pOp;
		// 	create help vector to apply diagonal
			size_t size = mat.num_rows();
			if(size != mat.num_cols())
			{
				UG_LOG("Square Matrix needed for GPUJacobi Iteration.\n");
				return false;
			}

		//	get damping in constant case to damp at once
			number damp = 1.0;
			if(damping()->constant_damping())
				damp = damping()->damping();

			if(CheckDiagonalInvertible(mat)==false)
				return false;
//			UG_ASSERT(CheckDiagonalInvertible(mat), "GPUJacobi: A has noninvertible diagonal");


			m_diagInv.resize(mat.num_rows());
		// 	invert diagonal and multiply by damping
			for(size_t i = 0; i < mat.num_rows(); ++i)
				m_diagInv[i] = damp/mat(i,i);

			CUDAHelper::get_instance();
		//	done
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
		// 	the computed correction is additive
			c.set_storage_type(PST_ADDITIVE);

		//	we make it consistent
			if(!c.change_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in 'JacobiPreconditioner::apply': "
						"Cannot change parallel status of correction to consistent.\n");
				return false;
			}
#endif

			const double *devDiagInv = m_diagInv.get_dev_ptr();
			double *devC = c.get_dev_ptr();
			const double *devD = d.get_dev_ptr();
			int N = c.size();
//			UG_LOG("devDiagInv = " << devDiagInv<< " " << m_diagInv.size()  << " devC = " << devC << " " << c.size() << " devD = "
//					<< devD << " " << d.size() << " N = " << N << "\n");
			CUDA_JacobiApply(devDiagInv, devC, devD, N);


			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		GPUVector<double> m_diagInv;

};

template<>
class Jacobi<GPUAlgebra> : public GPUJacobi
{
	using GPUJacobi::set_damp;
public:
	///	default constructor
		Jacobi() {this->set_damp(1.0);};

	///	constructor setting the damping parameter
		Jacobi(number damp) {this->set_damp(damp);};
};


} // end namespace ug

#endif
#endif
#endif // USE_CUDA
