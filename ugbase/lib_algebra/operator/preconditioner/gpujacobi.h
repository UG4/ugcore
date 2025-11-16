/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
		using TAlgebra = GPUAlgebra;
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = TAlgebra::matrix_type;

	///	Matrix Operator type
		using matrix_operator_type = IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;
		using base_type::damping;

	public:
	///	default constructor
		GPUJacobi() {this->set_damp(1.0);};

	///	constructor setting the damping parameter
		GPUJacobi(number damp) {this->set_damp(damp);};

	/// clone constructor
		GPUJacobi( const GPUJacobi &parent )
			: base_type(parent)
		{
		}

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new GPUJacobi(*this));
		}


	///	Destructor
		virtual ~GPUJacobi()
		{};

		virtual bool supports_parallel() const { return false; }
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

		void set_block(bool b)
		{
			// always scalar
		}
};


} // end namespace ug

#endif
#endif
#endif // USE_CUDA
