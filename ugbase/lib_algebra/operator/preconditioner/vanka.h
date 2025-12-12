/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__VANKA__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__VANKA__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"

#ifdef UG_PARALLEL
	//#include "pcl/pcl_util.h"
	//#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap.h"
#endif

namespace ug {

static constexpr size_t MAXBLOCKSIZE = 53;

template<typename Matrix_type, typename Vector_type>
bool Vanka_step(const Matrix_type &A, Vector_type &x, const Vector_type &b, number relax)
{
	DenseVector< VariableArray1<number> > s;
	DenseVector< VariableArray1<number> > localx;
	DenseMatrix< VariableArray2<number> > mat;
	
	size_t blockind[MAXBLOCKSIZE];
	
	for(size_t i=0; i < x.size(); i++)
    {
        x[i]=0;
    };

	for(size_t i=0; i < x.size(); i++)
	{
		if (A(i,i)!=0){
	/*		// do usual gauss-seidel (would be needed in case of Dirichlet pressure condition)
	        typename Vector_type::value_type def = b[i];

            for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) && it.index() < i; ++it)
                // s -= it.value() * x[it.index()];
                MatMultAdd(def, 1.0, def, -1.0, it.value(), x[it.index()]);
            // x[i] = s/A(i,i)
            InverseMatMult(x[i], relax, A(i,i), s);*/
			continue;
		};

		size_t blocksize=0;

		for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
			blockind[blocksize] = it.index();
			blocksize++;
		};
		if (blocksize>MAXBLOCKSIZE) UG_THROW("MAXBLOCKSIZE too small\n");
		mat.resize(blocksize,blocksize);
		s.resize(blocksize);
		localx.resize(blocksize);
		for (size_t j=0;j<blocksize;j++){
			// fill local block matrix
			for (size_t k=0;k<blocksize;k++){
				mat.subassign(j,k,A(blockind[j],blockind[k]));
			};
			// compute rhs
			typename Vector_type::value_type sj = b[blockind[j]];
			for(typename Matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
				MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
			};
			s.subassign(j,sj);
		};
		// solve block
		InverseMatMult(localx,1,mat,s);
		for (size_t j=0;j<blocksize;j++){
			x[blockind[j]] += relax*localx[j];
		};
	}

	return true;
}

// Diagonal Vanka block smoother:
// When setting up the local block matrix the side-diagonal entries are left away, except for the pressure.
// The local block matrix therefore has the form
//
// x 0 ...   0 x
// 0 x 0 ... 0 x
// 0 0 x 0 ... x
//      ...
// 0 ...   0 x x
// x x ... x x x
//
// The velocity off-diagonal entries are considered in the local defect vector.
// The local block system can be solved in O(n) time so that a step of this smoother is computationly cheaper than the
// full Vanka smoother.
template<typename Matrix_type, typename Vector_type>
bool Diag_Vanka_step(const Matrix_type &A, Vector_type &x, const Vector_type &b, number relax)
{
	using vector_block_type = typename Vector_type::value_type;
	DenseVector< VariableArray1<vector_block_type> > s;
	DenseVector< VariableArray1<vector_block_type> > localx;
	DenseMatrix< VariableArray2<number> > mat;
	s.resize(MAXBLOCKSIZE);
	using block_type = typename Matrix_type::value_type;

	size_t blockind[MAXBLOCKSIZE];

	for(size_t i=0; i < x.size(); i++)
    {
        x[i]=0;
    };

	for(size_t i=0; i < x.size(); i++)
	{
		if (A(i,i)!=0){
		/*	// do usual gauss-seidel (would be needed in case of Dirichlet pressure condition)
	        typename Vector_type::value_type def = b[i];

            for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) && it.index() < i; ++it)
                // s -= it.value() * x[it.index()];
                MatMultAdd(def, 1.0, def, -1.0, it.value(), x[it.index()]);
            // x[i] = s/A(i,i)
            InverseMatMult(x[i], relax, A(i,i), def);*/
			continue;
		};

		size_t blocksize=0;

		for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
			if (it.index()==i) continue;
			blockind[blocksize] = it.index();
			s[blocksize] = b[blockind[blocksize]];
			for(typename Matrix_type::const_row_iterator rowit = A.begin_row(blockind[blocksize]); rowit != A.end_row(blockind[blocksize]) ; ++rowit){
					if ((rowit.index()==blockind[blocksize])||(rowit.index()==i)) continue;
					// s[blocksize] -= a_ij*x_j
					MatMultAdd(s[blocksize], 1.0, s[blocksize], -1.0, rowit.value(), x[rowit.index()]);
			};
			blocksize++;
		};
		if (blocksize>MAXBLOCKSIZE) UG_THROW("MAXBLOCKSIZE too small\n");
		// remark: blocksize is without pressure variable, so actual blocksize is blocksize+1
		block_type a_ii = A(i,i);
		typename Vector_type::value_type s_i = b[i];
		// Gauss elimination on local block matrix
		for (size_t j=0;j<blocksize;j++){
			block_type a_q = A(i,blockind[j]);
			block_type a_jj =  A(blockind[j],blockind[j]);
			a_q /= a_jj;
			// s_i -= a_ij/a_jj*s_j
			MatMultAdd(s_i, 1.0, s_i, -1.0, a_q, s[j]);
			// a_ii -= a_ij/a_jj*a_ji
			a_ii-=a_q*A(blockind[j],i);
		}
		// solve diagonalized system
		// x[i] = s_i/a_ii
		InverseMatMult(x[i], 1.0, a_ii, s_i);
		for (size_t j=0;j<blocksize;j++){
			 // s_j-=a_ji*x_i
			 MatMultAdd(s[j], 1.0, s[j], -1.0, A(blockind[j],i), x[i]);
			 // x_j=1/a_jj*s_j
			 InverseMatMult(x[blockind[j]], relax, A(blockind[j],blockind[j]),s[j]);
		}
	}
	return true;
}


///	Vanka Preconditioner
template <typename TAlgebra>
class Vanka : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		Vanka() {m_relax=1;};

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			SmartPtr<Vanka > newInst(new Vanka());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

	///	Destructor
		~Vanka() override = default;

		bool supports_parallel() const override {return true;}

	/// set relaxation parameter
	public:
		void set_relax(number omega){m_relax=omega;};


	protected:
	///	Name of preconditioner
		const char* name() const override {return "Vanka";}

	///	Preprocess routine
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override {
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override {
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);

				if(!Vanka_step(m_A, c, dhelp, m_relax)) return false;

				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				if(!Vanka_step(*pOp, c, d, m_relax)) return false;

#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	///	Postprocess routine
		bool postprocess() override {return true;}

protected:
	number m_relax;
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

};

///	Diagvanka Preconditioner, description see above diagvanka_step function
template <typename TAlgebra>
class DiagVanka : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		DiagVanka() {m_relax=1;};

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			SmartPtr<DiagVanka > newInst(new DiagVanka());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

		[[nodiscard]] bool supports_parallel() const override {return true;}

	///	Destructor
		~DiagVanka() override = default;

		/// set relaxation parameter
	public:
		void set_relax(number omega){m_relax=omega;};

	protected:
	///	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "DiagVanka";}

	///	Preprocess routine
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override {
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override {
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);

				if(!Diag_Vanka_step(m_A, c, dhelp, m_relax)) return false;

				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{

				if(!Diag_Vanka_step(*pOp, c, d, m_relax)) return false;

#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	///	Postprocess routine
		bool postprocess() override {return true;}

protected:
	number m_relax;
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

};


} // end namespace ug

#endif
