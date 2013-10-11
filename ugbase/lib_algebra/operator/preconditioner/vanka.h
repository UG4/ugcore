/*
 * vanka.h
 *
 *  Created on: 24.07.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__VANKA__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__VANKA__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{

static const int MAXBLOCKSIZE = 53;

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
			x[it.index()] = 0;
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
			x[blockind[j]] = relax*localx[j];
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
	typedef typename Vector_type::value_type vector_block_type;
	DenseVector< VariableArray1<vector_block_type> > s;
	DenseVector< VariableArray1<vector_block_type> > localx;
	DenseMatrix< VariableArray2<number> > mat;
	s.resize(MAXBLOCKSIZE);
	typedef typename Matrix_type::value_type block_type;

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
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		Vanka() {m_relax=1;};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<Vanka<algebra_type> > newInst(new Vanka<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

	///	Destructor
		virtual ~Vanka()
		{};

	/// set relaxation parameter
	public:
		void set_relax(number omega){m_relax=omega;};

	protected:
		number m_relax;

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "Vanka";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
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

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
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
		virtual bool postprocess() {return true;}

	protected:
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
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		DiagVanka() {m_relax=1;};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<DiagVanka<algebra_type> > newInst(new DiagVanka<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

	///	Destructor
		virtual ~DiagVanka()
		{};

		/// set relaxation parameter
	public:
		void set_relax(number omega){m_relax=omega;};

	protected:
		number m_relax;

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "DiagVanka";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
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

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
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
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

};


} // end namespace ug

#endif
