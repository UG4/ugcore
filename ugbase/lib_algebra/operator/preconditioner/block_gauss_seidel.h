/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/algebra_common/core_smoothers.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

#include "lib_algebra/algebra_common/sparsematrix_util.h"
#include "lib_algebra/algebra_common/local_helper.h"
#include "lib_algebra/small_algebra/small_matrix/print.h"
#include "lib_algebra/operator/preconditioner/ilut.h"
#include "common/debug_print.h"
#include "common/progress.h"
#include "lib_algebra/small_algebra/additional_math.h"

#include "lib_algebra/common/graph/graph.h"
#include "lib_algebra/algebra_common/sparse_vector.h"



namespace ug{

//////////////////// FROM AMG //////////////////////////////

template<typename TMatrixType, typename TRowType>
void CopyOffDiagEntries(const TMatrixType &A, size_t i, TRowType &row, bool enforceNew=false)
{
	// copy matrix entries selected by rule
	for(typename TMatrixType::const_row_iterator connij = A.begin_row(i); connij != A.end_row(i); ++connij)
	{
		const size_t j=connij.index();
		row(connij.index()) = (i==j) ? 0.0 : connij.value();
	}
}

//! Adds 'strong negative connections' to graph.
/*! Criterion: \f$ -a_{ij} \ge \epsilon \max_{a_{ik}<0} |a_{ik}| \f$ */
class StrongNegativeConnectionsByBlockNorm
{
protected:
	double m_theta;
public:
	explicit StrongNegativeConnectionsByBlockNorm(double theta) : m_theta(theta){}

	template<typename TRowType>
	void find(const TRowType &Ci, size_t i, cgraph &graph)
	{
		double dmax = 0.0;
		using const_iterator = typename TRowType::const_iterator;

		for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			double d = BlockNorm(conn.value());
//			std::cout << d << std::endl;
			if(d > dmax) dmax = BlockNorm(conn.value());
		}

		for(const_iterator conn = Ci.begin(); conn != Ci.end(); ++conn)
		{
			if(conn.index() == i) continue; // skip diag
			if( BlockNorm(conn.value()) >= m_theta * dmax)
				graph.set_connection(i, conn.index());
		}
	}
};


template<typename matrix_type>
void CreateStrongConnectionGraphForSystems(const matrix_type &A, cgraph &graph, double theta)
{
	PROFILE_FUNC();
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		// Skip isolated node
		if(A.is_isolated(i)) continue;

		// copy all off-diag entries
		SparseVector<typename matrix_type::value_type> Ri(A.num_cols());
		CopyOffDiagEntries(A, i, Ri);

		// Find strong connections
		StrongNegativeConnectionsByBlockNorm strong(theta);
		strong.find(Ri, i, graph);
	}
}


/////////////////////////////////

/**
 * Calculates the gs correction by updating all unknowns in 'indices' at once
 * with the inverse of this block stored in AlocInv.
 * @param A					a SparseMatrix
 * @param x					solution to be GS-updated
 * @param b					right hand side b
 * @param AlocInv			cached inverse on a block/slice, calculated with GetSliceDenseInverse
 * @param indices			the indices of the block
 * @param tmp				a temporary vector
 */
// we cache the inverse! computing the inverse "on the fly" seems to be a lot slower.
template<typename TSparseMatrixType, typename TVectorType>
void GetBlockGSCorrection(const TSparseMatrixType &A, TVectorType &x, TVectorType &b,
		DenseMatrix<VariableArray2<double> > &AlocInv,
		std::vector<size_t> &indices, DenseVector<VariableArray1<double> > &tmp, DenseVector<VariableArray1<double> > &tmp2)
{
	// assume tmp is big enough
	size_t k;

	using const_row_iterator = typename TSparseMatrixType::const_row_iterator;
	using smallvec_type = typename TVectorType::value_type;

	k = 0;
	tmp.resize(indices.size()*GetSize(b[0]));
	tmp2.resize(indices.size()*GetSize(b[0]));
	for(size_t i=0; i<indices.size(); i++)
	{
		int j = indices[i];
		smallvec_type s = b[j];
		// calc b-Ax
		for(const_row_iterator it = A.begin_row(j); it != A.end_row(j); ++it)
			s -= it.value() * x[it.index()];
		for(size_t u=0; u<GetSize(s); u++)
			tmp[k++] = BlockRef(s, u);
	}
	MatMult(tmp2, 1.0, AlocInv, tmp);

	k=0;
	for(size_t i=0; i<indices.size(); i++)
	{
		smallvec_type &xi = x[indices[i]];
		for(size_t j=0; j<GetSize(xi); j++)
			BlockRef(xi, j) += tmp2[k++];
	}
}
template<typename TSparseMatrixType, typename TVectorType>
void GetBlockGSCorrection(const TSparseMatrixType &A, TVectorType &x, TVectorType &b,
		DenseMatrix<VariableArray2<double> > &AlocInv,
		std::vector<size_t> &indices)
{
	DenseVector<VariableArray1<double> > tmp;
	DenseVector<VariableArray1<double> > tmp2;
	GetBlockGSCorrection(A, x, b, AlocInv, indices, tmp, tmp2);
}

/**
 * Calculates the gs correction by updating all unknowns in 'indices' at once
 * with the inverse of this block implicitely in 'ilut'.
 * @param A					a SparseMatrix
 * @param x					solution to be GS-updated
 * @param b					right hand side b
 * @param ilut				cached ilut inverse type
 * @param indices			the indices of the block
 * @param tmp				a temporary vector
 * @param tmp2				another temporary vector
 */
template<typename TSparseMatrixType, typename TVectorType>
void GetBlockGSCorrectionILUT(const TSparseMatrixType &A, TVectorType &x, TVectorType &b,
		SmartPtr<ILUTPreconditioner<CPUAlgebra> > &ilut,
		std::vector<size_t> &indices, CPUAlgebra::vector_type &tmp, CPUAlgebra::vector_type &tmp2)
{
	size_t blockSize = GetSize(b[0]);
	using const_row_iterator = typename TSparseMatrixType::const_row_iterator;
	using smallvec_type = typename TVectorType::value_type;

	size_t k = 0;
	tmp.resize(indices.size()*blockSize);
	tmp2.resize(indices.size()*blockSize);
	for(size_t i=0; i<indices.size(); i++)
	{
		int j = indices[i];
		smallvec_type s = b[j];
		// calc b-Ax
		for(const_row_iterator it = A.begin_row(j); it != A.end_row(j); ++it)
			s -= it.value() * x[it.index()];
		for(size_t u=0; u<GetSize(s); u++)
			tmp[k++] = BlockRef(s, u);
	}
	ilut->solve(tmp2, tmp);

	k=0;
	for(size_t i=0; i<indices.size(); i++)
	{
		smallvec_type &xi = x[indices[i]];
		for(size_t j=0; j<GetSize(xi); j++)
			BlockRef(xi, j) += tmp2[k++];
	}
}

/**
 * @param A			a sparse matrix
 * @param indices	map local -> global indices
 * @param AlocInv	inverse on the indices to be used in @sa GetBlockGSCorrection
 */
template<typename TSparseMatrixType>
void GetSliceDenseInverse(const TSparseMatrixType &A, const std::vector<size_t> &indices,
		DenseMatrix<VariableArray2<double> > &AlocInv,
		DenseMatrix<VariableArray2<typename TSparseMatrixType::value_type> > &tmp)
{
	tmp.resize(indices.size(), indices.size());
	GetLocalMatrix(A, tmp, &indices[0], &indices[0]);
	BlockMatrixToDoubleMatrix(AlocInv, tmp);
	Invert(AlocInv);
}

template<typename TSparseMatrixType>
void GetSliceDenseInverse(const TSparseMatrixType &A, const std::vector<size_t> &indices,
		DenseMatrix<VariableArray2<double> > &AlocInv)
{
	DenseMatrix<VariableArray2<typename TSparseMatrixType::value_type> > L;
	GetSliceDenseInverse(A, indices, AlocInv, L);
}

/**
 * @param A			a sparse matrix
 * @param indices	map local -> global indices
 * @param R			sparse matrix slice in local indices
 */
template<typename TSparseMatrixType>
void GetSliceSparse(const TSparseMatrixType &A, const std::vector<size_t> &indices,
		CPUAlgebra::matrix_type &R)
{
	size_t blockSize = GetRows(A(0,0));
	size_t numRows = indices.size();
	size_t numCols = indices.size();

	R.resize_and_clear(numRows*blockSize, numCols*blockSize);
	for(size_t ri=0; ri<indices.size(); ri++)
	{
		size_t r=indices[ri];
		for(size_t ci=0; ci<indices.size(); ci++)
		{
			size_t c = indices[ci];
			if(A.has_connection(r, c))
			{
				const typename TSparseMatrixType::value_type &m = A(r, c);
				for(size_t j1=0; j1<blockSize; j1++)
					for(size_t j2=0; j2<blockSize; j2++)
						R(ri*blockSize + j1, ci*blockSize + j2) = BlockRef(m, j1, j2);
			}
		}
	}
	R.defragment();
}


template<typename TAlgebra>
class IBlockJacobiPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		using algebra_type = TAlgebra;

	//	Vector type
		using vector_type = typename TAlgebra::vector_type;

	//	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

		IBlockJacobiPreconditioner() = default;
		IBlockJacobiPreconditioner(const IBlockJacobiPreconditioner &parent) : IPreconditioner<TAlgebra>(parent)
		{ }

		~IBlockJacobiPreconditioner() override = default;


		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override {
			matrix_type &mat = *pOp;

#ifdef 	UG_PARALLEL
			A = mat;
			MatAddSlaveRowsToMasterRowOverlap0(A);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex,  mat.layouts()->slave());
			SetDirichletRow(A, vIndex);
			return block_preprocess(A);
#else
			return block_preprocess(mat);
#endif
		}

		//	Preprocess routine
		virtual bool block_preprocess(matrix_type &A) = 0;

	//	Postprocess routine
	bool postprocess() override {return true;}
	[[nodiscard]] bool supports_parallel() const override { return true; }

	//	Stepping routine
	bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override {
#ifdef UG_PARALLEL
			if(!m_spDtmp.valid() || m_spDtmp->size() != d.size())
				m_spDtmp = d.clone();
			else
				(*m_spDtmp) = d;
			m_spDtmp->change_storage_type(PST_UNIQUE);
			bool b = block_step(A, c, *m_spDtmp);
			c.set_storage_type(PST_ADDITIVE);
			c.change_storage_type(PST_CONSISTENT);
			return b;
#else
			return block_step(*pOp, c, d);
#endif
			return true;
		}

	virtual bool block_step(matrix_type &A, vector_type& c, const vector_type& d) = 0;

		//	Stepping routine
protected:

#ifdef 	UG_PARALLEL
	matrix_type A;
#endif

		SmartPtr<vector_type> m_spDtmp;

};

/**
 * BlockGaussSeidel
 * use the constructor or set_depth to set the depth
 * depth = 0 -> GS
 * depth = 1 -> Block i and neighbors of i
 * depth = 2 -> Block i and neighbors of i and their neighbors
 * depth = 3 ...
 */
template <typename TAlgebra, bool backward, bool forward>
class BlockGaussSeidel : public IBlockJacobiPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
	using algebra_type = TAlgebra;

	//	Vector type
	using vector_type = typename TAlgebra::vector_type;
	using smallvec_type = typename vector_type::value_type;
	//	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;
	using smallmat_type = typename matrix_type::value_type;
	///	Matrix Operator type
	using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
	using base_type = IBlockJacobiPreconditioner<TAlgebra>;

	using this_type = BlockGaussSeidel;

	using const_row_iterator = typename matrix_type::const_row_iterator;

protected:
	using block_type = typename matrix_type::value_type;

public:
	//	Constructor
		BlockGaussSeidel() {
			m_depth = 1;
		};

		explicit BlockGaussSeidel(int depth) {
			m_depth = depth;
		};

		BlockGaussSeidel(const this_type &parent) : base_type(parent)
		{
			m_depth = parent.m_depth;
		}
		~BlockGaussSeidel() override = default;

	// 	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new this_type(*this));
		}

		void set_depth(size_t d)
		{
			m_depth = d;
		}

	protected:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "BlockGaussSeidel";}

		//	Preprocess routine
		bool block_preprocess(matrix_type &A) override {
			size_t N = A.num_rows();
			DenseMatrix<VariableArray2<smallmat_type> > tmpMat;

			AlocInv.clear();
			AlocInv.resize(N);

			indices.resize(N);

			size_t maxSize = 0;
			PROGRESS_START(prog, N, "BlockGaussSeidel: compute blocks");
			for(size_t i=0; i<N; i++)
			{
				PROGRESS_UPDATE(prog, i);
				std::vector<bool> bVisited(N, false);
				indices[i].clear();
				GetNeighborhood(A, i, m_depth, indices[i], bVisited);
				//indices[i].push_back(i);

				GetSliceDenseInverse(A, indices[i], AlocInv[i], tmpMat);

				//PrintVector(indices[i], "indices");
//				UG_LOG("AlocInv " << i << ":\n" << JuliaString(AlocInv[i], "Aloc") << "\n");
//				UG_LOG("AlocInv " << i << ":\n" << JuliaString(AlocInv[i], "AlocInv") << "\n");

				if(AlocInv[i].num_rows() > maxSize) maxSize = AlocInv[i].num_rows();
			}
			PROGRESS_FINISH(prog);

			UG_LOG("Max Size = " << maxSize << "\n");
			return true;
		}

		using matrix_const_row_iterator = typename matrix_type::const_row_iterator;
		using matrix_row_iterator = typename matrix_type::row_iterator;

	//	Postprocess routine
		bool postprocess() override {return true;}
		[[nodiscard]] bool supports_parallel() const override { return true; }

		//	Stepping routine
		bool block_step(matrix_type &A, vector_type& c, const vector_type& d) override {
			PROFILE_BEGIN(BlockGaussSeidel_step);
			vector_type &x = c;
			x.set(0.0);
			vector_type b;
			b = d;

			DenseVector<VariableArray1<double> > tmp;
			DenseVector<VariableArray1<double> > tmp2;

			if(forward)
				for(size_t i=0; i<x.size(); i++)
				{
					// c = D^{-1}(b-Ax)
					// x = x + c
					GetBlockGSCorrection(A, x, b, AlocInv[i], indices[i], tmp, tmp2);
	//				do_correction_implicit(A, x, b, indices[i]);
				}
			if(backward)
				for(size_t i=x.size()-1; ; i--)
				{
					// c = D^{-1}(b-Ax)
					// x = x + c
					GetBlockGSCorrection(A, x, b, AlocInv[i], indices[i], tmp, tmp2);
					if(i==0) break;
	//				do_correction_implicit(A, x, b, indices[i]);
				}

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			c.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

		[[nodiscard]] std::string config_string() const override {
			std::stringstream ss ; ss << "BlockGaussSeidel(depth = " << m_depth << ")";
			return ss.str();
		}

public:
	std::vector< DenseMatrix<VariableArray2<double> > > AlocInv;
	size_t m_depth;
	std::vector<std::vector<size_t> > indices;

};


template <typename TAlgebra, bool backward, bool forward>
class BlockGaussSeidelIterative : public IBlockJacobiPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		using algebra_type = TAlgebra;

	//	Vector type
		using vector_type = typename TAlgebra::vector_type;
		using smallvec_type = typename vector_type::value_type;
	//	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;
		using smallmat_type = typename matrix_type::value_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IBlockJacobiPreconditioner<TAlgebra>;

		using this_type = BlockGaussSeidelIterative;

		using const_row_iterator = typename matrix_type::const_row_iterator;

	protected:
		using block_type = typename matrix_type::value_type;
		using matrix_const_row_iterator = typename matrix_type::const_row_iterator;
		using matrix_row_iterator = typename matrix_type::row_iterator;

	public:
	//	Constructor
		BlockGaussSeidelIterative() {
			m_depth = 1;
			m_nu = 2;
		};

		BlockGaussSeidelIterative(int depth, int nu) {
			m_depth = depth;
			m_nu = nu;
		};

		BlockGaussSeidelIterative(const this_type &parent) : base_type(parent)
		{
			m_depth = parent.m_depth;
			m_nu = parent.m_nu;
		}

		~BlockGaussSeidelIterative() override = default;

	// 	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new this_type(*this));
		}


		void set_depth(size_t d)
		{
			m_depth = d;
		}


	protected:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {
			if(backward&&forward) return "SymmetricBlockGaussSeidelIterative";
			else if(backward) return "BackwardBlockGaussSeidelIterative";
			return "BlockGaussSeidelIterative";
		}



		//	Preprocess routine
		bool block_preprocess(matrix_type &A) override {
			size_t N = A.num_rows();
			indices.resize(N);

			size_t maxSize = 0;
			for(size_t i=0; i<N; i++)
			{
				std::vector<bool> bVisited(N, false);
				indices[i].clear();
				GetNeighborhood(A, i, m_depth, indices[i], bVisited);
				maxSize = std::max(indices[i].size(), maxSize);
			}
			UG_LOG("Max Size = " << maxSize << "\n");
			return true;
		}



	//	Postprocess routine
		bool postprocess() override {return true;}
		[[nodiscard]] bool supports_parallel() const override { return true; }


		void correct(size_t i, const matrix_type &A, vector_type& x, const vector_type& b)
		{

			smallvec_type s = b[i];
			// calc b-Ax
			for(const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
				s -= it.value() * x[it.index()];
			smallvec_type c;
			InverseMatMult(c, 1.0, A(i,i), s);
			x[i] += c;

		}

		void correct_forward(size_t i, matrix_type &A, vector_type& x, const vector_type& b)
		{
			for(size_t k=0; k<m_nu; k++) {
				for(size_t j=0; j<indices[i].size(); j++) {
					correct(indices[i][j], A, x, b);
				}
			}
		}

		void correct_backward(size_t i, matrix_type &A, vector_type& x, const vector_type& b)
		{
			for(size_t k=0; k<m_nu; k++) {
				for(int j=(int)(indices[i].size())-1; j>=0 ; j--) {
					correct(indices[i][j], A, x, b);
				}
			}
		}



		//	Stepping routine
		bool block_step(matrix_type &A, vector_type& c, const vector_type& d) override {
			PROFILE_BEGIN(BlockGaussSeidelIterative_step);
			vector_type &x = c;
			x.set(0.0);
			vector_type b;
			b = d;

			if(forward)
				for(size_t i=0; i<x.size(); ++i)
				{
					correct_forward(i, A, x, b);
				}
			if(backward)
				for(size_t i=x.size()-1; ; --i)
				{
					correct_backward(i, A, x, b);
					if(i==0) break;
				}

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			c.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

	public:
			void set_iterative_steps(size_t nu)
			{
				m_nu=nu;
			}


		[[nodiscard]] std::string config_string() const override {
			std::stringstream ss ;
			if(backward&&forward) ss << "Symmetric";
			else if(backward) ss << "Backward";
			ss << "BlockGaussSeidelIterative(depth = " << m_depth << ", nu = " << m_nu << ")";
			return ss.str();
		}
public:
		std::vector< DenseMatrix<VariableArray2<double> > > AlocInv;
		size_t m_depth;
		std::vector<std::vector<size_t> > indices;
		size_t m_nu;
};


/**
 * SparseBlockGaussSeidel
 * experimental version
 * a) can use bigger stencils since it uses SparseLU for solving blocks
 * b) tries to use some overlapping blocks (BlockGaussSeidel always uses N blocks)
 */
template <typename TAlgebra, bool backward, bool forward>
class SparseBlockGaussSeidel : public IBlockJacobiPreconditioner<TAlgebra>
{
public:
	//	Algebra type
	using algebra_type = TAlgebra;

	//	Vector type
	using vector_type = typename TAlgebra::vector_type;
	using smallvec_type = typename vector_type::value_type;
	//	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;
	using smallmat_type = typename matrix_type::value_type;
	///	Matrix Operator type
	using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
	using base_type = IBlockJacobiPreconditioner<TAlgebra>;

	using const_row_iterator = typename matrix_type::const_row_iterator;
protected:
	using block_type = typename matrix_type::value_type;
	using this_type = SparseBlockGaussSeidel;

	public:
	//	Constructor
		SparseBlockGaussSeidel() {
			m_depth = 1;
		};

		explicit SparseBlockGaussSeidel(int depth) {
			m_depth = depth;
		};

		SparseBlockGaussSeidel(const this_type &parent) : base_type(parent)
		{
			m_depth = parent.m_depth;
		}

		~SparseBlockGaussSeidel() override = default;

	// 	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new this_type(*this));
		}

		void set_depth(size_t d)
		{
			m_depth = d;
		}


	protected:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "SparseBlockGaussSeidel";}

		//	Preprocess routine
		bool block_preprocess(matrix_type &A) override {
			size_t N = A.num_rows();
			DenseMatrix<VariableArray2<smallmat_type> > tmpMat;


			m_ilut.clear();

			indices.resize(N);

			size_t maxSize = 0;
			PROGRESS_START(prog, N, "SparseBlockGaussSeidel: compute blocks");



			std::vector<size_t> bVisited(N, 0);
			std::vector<bool> bVisited2(N, false);

			std::vector<size_t> levels;

			for(size_t i=0; i<N; i++)
			{
//				if(bVisited[i]>5) continue;
				PROGRESS_UPDATE(prog, i);

				indices[i].clear();

				GetNeighborhood(A, i, m_depth, indices[i], bVisited2);
				for(size_t j=0; j<indices[i].size(); j++)
					bVisited[indices[i][j]] ++;
				//for(size_t j=0; j<levels[m_depth-1]; j++)
					//bVisited2[indices[j]]=true;

				Aloc[i] = make_sp(new CPUAlgebra::matrix_type);

				GetSliceSparse(A, indices[i], *Aloc[i]);


				m_ilut[i] = make_sp(new ILUTPreconditioner<CPUAlgebra> (0.0));
				m_ilut[i]->preprocess_mat(*Aloc[i]);



				if(Aloc[i]->num_rows() > maxSize) maxSize = Aloc[i]->num_rows();
			}
			PROGRESS_FINISH(prog);

			UG_LOG("Max Size = " << maxSize << "\n");
			return true;
		}

		using matrix_const_row_iterator = typename matrix_type::const_row_iterator;
		using matrix_row_iterator = typename matrix_type::row_iterator;

	//	Postprocess routine
		bool postprocess() override {return true;}
		[[nodiscard]] bool supports_parallel() const override { return true; }

		//	Stepping routine
		bool block_step(matrix_type &A, vector_type& c, const vector_type& d) override {
			vector_type &x = c;
			x.set(0.0);
			vector_type b;
			b = d;

			DenseMatrix<VariableArray2<double> > Adense;
			DenseMatrix<VariableArray2<smallmat_type> > Atmp;
			CPUAlgebra::vector_type tmp, tmp2;
			PROGRESS_START(prog, x.size(), "SparseBlockGaussSeidel: step");
			if(forward)
				for(size_t i=0; i<x.size(); ++i)
				{
					PROGRESS_UPDATE(prog, i);
					// c = D^{-1}(b-Ax)
					// x = x + c
					if(indices[i].size() != 0)
						GetBlockGSCorrectionILUT(A, x, b, m_ilut[i], indices[i], tmp, tmp2);
				}
			if(backward)
				for(size_t i=x.size()-1; ; --i)
				{
					PROGRESS_UPDATE(prog, i);
					// c = D^{-1}(b-Ax)
					// x = x + c
					if(indices[i].size() != 0)
						GetBlockGSCorrectionILUT(A, x, b, m_ilut[i], indices[i], tmp, tmp2);
					if(i==0) break;
				}
			PROGRESS_FINISH(prog);

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			c.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

		[[nodiscard]] std::string config_string() const override {
			std::stringstream ss ;
			if(backward&&forward) ss << "Symmetric";
			else if(backward) ss << "Backward";
			ss << "SparseBlockGaussSeidel(depth = " << m_depth << ")";
			return ss.str();
		}

		std::map<size_t, SmartPtr<ILUTPreconditioner<CPUAlgebra> > > m_ilut;

		size_t m_depth;
		std::vector<std::vector<size_t> > indices;
		std::map<size_t, SmartPtr<CPUAlgebra::matrix_type> > Aloc;

};


/**
 * SparseBlockGaussSeidel
 * experimental version
 * a) can use bigger stencils since it uses SparseLU for solving blocks
 * b) tries to use some overlapping blocks (BlockGaussSeidel always uses N blocks)
 */
template <typename TAlgebra, bool backward, bool forward>
class SparseBlockGaussSeidel2 : public IBlockJacobiPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		using algebra_type = TAlgebra;

	//	Vector type
		using vector_type = typename TAlgebra::vector_type;
		using smallvec_type = typename vector_type::value_type;
	//	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;
		using smallmat_type = typename matrix_type::value_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IBlockJacobiPreconditioner<TAlgebra>;
		using this_type = SparseBlockGaussSeidel2;

		using const_row_iterator = typename matrix_type::const_row_iterator;

	protected:
		using block_type = typename matrix_type::value_type;
		using matrix_const_row_iterator = typename matrix_type::const_row_iterator;
		using matrix_row_iterator = typename matrix_type::row_iterator;

	public:
	//	Constructor
		SparseBlockGaussSeidel2() {
			m_depth = 1;
		};

		explicit SparseBlockGaussSeidel2(int depth) {
			m_depth = depth;
		};

		SparseBlockGaussSeidel2(const this_type &parent) : base_type(parent)
		{
			m_depth = parent.m_depth;
		}

		~SparseBlockGaussSeidel2() override = default;

	// 	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new this_type(*this));
		}



		void set_depth(size_t d)
		{
			m_depth = d;
		}


	protected:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "SparseBlockGaussSeidel2";}

		//	Preprocess routine
		bool block_preprocess(matrix_type &A) override {
			size_t N = A.num_rows();
			DenseMatrix<VariableArray2<smallmat_type> > tmpMat;
			m_ilut.clear();

			indices.resize(N);

			size_t maxSize = 0;
			PROGRESS_START(prog, N, "SparseBlockGaussSeidel: compute blocks");

			cgraph G;
			{
				cgraph graph;
				CreateStrongConnectionGraphForSystems(A, graph, 0.3);
				G.resize(N);
				for(size_t i=0; i<graph.size(); i++)
					for(auto it = graph.begin_row(i); it != graph.end_row(i); ++it)
					{
						G.set_connection(*it, i);
						G.set_connection(i, *it);
					}
			}
			std::vector<int> iComponent(N, -1);
			std::vector<std::vector<int> > components;
//			G.print();
			for(size_t i=0; i<G.size(); i++)
			{
				if(iComponent[i] != -1) continue;
				int myComponent = iComponent[i] = components.size();
				components.resize(components.size()+1);
				std::vector<int> &myComponents = components[myComponent];
				myComponents.push_back(i);
				for(size_t c=0; c<myComponents.size(); c++)
				{
					int j = myComponents[c];
					for(auto it = G.begin_row(j); it != G.end_row(j); ++it)
					{
						if(iComponent[*it] == myComponent) continue;
						myComponents.push_back(*it);
						UG_COND_THROW(iComponent[*it] != -1, i );
						iComponent[*it] = myComponent;
					}
				}
			}
			for(size_t i=0; i<N; i++)
				indices[i].clear();
			for(size_t c=0; c<components.size(); c++)
			{
				int i=components[c][0];
				indices[i].insert(indices[i].begin(), components[c].begin(), components[c].end());
				PRINT_VECTOR(indices[i], i);

				Aloc[i] = make_sp(new CPUAlgebra::matrix_type);
				GetSliceSparse(A, indices[i], *Aloc[i]);

				m_ilut[i] = make_sp(new ILUTPreconditioner<CPUAlgebra> (0.0));
				m_ilut[i]->preprocess_mat(*Aloc[i]);
				if(Aloc[i]->num_rows() > maxSize) maxSize = Aloc[i]->num_rows();
			}
			PROGRESS_FINISH(prog);

			UG_LOG("Max Size = " << maxSize << "\n");
			return true;
		}


	//	Postprocess routine
		bool postprocess() override {return true;}
		[[nodiscard]] bool supports_parallel() const override { return true; }

		//	Stepping routine
		bool block_step(matrix_type &A, vector_type& c, const vector_type& d) override {
			vector_type &x = c;
			x.set(0.0);
			vector_type b;
			b = d;

			DenseMatrix<VariableArray2<double> > Adense;
			DenseMatrix<VariableArray2<smallmat_type> > Atmp;
			CPUAlgebra::vector_type tmp, tmp2;
			PROGRESS_START(prog, x.size(), "SparseBlockGaussSeidel: step");
			if(forward)
				for(size_t i=0; i<x.size(); ++i)
				{
					PROGRESS_UPDATE(prog, i);
					// c = D^{-1}(b-Ax)
					// x = x + c
					if(indices[i].size() != 0)
						GetBlockGSCorrectionILUT(A, x, b, m_ilut[i], indices[i], tmp, tmp2);
				}
			if constexpr (backward)
				for(size_t i=x.size()-1; ; --i)
				{
					PROGRESS_UPDATE(prog, i);
					// c = D^{-1}(b-Ax)
					// x = x + c
					if(indices[i].size() != 0)
						GetBlockGSCorrectionILUT(A, x, b, m_ilut[i], indices[i], tmp, tmp2);
					if(i==0) break;
				}
			PROGRESS_FINISH(prog);

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			c.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

		[[nodiscard]] std::string config_string() const override {
			std::stringstream ss ;
			if(backward&&forward) ss << "Symmetric";
			else if(backward) ss << "Backward";
			ss << "SparseBlockGaussSeidel(depth = " << m_depth << ")";
			return ss.str();
		}
public:

	std::map<size_t, SmartPtr<ILUTPreconditioner<CPUAlgebra> > > m_ilut;

	size_t m_depth;
	std::vector<std::vector<size_t> > indices;
	std::map<size_t, SmartPtr<CPUAlgebra::matrix_type> > Aloc;

};


} // end namespace ug

#endif