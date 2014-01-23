/*
 * gauss.h
 *
 *  Created on: 14.07.2010
 *      Author: Martin Rupp
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


namespace ug{


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
		std::vector<size_t> &indices, DenseVector<VariableArray1<double> > &tmp)
{
	// assume tmp is big enough
	size_t k;

	typedef typename TSparseMatrixType::const_row_iterator const_row_iterator;
	typedef typename TVectorType::value_type smallvec_type;

	k = 0;
	tmp.resize(indices.size()*GetSize(b[0]));
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
	tmp = AlocInv*tmp;

	k=0;
	for(size_t i=0; i<indices.size(); i++)
	{
		smallvec_type &xi = x[indices[i]];
		for(size_t j=0; j<GetSize(xi); j++)
			BlockRef(xi, j) += tmp[k++];
	}
}
template<typename TSparseMatrixType, typename TVectorType>
void GetBlockGSCorrection(const TSparseMatrixType &A, TVectorType &x, TVectorType &b,
		DenseMatrix<VariableArray2<double> > &AlocInv,
		std::vector<size_t> &indices)
{
	DenseVector<VariableArray1<double> > tmp;
	GetBlockGSCorrection(A, x, b, AlocInv, indices);
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
	size_t k;
	typedef typename TSparseMatrixType::const_row_iterator const_row_iterator;
	typedef typename TVectorType::value_type smallvec_type;

	k = 0;
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
	typedef typename TSparseMatrixType::const_row_iterator const_row_iterator;
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



/**
 * BlockGaussSeidel
 * use the constructor or set_depth to set the depth
 * depth = 0 -> GS
 * depth = 1 -> Block i and neighbors of i
 * depth = 2 -> Block i and neighbors of i and their neighbors
 * depth = 3 ...
 */
template <typename TAlgebra>
class BlockGaussSeidel : public IPreconditioner<TAlgebra>
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
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		typedef typename matrix_type::value_type block_type;

		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	//	Constructor
		BlockGaussSeidel() {
			m_depth = 1;
		};

		BlockGaussSeidel(int depth) {
			m_depth = depth;
		};

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<BlockGaussSeidel<algebra_type> > newInst(new BlockGaussSeidel<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_depth(m_depth);
			return newInst;
		}

		typedef typename matrix_type::value_type smallmat_type;
		typedef typename vector_type::value_type smallvec_type;
		typedef typename matrix_type::const_row_iterator const_row_iterator;

		std::vector< DenseMatrix<VariableArray2<double> > > AlocInv;
		size_t m_depth;
		std::vector<std::vector<size_t> > indices;


		void set_depth(size_t d)
		{
			m_depth = d;
		}


	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "BlockGaussSeidel";}


		//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type &A = *pOp;

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

		typedef typename matrix_type::const_row_iterator matrix_const_row_iterator;
		typedef typename matrix_type::row_iterator matrix_row_iterator;

	//	Postprocess routine
		virtual bool postprocess() {return true;}
		virtual bool supports_parallel() const { return true; }

		//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			matrix_type &A = *pOp;
			vector_type &x = c;
			x.set(0.0);
			vector_type b;
			b = d;

			DenseVector<VariableArray1<double> > tmp;
			PROGRESS_START(prog, x.size(), "BlockGaussSeidel: step");
			for(size_t i=0; i<x.size(); i++)
			{
				PROGRESS_UPDATE(prog, i);
				// c = D^{-1}(b-Ax)
				// x = x + c
				GetBlockGSCorrection(A, x, b, AlocInv[i], indices[i], tmp);
//				do_correction_implicit(A, x, b, indices[i]);
			}
			PROGRESS_FINISH(prog);

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			c.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss ; ss << "BlockGaussSeidel(depth = " << m_depth << ")";
			return ss.str();
		}

};


/**
 * BlockGaussSeidel2
 * experimental version
 * a) can use bigger stencils since it uses SparseLU for solving blocks
 * b) tries to use some overlapping blocks (BlockGaussSeidel always uses N blocks)
 */
template <typename TAlgebra>
class BlockGaussSeidel2 : public IPreconditioner<TAlgebra>
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
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		typedef typename matrix_type::value_type block_type;

		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	//	Constructor
		BlockGaussSeidel2() {
			m_depth = 1;
		};

		BlockGaussSeidel2(int depth) {
			m_depth = depth;
		};

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<BlockGaussSeidel2<algebra_type> > newInst(new BlockGaussSeidel2<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_depth(m_depth);
			return newInst;
		}

		typedef typename matrix_type::value_type smallmat_type;
		typedef typename vector_type::value_type smallvec_type;
		typedef typename matrix_type::const_row_iterator const_row_iterator;

		std::map<size_t, SmartPtr<ILUTPreconditioner<CPUAlgebra> > > m_ilut;

		size_t m_depth;
		std::vector<std::vector<size_t> > indices;
		std::map<size_t, SmartPtr<CPUAlgebra::matrix_type> > Aloc;


		void set_depth(size_t d)
		{
			m_depth = d;
		}


	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "BlockGaussSeidel2";}


		template<typename TMatrix>
		void GetNeighborhood_worker(const TMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices, std::vector<bool> &bVisited)
		{
			if(depth==0) return;
			size_t iSizeBefore = indices.size();
			for(typename TMatrix::const_row_iterator it = A.begin_row(node); it != A.end_row(node); ++it)
			{
				if(it.value() == 0) continue;
				if(bVisited[it.index()] == false)
				{

					bVisited[it.index()] = true;
					indices.push_back(it.index());
				}
			}

			if(depth==1) return;
			size_t iSizeAfter = indices.size();
			for(size_t i=iSizeBefore; i<iSizeAfter; i++)
				GetNeighborhood_worker(A, indices[i], depth-1, indices, bVisited);
		}

		template<typename TMatrix>
		void GetNeighborhood2(const TMatrix &A, size_t node, size_t depth,
				std::vector<size_t> &indices,
				std::vector<size_t> &levels,
				std::vector<bool> &bVisited,
				bool bResetVisitedFlags=true)
		{
			PROFILE_FUNC_GROUP("algebra");
			levels.clear();
			levels.push_back(0);

			indices.clear();

			bVisited[node] = true;
			indices.push_back(node);
			levels.push_back(indices.size());

			for(size_t d = 1; d < depth ; d++)
			{
				for(size_t ind = levels[d-1]; ind < levels[d]; ind++)
				{
					size_t i = indices[i];
					for(typename TMatrix::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
					{
						if(it.value() == 0) continue;
						if(bVisited[it.index()] == false)
						{

							bVisited[it.index()] = true;
							indices.push_back(it.index());
						}
					}
				}
				levels.push_back(indices.size());
			}

			if(bResetVisitedFlags)
				for(size_t i=0; i<indices.size(); i++)
					bVisited[indices[i]] = false;
		}

		//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type &A = *pOp;

			size_t N = A.num_rows();
			DenseMatrix<VariableArray2<smallmat_type> > tmpMat;


			m_ilut.clear();

			indices.resize(N);

			size_t maxSize = 0;
			PROGRESS_START(prog, N, "BlockGaussSeidel2: compute blocks");


			std::vector<bool> bVisited(N, false);
			std::vector<bool> bVisited2(N, false);

			std::vector<size_t> levels;
			for(size_t i=0; i<N; i++)
			{
				if(bVisited[i]) continue;
				PROGRESS_UPDATE(prog, i);

				indices[i].clear();

				GetNeighborhood(A, i, m_depth, indices[i], bVisited2);
				for(size_t j=0; j<indices[i].size(); j++)
					bVisited[indices[i][j]]=true;
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

		typedef typename matrix_type::const_row_iterator matrix_const_row_iterator;
		typedef typename matrix_type::row_iterator matrix_row_iterator;

	//	Postprocess routine
		virtual bool postprocess() {return true;}
		virtual bool supports_parallel() const { return true; }

		//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			matrix_type &A = *pOp;
			vector_type &x = c;
			x.set(0.0);
			vector_type b;
			b = d;

			DenseMatrix<VariableArray2<double> > Adense;
			DenseMatrix<VariableArray2<smallmat_type> > Atmp;
			CPUAlgebra::vector_type tmp, tmp2;
			PROGRESS_START(prog, x.size(), "BlockGaussSeidel2: step");
			for(size_t i=0; i<x.size(); i++)
			{
				PROGRESS_UPDATE(prog, i);
				// c = D^{-1}(b-Ax)
				// x = x + c
				if(indices[i].size() != 0)
					GetBlockGSCorrectionILUT(A, x, b, m_ilut[i], indices[i], tmp, tmp2);

			}
			PROGRESS_FINISH(prog);

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			c.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss ; ss << "BlockGaussSeidel2(depth = " << m_depth << ")";
			return ss.str();
		}

};


} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
