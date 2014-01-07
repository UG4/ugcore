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


#include "../plugins/experimental/amg/src/rsamg/strong_connections.h"

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

template<typename matrix_type>
void GetStrongConnectedGraph(matrix_type &A, cgraph &g, double theta)
{

	cgraph graph;
	graph.resize(A.num_rows());

	for(size_t i=0; i< A.num_rows(); i++)
	{
		// Skip isolated node
		if(A.is_isolated(i)){
		//	UG_LOG(i << "is A-isolated!" << std::endl);
			continue;
		}

		// Eliminate positive row entries
		SparseVector<typename matrix_type::value_type> Ri(A.num_cols());
		PositiveOffDiagCheck<typename matrix_type::value_type> positiveEntries;
		EliminateEntries(positiveEntries, A, i, Ri);

		// Find strong connections
		FindStrongNegativeConnections(Ri, i, theta, graph);
	}

	g.resize(A.num_rows());
	std::vector<bool> b;
	std::vector<size_t> toReset;
	b.resize(graph.size(), false);
	for(size_t i=0; i< graph.size(); i++)
	{
		toReset.clear();
		for(cgraph::row_iterator it = graph.begin_row(i); it != graph.end_row(i); ++it)
		{
			int j = *it;
			g.set_connection(i, j);
			for(cgraph::row_iterator jt = graph.begin_row(j); jt != graph.end_row(j); ++jt)
			{
				size_t k = *jt;
				if(A.has_connection(i, k) == false) continue;
				if(b[k] == true) g.set_connection(i, k);
				else
				{
					b[k] = true;
					toReset.push_back(k);
				}
			}
		}
		for(size_t k=0; k<toReset.size(); k++)
			b [ toReset[k] ] = false;
	}
}
/*
bool DoesInterpolationOnlyFromStrongConnected
	(cgraph &strong, matrix_type &A, std::map<size_t, size_t> &newToOld, matrix_type &P, size_t i)
{
	size_t newIndex;
	if(P.num_connections(i) == 1) return false;

	std::vector<size_t> parentNodes;
	for(matrix_row_iterator it = P.begin_row(i); it != P.end_row(i); ++it)
	{
		size_t newIndex = it.index();
		size_t oldIndex = newToOld[newIndex];
		parentNodes.push_back(oldIndex);
	}

	for(size_t j=0; j<parentNodes.size(); j++)
	{
		size_t k = parentNodes[j];
		if(strong.has_connection(i, k) == false)
			return false;
	}

	return true;
}*/

template<typename matrix_type>
void GetStrongConnectedComponents(matrix_type &A, std::vector<int> &inComp, size_t maxCompSize)
{
	cgraph strongG;
	GetStrongConnectedGraph(A, strongG, 0.1);

	std::vector<bool> bvisited(A.num_rows(), false);

	inComp.resize(A.num_rows());
	std::vector<std::vector<size_t> > allComps;

	//if(DoesInterpolationOnlyFromStrongConnected) continue;
	//std::vector<bool> bDoesInterpolationOnlyFromStrongConnected;


//	for(size_t i=0; i<A.num_rows(); i++)
//	{
//		if(bDoesInterpolationOnlyFromStrongConnected[i]) { inComp[i] = -1; }
//		if(A.is_isolated(i)) { inComp[i] = -2; }
//	}

	size_t components = 0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(A.is_isolated(i)) continue;
		if(bvisited[i]) continue;

		size_t pos =0;
		std::vector<size_t> comp;
		comp.push_back(i);
//		std::set<size_t> connectedwithComponent;

		while(pos < comp.size())
		{
			size_t k = comp[pos];
			pos++;
			if(bvisited[k]) { assert(inComp[k] == components); continue; }

			bvisited[k] = true;
			inComp[k] = components;

			//if(bDoesInterpolationOnlyFromStrongConnected[k]) continue;


			for(cgraph::row_iterator it = strongG.begin_row(i); it != strongG.end_row(i); ++it)
			{
				size_t j=*it;

				if(bvisited[j])
				{
					// -1 is special group of node for bDoesInterpolationOnlyFromStrongConnected nodes
					if(inComp[j] == -1)
					{
						bvisited[j] = false;
						comp.push_back(k);
					}
					continue;

					/*size_t otherComp = inComp[j];
					// check if other component
					if(otherComp == components) // same component
						continue;
					else
						connectedwithComponent.insert(otherComp);*/
				}
				else
					comp.push_back(k);
			}
			if(comp.size() > maxCompSize)
			{
				// we hit the component limit. mark rest as visited.
				while(pos < comp.size())
					bvisited[comp[pos++]] = true;
				break;
			}
		}

		allComps.push_back(comp);
		components++;
	}

	for(size_t i=0; i<components; i++)
	{
		UG_LOG("COMPONENT " << i << "\n");
		for(size_t j=0; j<allComps[i].size(); j++)
		{
			UG_LOG(allComps[i][j] << " ");
		}
		UG_LOG("\n");
	}
}

/*
void bla(matrix_type &A, matrix_type &P)
{
	std::map<size_t, size_t> newToOld;
	for(size_t i=0; i<P.num_rows(); i++)
		if(P.num_connections() == 1)
		{
			size_t j = P.begin_row(i).index();
			newToOld[j] = i;
		}


	cgraph strongG;
	GetStrongConnectedGraph(A, strongG);


	std::vector<bool> diofgc(A.num_rows(), false);
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(A.is_isolated(i)) continue;

		diofgc[i] = DoesInterpolationOnlyFromStrongConnected(strongG, A, newToOld, P, i);
	}
}*/




/*void bla(matrix_type &A)
{
	cgraph graph;
	GetStrongConnectedGraph(A, graph);



}*/

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

			for(size_t i=0; i<N; i++)
			{
				if(bVisited[i]) continue;
				PROGRESS_UPDATE(prog, i);

				indices[i].clear();

				GetNeighborhood(A, i, m_depth, indices[i], bVisited2);
				for(size_t j=0; j<indices[i].size(); j++)
					bVisited[j]=true;
				//indices[i].push_back(i);

				Aloc[i] = new CPUAlgebra::matrix_type;

				GetSliceSparse(A, indices[i], *Aloc[i]);

//				DenseMatrix<VariableArray2<smallmat_type> > A2;
//				A2.resize(indices[i].size(), indices[i].size());

//				GetLocalMatrix(A, A2, &indices[i][0], &indices[i][0]);

//				UG_LOG("\n\nindices " << i << " ( " << indices[i].size() << ")\n");
//				PrintVector(indices[i], "indices");
//				UG_LOG("A2 " << i << ":\n" << JuliaString(A2, "A2") << "\n");


//				Aloc[i]->print("Aloc");

				m_ilut[i] = new ILUTPreconditioner<CPUAlgebra> (0.0);
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
