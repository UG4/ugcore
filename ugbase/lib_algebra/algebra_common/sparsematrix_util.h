/**
 * \file sparsematrix_util.h
 *
 * \author Martin Rupp
 *
 * \date 11.06.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX_UTIL__
#define __H__UG__CPU_ALGEBRA__SPARSEMATRIX_UTIL__

#include "../small_algebra/small_algebra.h"
#include "common/profiler/profiler.h"

namespace ug
{

/// \addtogroup lib_algebra
///	@{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateAsMultiplyOf:
//-------------------------
/**
 * \brief Calculates M = A*B*C.
 * \param M (out) Matrix M, M = A*B*C$
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * \param C (in) Matrix C
 *
 * Complete formula for calculating M=A*B*C:
 *  \f[
 *  	 M_{ij} = \sum_{kl} A_{ik} * B_{kl} * C_{lj}
 *  \f]
 * Calculation is done on row-basis without
 * a temporary BC or AB matrix. This has shown to be much faster than
 * implementations with temporary matrices due to cache effects.
 * We also added an improved way of storing the results of the calculation:
 *  when we go through connections of B and C and want to add the connection
 *  (i, j) to M, we need to know if this connection already exists. For this we have
 *  an array posInConnections, needs n=A.num_rows() memory.
 *  posInConnections[i]: index in the connections for current row (if not in row: -1)
 *  tried this also with std::map, but took 1511.53 ms instead of 393.972 ms
 *  searching in the connections array is also slower
 */
template<typename ABC_type, typename A_type, typename B_type, typename C_type>
void CreateAsMultiplyOf(ABC_type &M, const A_type &A, const B_type &B, const C_type &C, double epsilonTruncation=0.0)
{
	PROFILE_FUNC_GROUP("algebra");
	UG_ASSERT(C.num_rows() == B.num_cols() && B.num_rows() == A.num_cols(), "sizes must match");

	// create output matrix M
	M.resize_and_clear(A.num_rows(), C.num_cols());


	std::vector<int> posInConnections(C.num_cols(), -1);

	// types
	std::vector<typename ABC_type::connection > con; con.reserve(16);
	std::vector<typename ABC_type::connection > con2; con2.reserve(16);

	typedef typename A_type::value_type avalue;
	typename block_multiply_traits<typename A_type::value_type, typename B_type::value_type>::ReturnType ab;
	//typename C_type::value_type cvalue;

	typename ABC_type::connection c;

	typedef typename A_type::const_row_iterator cAiterator;
	typedef typename B_type::const_row_iterator cBiterator;
	typedef typename C_type::const_row_iterator cCiterator;

	// do
	// M_{ij} = \sum_kl A_{ik} * B_{kl} * C_{lj}
	for(size_t i=0; i < A.num_rows(); i++)
	{
		con.clear();
		for(cAiterator itAik = A.begin_row(i); itAik != A.end_row(i); ++itAik)
		{
			if(itAik.value() == 0.0) continue;

			size_t k = itAik.index();
			for(cBiterator itBkl = B.begin_row(k); itBkl != B.end_row(k); ++itBkl)
			{
				if(itBkl.value() == 0.0) continue;
				size_t l = itBkl.index();
				// ab = A_{ik} * B_{kl}
				AssignMult(ab, itAik.value(), itBkl.value());

				for(cCiterator itClj = C.begin_row(l); itClj != C.end_row(l); ++itClj)
				{
					if(itClj.value() == 0.0) continue;
					size_t j = itClj.index();

					if(posInConnections[j] == -1)
					{
						// we havent visited node <indexTo>
						// so we need to add a Connection to the row
						// save the index of the connection in the row
						posInConnections[j] = con.size();
						c.iIndex = j;
						AssignMult(c.dValue, ab, itClj.value());
						con.push_back(c);
					}
					else
					{
						// we have visited this node before,
						// so we know the index of the connection
						// -> add a*b*c
						AddMult(con[posInConnections[j]].dValue, ab, itClj.value());
					}

				}
			}
		}

		// reset posInConnections to -1
		for(size_t j=0; j<con.size(); j++) posInConnections[con[j].iIndex] = -1;
		if(epsilonTruncation != 0.0)
		{
			double m=0;
			for(size_t j=0; j<con.size(); j++)
			{
				double d = BlockNorm(con[j].dValue);
				if(d > m) m = d;
			}
			m *= epsilonTruncation;
			con2.clear();
			for(size_t j=0; j<con.size(); j++)
				if( BlockNorm(con[j].dValue) > m )
					con2.push_back(con[j]);
			M.set_matrix_row(i, &con2[0], con2.size());
		}
		else
			M.set_matrix_row(i, &con[0], con.size());
	}

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateAsMultiplyOf:
//-------------------------
/**
 * \brief Calculates M = A*B.
 * \param M (out) Matrix M, M = A*B$
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * \f$ M_{ij} = \sum_k A_{ik} * B_{kj} \f$
 * For implementation details, see also CreateAsMultiplyOf(M, A, B, C).
 */
template<typename AB_type, typename A_type, typename B_type>
void CreateAsMultiplyOf(AB_type &M, const A_type &A, const B_type &B)
{
	PROFILE_FUNC_GROUP("algebra");
	UG_ASSERT(B.num_rows() == A.num_cols(), "sizes must match");

	// create output matrix M
	M.resize_and_clear(A.num_rows(), B.num_cols());

	std::vector<int> posInConnections(B.num_cols(), -1);

	// types
	std::vector<typename AB_type::connection > con; con.reserve(255);
	typename AB_type::connection c;
	typedef typename A_type::const_row_iterator cAiterator;
	typedef typename B_type::const_row_iterator cBiterator;

	// M_{ij} = \sum_k A_{ik} * B_{kj}
	for(size_t i=0; i < A.num_rows(); i++)
	{
		con.clear();
		for(cAiterator itAik = A.begin_row(i); itAik != itAik.end_row(i); ++itAik)
		{
			if(itAik.value() == 0.0) continue;
			size_t k = itAik.index();


			for(cBiterator itBkj = B.begin(k); itBkj != B.end_row(k); ++itBkj)
			{
				if(itBkj.value() == 0.0) continue;
				size_t j = itBkj.index();

				if(posInConnections[j] == -1)
				{
					// we havent visited node <indexTo>
					// so we need to add a Connection to the row
					// save the index of the connection in the row
					posInConnections[j] = con.size();
					c.iIndex = j;
					AssignMult(c.dValue, itAik.value(), itBkj.value());
					con.push_back(c);
				}
				else
				{
					// we have visited this node before,
					// so we know the index of the connection
					// -> add a*b*c
					AddMult(con[posInConnections[j]].dValue, itAik.value(), itBkj.value());
				}
			}
		}

		// reset posInConnections to -1
		for(size_t l=0; l<con.size(); l++) posInConnections[con[l].iIndex] = -1;
		// set Matrix_type Row in AH
		M.set_matrix_row(i, &con[0], con.size());
	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MatAdd:
//-------------------------
/**
 * \brief Calculates M = A + B
 * \param M (out) Matrix M, M = A + B
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * note: A and/or B may be equal to M.
 */
template<typename matrix_type>
void MatAdd(matrix_type &M, number &alpha1, const matrix_type &A, number &alpha2, const matrix_type &B)
{
	PROFILE_FUNC_GROUP("algebra");
	UG_ASSERT(A.num_rows() == B.num_rows() && A.num_cols() == B.num_cols(), "sizes must match");
	typedef typename matrix_type::const_row_iterator criterator;

	// create output matrix M
	if(&M != &A)
		M.resize_and_clear(A.num_rows(), A.num_cols());

	// types
	std::vector<typename matrix_type::connection > con; con.reserve(10);

	typename matrix_type::connection c;
	for(size_t i=0; i < A.num_rows(); i++)
	{
		con.clear();
		criterator itA = A.begin_row(i), endA = A.end_row(i);
		criterator itB = B.begin_row(i), endB = B.end_row(i);

		while(itA != endA && itB != endB)
		{
			if(itA.index() == itB.index())
			{
				c.dValue = alpha1 * itA.value() + alpha2 * itB.value();
				c.iIndex = itA.index();
				++itA; ++itB;
			}
			else if (itA.index() < itB.index())
			{
				c.dValue = itA.value();
				c.dValue *= alpha1;
				c.iIndex = itA.index();
				++itA;
			}
			else
			{
				c.dValue = itB.value();
				c.dValue *= alpha2;
				c.iIndex = itB.index();
				++itB;
			}
			con.push_back(c);
		}
		while(itA != endA)
		{
			c.dValue = itA.value();
			c.dValue *= alpha1;
			c.iIndex = itA.index();
			++itA;
			con.push_back(c);
		}
		while(itB != endB)
		{
			c.dValue = itB.value();
			c.dValue *= alpha2;
			c.iIndex = itB.index();
			++itB;
			con.push_back(c);
		}

		M.set_matrix_row(i, &con[0], con.size());
	}
	M.defragment();
}

template<typename TMatrix>
void GetNeighborhood_worker(const TMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices, std::vector<bool> &bVisited)
{
	if(depth==0) return;
	size_t iSizeBefore = indices.size();
	for(typename TMatrix::const_row_iterator it = A.begin_row(node); A.end_row(node); ++it)
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
void GetNeighborhood(const TMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices, std::vector<bool> &bVisited, bool bResetVisitedFlags=true)
{
	PROFILE_FUNC_GROUP("algebra");
	indices.clear();

	if(bVisited[node] == false)
	{
		bVisited[node] = true;
		indices.push_back(node);
	}
	GetNeighborhood_worker(A, node, depth, indices, bVisited);

	if(bResetVisitedFlags)
		for(size_t i=0; i<indices.size(); i++)
			bVisited[i] = false;
}

template<typename TMatrix>
void GetNeighborhood(const TMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices)
{
	PROFILE_FUNC_GROUP("algebra");
	std::vector<bool> bVisited(max(A.num_cols(), A.num_rows()), false);
	GetNeighborhood(A, node, depth, indices, bVisited, false);
}


template<typename TMatrix>
void MarkNeighbors(const TMatrix &A, size_t node, size_t depth, std::vector<bool> &bVisited)
{
	for(typename TMatrix::const_row_iterator it = A.begin_row(node); it != A.end_row(node); ++it)
	{
		if(it.value() == 0) continue;
		bVisited[it.index()] = true;
		MarkNeighbors(A, it.index(), depth, bVisited);
	}
}

template<typename TMatrix>
void GetNeighborhoodHierachy_worker(const TMatrix &A, size_t node, size_t depth, size_t maxdepth, std::vector< std::vector<size_t> > &indices, std::vector<bool> &bVisited)
{
	size_t iSizeBefore = indices[depth].size();
	for(typename TMatrix::const_row_iterator it = A.begin_row(node); it != A.end_row(node); ++it)
	{
		if(it.value() == 0) continue;
		if(bVisited[it.index()] == false)
		{
			bVisited[it.index()] = true;
			indices[depth].push_back(it.index());
		}
	}

	if(depth==maxdepth) return;
	size_t iSizeAfter = indices[depth].size();
	for(size_t i=iSizeBefore; i<iSizeAfter; i++)
		GetNeighborhoodHierachy_worker(A, indices[i], depth+1, maxdepth, indices, bVisited);
}

template<typename TMatrix>
void GetNeighborhoodHierachy(const TMatrix &A, size_t node, size_t depth, std::vector< std::vector<size_t> > &indices, std::vector<bool> &bVisited,
		bool bResetVisitedFlags=true)
{
	PROFILE_FUNC_GROUP("algebra");
	if(indices.size() != depth+1)
		indices.resize(depth+1);
	for(size_t i=0; i < depth+1; i++)
		indices[i].clear();

	bVisited[node] = true;
	indices[0].push_back(node);

	if(depth==0) return;

	for(size_t d = 0; d < depth; d++)
	{
		for(size_t i=0; i<indices[d].size(); i++)
		{
			size_t k = indices[d][i];
			for(typename TMatrix::const_row_iterator it = A.begin_row(k); it != A.end_row(k); ++it)
			{
				if(it.value() == 0) continue;
				if(bVisited[it.index()] == false)
				{
					bVisited[it.index()] = true;
					indices[d+1].push_back(it.index());
				}
			}

		}
	}

	if(bResetVisitedFlags)
		for(size_t i=0; i < depth+1; i++)
			for(size_t j=0; i<indices[j].size(); j++)
				bVisited[j] = false;
}


template<typename TMatrix>
void GetNeighborhoodHierachy(const TMatrix &A, size_t node, size_t depth, std::vector< std::vector<size_t> > &indices)
{
	std::vector<bool> bVisited(std::max(A.num_cols(), A.num_rows()), false);
	GetNeighborhoodHierachy(A, node, depth, indices, bVisited, false);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GetNeighborhood:
//-------------------------
/**
 * \brief gets the neighborhood of a node in the connectivity graph of a SparseMatrix.
 * \param A (in) Matrix A
 * \param node (in) the node where to start
 * \param depth (in) the depth of neighborhood. 0 = empty.
 * \param indices (out) the indices of the neighbors
 * \param posInConnections array to speed up computation. Has to be posInConnections[i] = 0 for all i=0..A.num_rows(). Can be NULL.
 * \note the node itself is only included if there is a connection from node to node.
  */
#if 0
template<typename T>
void GetNeighborhood(SparseMatrix<T> &A, size_t node, size_t depth, std::vector<size_t> &indices, int *posInConnections=NULL)
{
	// perhaps this is better with recursion
	indices.clear();



	vector<typename SparseMatrix<T>::const_row_iterator> iterators;
	iterators.reserve(depth);

	iterators.push_back( A.begin_row(node) );

	while(iterators.size() != 0)
	{
		if(iterators.back().isEnd())
			iterators.pop_back();
		else
		{
			size_t index = iterators.back().index();
			++iterators.back();
			if(iterators.size() < depth)
				iterators.push_back( A.begin_row(index) );
			else
			{
				size_t pos;
				if(posInConnections == NULL)
				{
					for(pos=0; pos<indices.size(); pos++)
						if(indices[pos] == index)
							break;
					if(pos == indices.size())
						indices.push_back(index);
				}
				else
				{
					pos = posInConnections[index];
					if(pos == -1)
					{
						pos = posInConnections[index] = indices.size();
						indices.push_back(index);
					}
				}
				// else (count etc.)
			}
		}
	}

	// reset posInConnections
	if(posInConnections)
	{
		for(size_t i=0; i<indices.size(); i++)
			posInConnections[indices[i]] = -1;
	}

	// sort indices
	sort(indices.begin(), indices.end());
}
#endif



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IsCloseToBoundary:
//-------------------------
/**
 * \brief determines if a node is close to a unconnected node in the connectivity graph of a SparseMatrix.
 * \param A (in) Matrix A
 * \param node (in) the node where to start
 * \param distance (in) up to which distance "close" is.
 * \return if there is a distance long path in graph(A) to an unconnected node, true. otherwise false.
  */
template<typename T>
bool IsCloseToBoundary(const T &A, size_t node, size_t distance)
{
	if(distance == 0) return A.is_isolated(node);
	bool bFound = false;
	for(typename T::const_row_iterator itA = A.begin_row(node); itA != A.end_row(node) && !bFound; ++itA)
		bFound = IsCloseToBoundary(A, itA.index(), distance-1);

	return bFound;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetDirichletRow:
//-------------------------
/**
 * set Dirichlet row for entry (i,alpha).
 * \param A (in) Matrix A
 * \param i (in) row to set dirichlet, that is A(i,i)(alpha, alpha) = 1.0, A(i,k)(alpha, beta) = 0.0 for all (k, beta) != (i, alpha)$.
 * \param alpha the alpha index
 */
template <typename T>
void SetDirichletRow(T &A, size_t i, size_t alpha)
{
	BlockRef(A(i,i), alpha, alpha) = 1.0;
	for(typename T::row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		typename T::value_type& block = conn.value();
		// block : 1x1 (CPU=1), 3x3 (CPU=3)
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
		{
			if(conn.index() != i) BlockRef(block, alpha, beta) = 0.0;
			else if(beta != alpha) BlockRef(block, alpha, beta) = 0.0;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetDirichletRow:
//-------------------------
/**
 * set Dirichlet row for block i.
 * \param A (in) Matrix A
 * \param i (in) row to set dirichlet, that is A(i,i) = 1.0, A(i,k) = 0.0 for all k != i.
 */
template <typename T>
void SetDirichletRow(T& A, size_t i)
{
	A(i,i) = 1.0;
	for(typename T::row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
	{
		typename T::value_type& block = conn.value();
		if(conn.index() != i) block = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetDirichletRow:
//-------------------------
/**
 * set Dirichlet row for block i.
 * \param[in/out] A Matrix A
 * \param[in] vIndex vector of row indices to set dirichlet, that is A(i,i) = 1.0, A(i,k) = 0.0 for all k != i.
 */
template <typename T>
void SetDirichletRow(T& A, const std::vector<size_t> vIndex)
{
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < A.num_rows(), "Index to large in index set.");

		A(i,i) = 1.0;
		for(typename T::row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			typename T::value_type& block = conn.value();
			if(conn.index() != i) block = 0.0;
		}
	}
}

template<typename T, class TOStream>
void SerializeMatrix(TOStream &buf, const T &A)
{
	Serialize(buf, A.num_rows());
	Serialize(buf, A.num_cols());

	for(size_t i=0; i < A.num_rows(); i++)
	{
		size_t num_connections = A.num_connections(i);

		// serialize number of connections
		Serialize(buf, num_connections);

		for(typename T::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			// serialize connection
			Serialize(buf, conn.index());
			Serialize(buf, conn.value());
		}
	}
}

template <typename T, class TIStream>
void DeserializeMatrix(TIStream& buf, T &A)
{
	size_t numRows, numCols, num_connections;

	Deserialize(buf, numRows);
	Deserialize(buf, numCols);
	A.resize_and_clear(numRows, numCols);

	std::vector<typename T::connection> con; con.reserve(16);

	for(size_t i=0; i < A.num_rows; i++)
	{
		Deserialize(buf, num_connections);

		con.resize(num_connections);

		for(size_t j=0; j<num_connections; j++)
		{
			Deserialize(buf, con[j].iIndex);
			Deserialize(buf, con[j].dValue);
		}
		A.set_matrix_row(i, &con[0], num_connections);
	}
	A.defragment();
}

template<typename T>
void ScaleSparseMatrixCommon(T &A, double d)
{
	for(size_t r=0; r < A.num_rows(); r++)
		for(typename T::row_iterator it = A.begin_row(r); it != A.end_row(r); ++it)
			it.value() *= d;
}


template<typename T, typename vector_t>
bool AxpyCommonSparseMatrix(const T &A, vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1)
{
//	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
//	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);

	if(alpha1 == 0.0)
	{
		for(size_t i=0; i < A.num_rows(); i++)
		{
			typename T::const_row_iterator conn = A.begin_row(i);
			if(conn  == A.end_row(i)) continue;
			MatMult(dest[i], beta1, conn.value(), w1[conn.index()]);
			for(++conn; conn != A.end_row(i); ++conn)
				// res[i] += conn.value() * x[conn.index()];
				MatMultAdd(dest[i], 1.0, dest[i], beta1, conn.value(), w1[conn.index()]);
		}
	}
	else if(&dest == &v1)
	{
		if(alpha1 != 1.0)
			for(size_t i=0; i < A.num_rows(); i++)
			{
				dest[i] *= alpha1;
				A.mat_mult_add_row(i, dest[i], beta1, w1);
			}
		else
			for(size_t i=0; i < A.num_rows(); i++)
				A.mat_mult_add_row(i, dest[i], beta1, w1);

	}
	else
	{
		for(size_t i=0; i < A.num_rows(); i++)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			A.mat_mult_add_row(i, dest[i], beta1, w1);
		}
	}

	return true;
}
	

template<typename T, typename vector_t>
bool Axpy_transposedCommonSparseMatrix(const T &A, vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1)
{
//	UG_ASSERT(rows == x.size(), "x: " << x << " has wrong length (should be " << rows << "). A: " << *this);
//	UG_ASSERT(cols == res.size(), "res: " << x << " has wrong length (should be " << cols << "). A: " << *this);

	if(&dest == &v1) {
		if(alpha1 == 0.0)
			dest.set(0.0);
		else if(alpha1 != 1.0)
			dest *= alpha1;
	}
	else if(alpha1 == 0.0)
		dest.set(0.0);
	else
		VecScaleAssign(dest, alpha1, v1);

	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(typename T::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(conn.value() != 0.0)
				// dest[conn.index()] += beta1 * conn.value() * w1[i];
				MatMultTransposedAdd(dest[conn.index()], 1.0, dest[conn.index()], beta1, conn.value(), w1[i]);
		}
	}
	return true;
}

/// @}
} // end namespace ug


#endif
