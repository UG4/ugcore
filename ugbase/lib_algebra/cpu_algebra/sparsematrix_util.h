/**
 * \file sparsematrix_util.h
 *
 * \author Martin Rupp
 *
 * \date 11.06.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__MARTIN_ALGEBRA__SPARSEMATRIX_UTIL__
#define __H__UG__MARTIN_ALGEBRA__SPARSEMATRIX_UTIL__

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
 */
template<typename ABC_type, typename A_type, typename B_type, typename C_type>
void CreateAsMultiplyOf(ABC_type &M, const A_type &A, const B_type &B, const C_type &C)
{
	UG_ASSERT(C.num_rows() == B.num_cols() && B.num_rows() == A.num_cols(), "sizes must match");

	// create output matrix M
	M.create(A.num_rows(), C.num_cols());

	// speedup with array posInConnections, needs n memory
	// posInConnections[i]: index in the connections for current row (if not in row: -1)
	// tried this also with std::map, but took 1511.53 ms instead of 393.972 ms
	// searching in the connections is also slower

	int *posInConnections = new int[C.num_cols()];
	for(size_t i=0; i<C.num_cols(); i++) posInConnections[i] = -1;

	// types
	vector<typename ABC_type::connection > con(255);

	typename A_type::entry_type a;
	typename block_multiply_traits<typename A_type::entry_type, typename B_type::entry_type>::ReturnType ab;
	typename C_type::entry_type cvalue;

	typename ABC_type::connection c;

	// do
	for(size_t i=0; i < A.num_rows(); i++)
	{
		con.clear();
		for(typename A_type::cRowIterator itA(A, i); !itA.isEnd(); ++itA)
		{
			if((*itA).dValue == 0.0) continue;
			a = (*itA).dValue;

			for(typename B_type::cRowIterator itB(B, (*itA).iIndex); !itB.isEnd(); ++itB)
			{
				if((*itB).dValue == 0.0) continue;
				AssignMult(ab, a, (*itB).dValue);

				for(typename C_type::cRowIterator itC(C, (*itB).iIndex); !itC.isEnd(); ++itC)
				{
					cvalue = (*itC).dValue;
					if(cvalue == 0.0) continue;
					size_t indexTo = (*itC).iIndex;

					if(posInConnections[indexTo] == -1)
					{
						// we havent visited node <indexTo>
						// so we need to add a Connection to the row
						// save the index of the connection in the row
						posInConnections[indexTo] = con.size();
						c.iIndex = indexTo;
						AssignMult(c.dValue, ab, cvalue);
						con.push_back(c);
					}
					else
					{
						// we have visited this node before,
						// so we know the index of the connection
						// -> add a*b*c
						AddMult(con[posInConnections[indexTo]].dValue, ab, cvalue);
					}

				}
			}
		}

		// reset posInConnections to -1
		for(size_t j=0; j<con.size(); j++) posInConnections[con[j].iIndex] = -1;
		// set Matrix_type Row in AH
		M.set_matrix_row(i, &con[0], con.size());
	}

	delete[] posInConnections;
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
template<typename T>
void GetNeighborhood(SparseMatrix<T> &A, size_t node, size_t depth, vector<size_t> &indices, int *posInConnections=NULL)
{
	// perhaps this is better with recursion
	indices.clear();

	vector<typename SparseMatrix<T>::cRowIterator> iterators;
	iterators.reserve(depth);

	iterators.push_back( A.beginRow(node) );

	while(iterators.size() != 0)
	{
		if(iterators.back().isEnd())
			iterators.pop_back();
		else
		{
			size_t index = (*iterators.back()).iIndex;
			++iterators.back();
			if(iterators.size() < depth)
				iterators.push_back( A.beginRow(index) );
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
bool IsCloseToBoundary(const SparseMatrix<T> &A, size_t node, size_t distance)
{
	if(distance == 0) return A.isUnconnected(node);
	bool bFound = false;
	for(typename SparseMatrix<T>::cRowIterator itA = A.beginRow(node); !itA.isEnd() && !bFound; ++itA)
		bFound = IsCloseToBoundary(A, (*itA).iIndex, distance-1);

	return bFound;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetDirichletRow:
//-------------------------
/**
 * set Dirichlet row for entry (i,alpha).
 * \param A (in) Matrix A
 * \param i (in) row to set dirichlet, that is A(i,i)(alpha, alpha) = 1.0, A(i,k)(alpha, beta) = 0.0 for all (k, beta) != (i, alpha)$.
   */
template <typename T>
void SetDirichletRow(SparseMatrix<T>& A, size_t i, size_t alpha)
{
	BlockRef(A(i,i), alpha, alpha) = 1.0;
	for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
	{
		typename SparseMatrix<T>::entry_type& block = (*conn).dValue;
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
		{
			if((*conn).iIndex != i) BlockRef(block, alpha, beta) = 0.0;
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
void SetDirichletRow(SparseMatrix<T>& A, size_t i)
{
	A(i,i) = 1.0;
	for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
	{
		typename SparseMatrix<T>::entry_type& block = (*conn).dValue;
		if((*conn).iIndex != i) block = 0.0;
	}
}


/// @}
} // end namespace ug


#endif
