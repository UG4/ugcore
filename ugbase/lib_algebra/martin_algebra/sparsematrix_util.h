/*
 *  sparsematrix_util.h
 *  flexamg
 *
 *  Created by Martin Rupp on 11.06.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__SPARSEMATRIX_UTIL__
#define __H__UG__MARTIN_ALGEBRA__SPARSEMATRIX_UTIL__

namespace ug
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createAsMultiplyOf:
//-------------------------
//! Calculates *this = A B C. posInConnections only needed for speedup (has to be -1 forall i).
//! other possibility: search in vector con for con[i].iIndex == indexTo about 3-4x slower (3d)
//! @param  A
//! @param  B
//! @param  C
//! @param posInConnections		array of size B.getLength() for speedup of neighbor-neighbor-calculation inited with -1.
template<typename ABC_type, typename A_type, typename B_type, typename C_type>
void CreateAsMultiplyOf(ABC_type &M, const A_type &A, const B_type &B, const C_type &C, int *posInConnections)
{
	UG_ASSERT(C.num_rows() == B.num_cols() && B.num_rows() == A.num_cols(), "sizes must match");
	//cout << endl << " Creating Galerkin Matrix_type..." << endl;

	// create output matrix M
	M.create(A.num_rows(), C.num_cols());

	// speedup with array posInConnections, needs n memory
	// posInConnections[i]: index in the connections for current row.
	// has to be -1 for all nodes

	// tried this also with std::map, but took 1511.53 ms instead of 393.972 ms

	bool bOwnMem = false;
	if(posInConnections == NULL)
	{
		posInConnections = new int[C.num_cols()];
		for(size_t i=0; i<C.num_cols(); i++) posInConnections[i] = -1; // memset(posInConnections, -1, sizeof(int)*C.num_cols());
		bOwnMem = true;
		//assert(0);
	}
#ifndef NDEBUG
	else
	{
		for(size_t i=0; i<C.num_cols(); i++)
			UG_ASSERT(posInConnections[i] == -1, "posInConnections[" << i << "] has to be -1, but is " << posInConnections[i] << ".");
	}
#endif

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
						//TODO
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

	if(bOwnMem)
		delete[] posInConnections;
}





//!
//! @param node
//! @param depth
template<typename T>
void GetNeighborhood(SparseMatrix<T> &A, size_t node, size_t depth, vector<size_t> &indices, int *posInConnections)
{
	// do this with a map???
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

	if(posInConnections)
	{
		for(size_t i=0; i<indices.size(); i++)
			posInConnections[indices[i]] = -1;
	}

	sort(indices.begin(), indices.end());
}




template<typename T>
bool IsCloseToBoundary(const SparseMatrix<T> &A, size_t node, size_t distance)
{
	if(distance == 0) return A.isUnconnected(node);
	bool bFound = false;
	for(typename SparseMatrix<T>::cRowIterator itA = A.beginRow(node); !itA.isEnd() && !bFound; ++itA)
		bFound = IsCloseToBoundary(A, (*itA).iIndex, distance-1);

	return bFound;
}


/// set Dirichlet row for entry (i,alpha)
/// assumes, that diagonal block (i,i) exists
// TODO: This should also work, if diag block does not exist
template <typename T>
void SetDirichletRow(SparseMatrix<T>& A, size_t i, size_t alpha)
{
	for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
	{
		typename SparseMatrix<T>::entry_type& block = (*conn).dValue;
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
		{
			if((*conn).iIndex == i && beta == alpha) BlockRef(block, alpha, beta) = 1.0;
			else BlockRef(block, alpha, beta) = 0.0;
		}
	}
}

/// set Dirichlet row for block i
/// assumes, that diagonal block (i,i) exists
// TODO: This should also work, if diag block does not exist
template <typename T>
void SetDirichletRow(SparseMatrix<T>& A, size_t i)
{
	for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(i); !conn.isEnd(); ++conn)
	{
		typename SparseMatrix<T>::entry_type& block = (*conn).dValue;
		if((*conn).iIndex == i) block = 1.0;
		else block = 0.0;
	}
}


} // end namespace ug


#endif
