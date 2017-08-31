/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX_UTIL__
#define __H__UG__CPU_ALGEBRA__SPARSEMATRIX_UTIL__

#include "common/profiler/profiler.h"
#include "unsorted_sparse_vector.h"
#include "../small_algebra/small_algebra.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallel_matrix.h"
#include "lib_algebra/parallelization/parallel_vector.h"
#endif

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


	typename block_multiply_traits<typename A_type::value_type, typename B_type::value_type>::ReturnType ab;

	typedef UnsortedSparseVector<typename ABC_type::value_type> RowType;
	typedef typename RowType::iterator RowIterator;
	RowType row(C.num_cols());;

	std::vector<typename ABC_type::connection> con2;
	//typename C_type::value_type cvalue;

	typename ABC_type::connection c;

	typedef typename A_type::const_row_iterator cAiterator;
	typedef typename B_type::const_row_iterator cBiterator;
	typedef typename C_type::const_row_iterator cCiterator;

	// do
	// M_{ij} = \sum_kl A_{ik} * B_{kl} * C_{lj}
	for(size_t i=0; i < A.num_rows(); i++)
	{
		row.clear();
		cAiterator itAikEnd = A.end_row(i);
		for(cAiterator itAik = A.begin_row(i); itAik != itAikEnd; ++itAik)
		{
			if(itAik.value() == 0.0) continue;

			size_t k = itAik.index();
			cBiterator itBklEnd = B.end_row(k);
			for(cBiterator itBkl = B.begin_row(k); itBkl != itBklEnd; ++itBkl)
			{
				if(itBkl.value() == 0.0) continue;
				size_t l = itBkl.index();
				// ab = A_{ik} * B_{kl}
				AssignMult(ab, itAik.value(), itBkl.value());

				cCiterator itCljEnd = C.end_row(l);
				for(cCiterator itClj = C.begin_row(l); itClj != itCljEnd; ++itClj)
				{
					if(itClj.value() == 0.0) continue;
					AddMult( row (itClj.index() ), ab, itClj.value());
				}
			}
		}

		if(epsilonTruncation != 0.0)
		{
			double m=0;
			for(RowIterator it = row.begin(); it != row.end(); ++it)
			{
				double d = BlockNorm(it->value());
				if(d > m) m = d;
			}
			m *= epsilonTruncation;
			con2.clear();
			for(RowIterator it = row.begin(); it != row.end(); ++it)
				if( BlockNorm(it->value()) > m )
					con2.push_back(*it);
			M.set_matrix_row(i, &con2[0], con2.size());
		}
		else
			M.set_matrix_row(i, row.unsorted_raw_ptr(), row.num_connections());
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AddMultiplyOf:
//-------------------------
/**
 * \brief Calculates M += A*B*C.
 * \param M (in/out) Matrix M, M += A*B*C$
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * \param C (in) Matrix C
 *
 * Complete formula for calculating M=A*B*C:
 *  \f[
 *  	 M_{ij} += \sum_{kl} A_{ik} * B_{kl} * C_{lj}
 *  \f]
 */
template<typename ABC_type, typename A_type, typename B_type, typename C_type>
void AddMultiplyOf(ABC_type &M, const A_type &A, const B_type &B, const C_type &C, double epsilonTruncation=0.0)
{
	PROFILE_FUNC_GROUP("algebra");
	UG_ASSERT(C.num_rows() == B.num_cols() && B.num_rows() == A.num_cols(), "sizes must match");

	// check
	if(M.num_rows() != A.num_rows())
		UG_THROW("AddMultiplyOf: row sizes mismatch: M.num_rows = "<<
		         M.num_rows()<<", A.num_rows = "<<A.num_rows());
	if(M.num_cols() != C.num_cols())
		UG_THROW("AddMultiplyOf: column sizes mismatch: M.num_cols = "<<
		         M.num_cols()<<", C.num_cols = "<<C.num_cols());

	typename block_multiply_traits<typename A_type::value_type, typename B_type::value_type>::ReturnType ab;

	typedef UnsortedSparseVector<typename ABC_type::value_type> RowType;
	typedef typename RowType::iterator RowIterator;
	RowType row(C.num_cols());;

	std::vector<typename ABC_type::connection> con2;
	//typename C_type::value_type cvalue;

	typename ABC_type::connection c;

	typedef typename A_type::const_row_iterator cAiterator;
	typedef typename B_type::const_row_iterator cBiterator;
	typedef typename C_type::const_row_iterator cCiterator;

	// do
	// M_{ij} = \sum_kl A_{ik} * B_{kl} * C_{lj}
	for(size_t i=0; i < A.num_rows(); i++)
	{
		row.clear();
		cAiterator itAikEnd = A.end_row(i);
		for(cAiterator itAik = A.begin_row(i); itAik != itAikEnd; ++itAik)
		{
			if(itAik.value() == 0.0) continue;

			size_t k = itAik.index();
			cBiterator itBklEnd = B.end_row(k);
			for(cBiterator itBkl = B.begin_row(k); itBkl != itBklEnd; ++itBkl)
			{
				if(itBkl.value() == 0.0) continue;
				size_t l = itBkl.index();
				// ab = A_{ik} * B_{kl}
				AssignMult(ab, itAik.value(), itBkl.value());

				cCiterator itCljEnd = C.end_row(l);
				for(cCiterator itClj = C.begin_row(l); itClj != itCljEnd; ++itClj)
				{
					if(itClj.value() == 0.0) continue;
					AddMult( row (itClj.index() ), ab, itClj.value());
				}
			}
		}

		if(epsilonTruncation != 0.0)
		{
			double m=0;
			for(RowIterator it = row.begin(); it != row.end(); ++it)
			{
				double d = BlockNorm(it->value());
				if(d > m) m = d;
			}
			m *= epsilonTruncation;
			con2.clear();
			for(RowIterator it = row.begin(); it != row.end(); ++it)
				if( BlockNorm(it->value()) > m )
					con2.push_back(*it);
			M.add_matrix_row(i, &con2[0], con2.size());
		}
		else
			M.add_matrix_row(i, row.unsorted_raw_ptr(), row.num_connections());
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

	// types
	typedef typename A_type::const_row_iterator cAiterator;
	typedef typename B_type::const_row_iterator cBiterator;
	typedef UnsortedSparseVector<typename AB_type::value_type> RowType;

	RowType row(B.num_cols());

	// M_{ij} = \sum_k A_{ik} * B_{kj}
	for(size_t i=0; i < A.num_rows(); i++)
	{
		row.clear();
		cAiterator itAikEnd = A.end_row(i);
		for(cAiterator itAik = A.begin_row(i); itAik != itAikEnd; ++itAik)
		{
			if(itAik.value() == 0.0) continue;
			size_t k = itAik.index();

			cBiterator itBklEnd = B.end_row(k);
			for(cBiterator itBkj = B.begin(k); itBkj != itBklEnd; ++itBkj)
			{
				if(itBkj.value() == 0.0) continue;
				size_t j = itBkj.index();
				AddMult( row(j), itAik.value(), itBkj.value());
			}
		}

		M.set_matrix_row(i, row.unsorted_raw_ptr(), row.num_connection());
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
void MatAdd(matrix_type &M, number alpha1, const matrix_type &A, number alpha2, const matrix_type &B)
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

/**
 * \brief Calculates M = A + B
 * \param M (out) Matrix M, M = A + B
 * \param A (in) Matrix A
 * \param B (in) Matrix B
 * note: A and/or B may be equal to M.
 */
template<typename matrix_type>
void MatAddNonDirichlet(matrix_type &M, number alpha1, const matrix_type &A, number alpha2, const matrix_type &B)
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
		{
		criterator itA = A.begin_row(i), endA = A.end_row(i);
		criterator itB = B.begin_row(i), endB = B.end_row(i);

		// copy only A for dirichlet rows
		// -> create a pattern Pii
		typedef typename matrix_type::value_type value_type;
		const value_type &Aii = A(i,i);
		value_type Pii = 1.0;


		UG_ASSERT (GetRows(Aii)==GetRows(Pii), "Huhh: Numbers of rows does not match!");
 		for(size_t alpha = 0; alpha < (size_t) GetRows(Aii); ++alpha)
 		{
 			if (IsDirichletRow(A, i, alpha)) BlockRef(Pii, alpha, alpha) =  0.0;
 		}


		// proceed as usual
		while(itA != endA && itB != endB)
		{
			// add this value: Bij:=Pii*Bij
			value_type Bij(Pii*itB.value());
			// UG_LOG("Bij =" << Bij << "," << "Pii=" << Pii);

			if(itA.index() == itB.index())
			{
				c.dValue = alpha1 * itA.value() + alpha2 * Bij;
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
				c.dValue = Bij;
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
			value_type Bij(Pii*itB.value());

			c.dValue = Bij;
			c.dValue *= alpha2;
			c.iIndex = itB.index();
			++itB;
			con.push_back(c);
		}
	}
		M.set_matrix_row(i, &con[0], con.size());
	}
	M.defragment();
}

template<typename TSparseMatrix>
void GetNeighborhood_worker(const TSparseMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices, std::vector<bool> &bVisited)
{
	if(depth==0) return;
	size_t iSizeBefore = indices.size();
	typename TSparseMatrix::const_row_iterator itEnd = A.end_row(node);
	for(typename TSparseMatrix::const_row_iterator it = A.begin_row(node); it != itEnd; ++it)
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

template<typename TSparseMatrix>
void GetNeighborhood(const TSparseMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices, std::vector<bool> &bVisited, bool bResetVisitedFlags=true)
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
			bVisited[indices[i]] = false;
}

template<typename TSparseMatrix>
void GetNeighborhood(const TSparseMatrix &A, size_t node, size_t depth, std::vector<size_t> &indices)
{
	PROFILE_FUNC_GROUP("algebra");
	std::vector<bool> bVisited(max(A.num_cols(), A.num_rows()), false);
	GetNeighborhood(A, node, depth, indices, bVisited, false);
}


template<typename TSparseMatrix>
void MarkNeighbors(const TSparseMatrix &A, size_t node, size_t depth, std::vector<bool> &bVisited)
{
	typename TSparseMatrix::const_row_iterator itEnd = A.end_row(node);
	for(typename TSparseMatrix::const_row_iterator it = A.begin_row(node); it != itEnd; ++it)
	{
		if(it.value() == 0) continue;
		bVisited[it.index()] = true;
		MarkNeighbors(A, it.index(), depth, bVisited);
	}
}

template<typename TSparseMatrix>
void GetNeighborhoodHierachy_worker(const TSparseMatrix &A, size_t node, size_t depth, size_t maxdepth, std::vector< std::vector<size_t> > &indices, std::vector<bool> &bVisited)
{
	size_t iSizeBefore = indices[depth].size();
	typename TSparseMatrix::const_row_iterator itEnd = A.end_row(node);
	for(typename TSparseMatrix::const_row_iterator it = A.begin_row(node); it != itEnd; ++it)
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

template<typename TSparseMatrix>
void GetNeighborhoodHierachy(const TSparseMatrix &A, size_t node, size_t depth, std::vector< std::vector<size_t> > &indices, std::vector<bool> &bVisited,
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
			typename TSparseMatrix::const_row_iterator itEnd = A.end_row(k);
			for(typename TSparseMatrix::const_row_iterator it = A.begin_row(k); it != itEnd; ++it)
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


template<typename TSparseMatrix>
void GetNeighborhoodHierachy(const TSparseMatrix &A, size_t node, size_t depth, std::vector< std::vector<size_t> > &indices)
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
template<typename TSparseMatrix>
bool IsCloseToBoundary(const TSparseMatrix &A, size_t node, size_t distance)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
	if(distance == 0) return A.is_isolated(node);
	bool bFound = false;
	iterator itEnd = A.end_row(node);
	for(iterator itA = A.begin_row(node); itA != itEnd && !bFound; ++itA)
		bFound = IsCloseToBoundary(A, itA.index(), distance-1);

	return bFound;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * set value for row for entry (i,alpha).
 * \param A (in) Matrix A
 * \param i (in) row to set
 * \param alpha the alpha index
 * \param val the value to be set
 */
template <typename TSparseMatrix>
void SetRow(TSparseMatrix &A, size_t i, size_t alpha, number val = 0.0)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	typedef typename TSparseMatrix::value_type value_type;
	iterator itEnd = A.end_row(i);
	for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
	{
		value_type& block = conn.value();
		// block : 1x1 (CPU=1), 3x3 (CPU=3)
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
		{
			BlockRef(block, alpha, beta) = val;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * set value for col for entry (i,alpha).
 * \param A (in) Matrix A
 * \param i (in) col to set
 * \param alpha the alpha index
 * \param val the value to be set
 */
template <typename TSparseMatrix>
void SetCol(TSparseMatrix &A, size_t i, size_t alpha, number val = 0.0)
{
	typedef typename TSparseMatrix::value_type value_type;
	for(size_t row = 0; row != A.num_rows(); ++row)
	{
		value_type& block = A(row, i);
		// block : 1x1 (CPU=1), 3x3 (CPU=3)
		for(size_t beta = 0; beta < (size_t) GetRows(block); ++beta)
		{
			BlockRef(block, beta, alpha) = val;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetRow:
//-------------------------
/**
 * set value for (block-)row i.
 * \param A (in) Matrix A
 * \param i (in) row to scales
 * \param val (in) value to be set
 */
template <typename TSparseMatrix>
void SetRow(TSparseMatrix& A, size_t i, number val = 0.0)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	iterator itEnd = A.end_row(i);
	for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
		conn.value() = val;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ScaleRow:
//-------------------------
/**
 * scales (block-)row i.
 * \param A (in) Matrix A
 * \param i (in) row to scales
 * \param fac (in) Scaling factor
 */
template <typename TSparseMatrix>
void ScaleRow(TSparseMatrix& A, size_t i, number fac)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	iterator itEnd = A.end_row(i);
	for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
		conn.value() *= fac;
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
template <typename TSparseMatrix>
void SetDirichletRow(TSparseMatrix &A, size_t i, size_t alpha)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	typedef typename TSparseMatrix::value_type value_type;
	BlockRef(A(i,i), alpha, alpha) = 1.0;

	iterator itEnd = A.end_row(i);
	for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
	{
		value_type& block = conn.value();
		// block : 1x1 (CPU=1), 3x3 (CPU=3)
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
		{
			if(conn.index() != i) BlockRef(block, alpha, beta) = 0.0;
			else if(beta != alpha) BlockRef(block, alpha, beta) = 0.0;
		}
	}
}

//! Evaluates 'true', iff corresponding row is Dirichlet
template <typename TSparseMatrix>
bool IsDirichletRow(const TSparseMatrix &A, size_t i, size_t alpha)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
	typedef typename TSparseMatrix::value_type value_type;

	// no Dirichlet row,
	if (BlockRef(A(i,i), alpha, alpha) != 1.0) return false;

	// check, if row sum equals 1
	number sum=0.0;
	iterator itEnd = A.end_row(i);
	for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
	{
		const value_type& block = conn.value();
		// block : 1x1 (CPU=1), 3x3 (CPU=3)
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta)
		{
			sum += BlockRef(block, alpha, beta) * BlockRef(block, alpha, beta);
		}
	}
	return (sum==1.0);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SetDirichletRow:
//-------------------------
/**
 * set Dirichlet row for block i.
 * \param A (in) Matrix A
 * \param i (in) row to set dirichlet, that is A(i,i) = 1.0, A(i,k) = 0.0 for all k != i.
 */
template <typename TSparseMatrix>
void SetDirichletRow(TSparseMatrix& A, size_t i)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	typedef typename TSparseMatrix::value_type value_type;
	A(i,i) = 1.0;
	iterator itEnd = A.end_row(i);
	for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
	{
		value_type& block = conn.value();
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
template <typename TSparseMatrix>
void SetDirichletRow(TSparseMatrix& A, const std::vector<size_t> vIndex)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	typedef typename TSparseMatrix::value_type value_type;
	std::vector<size_t>::const_iterator iter = vIndex.begin();
	std::vector<size_t>::const_iterator iterEnd = vIndex.end();

	for(; iter < iterEnd; ++iter)
	{
		const size_t i = *iter;
		UG_ASSERT(i < A.num_rows(), "Index to large in index set.");

		A(i,i) = 1.0;
		iterator itEnd = A.end_row(i);
		for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
		{
			value_type& block = conn.value();
			if(conn.index() != i) block = 0.0;
		}
	}
}

template<typename TSparseMatrix, class TOStream>
void SerializeMatrix(TOStream &buf, const TSparseMatrix &A)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
	Serialize(buf, A.num_rows());
	Serialize(buf, A.num_cols());

	for(size_t i=0; i < A.num_rows(); i++)
	{
		size_t num_connections = A.num_connections(i);

		// serialize number of connections
		Serialize(buf, num_connections);

		iterator itEnd = A.end_row(i);
		for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
		{
			// serialize connection
			Serialize(buf, conn.index());
			Serialize(buf, conn.value());
		}
	}
}

template <typename TSparseMatrix, class TIStream>
void DeserializeMatrix(TIStream& buf, TSparseMatrix &A)
{
	size_t numRows, numCols, num_connections;

	Deserialize(buf, numRows);
	Deserialize(buf, numCols);
	A.resize_and_clear(numRows, numCols);

	std::vector<typename TSparseMatrix::connection> con; con.reserve(16);

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

template<typename TSparseMatrix>
void ScaleSparseMatrixCommon(TSparseMatrix &A, double d)
{
	typedef typename TSparseMatrix::row_iterator iterator;
	for(size_t r=0; r < A.num_rows(); r++)
	{
		iterator itEnd = A.end_row(r);
		for(iterator it = A.begin_row(r); it != itEnd; ++it)
			it.value() *= d;
	}
}


template<typename TSparseMatrix, typename vector_t>
bool AxpyCommonSparseMatrix(const TSparseMatrix &A, vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
//	UG_ASSERT(cols == x.size(), "x: " << x << " has wrong length (should be " << cols << "). A: " << *this);
//	UG_ASSERT(rows == res.size(), "res: " << x << " has wrong length (should be " << rows << "). A: " << *this);

	if(alpha1 == 0.0)
	{
		for(size_t i=0; i < A.num_rows(); i++)
		{
			iterator conn = A.begin_row(i);
			iterator itEnd = A.end_row(i);
			if(conn  == itEnd) continue;
			MatMult(dest[i], beta1, conn.value(), w1[conn.index()]);
			for(++conn; conn != itEnd; ++conn)
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
	

template<typename TSparseMatrix, typename vector_t>
bool Axpy_transposedCommonSparseMatrix(const TSparseMatrix &A, vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
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
		iterator itEnd = A.end_row(i);
		for(iterator conn = A.begin_row(i); conn != itEnd; ++conn)
		{
			if(conn.value() != 0.0)
				// dest[conn.index()] += beta1 * conn.value() * w1[i];
				MatMultTransposedAdd(dest[conn.index()], 1.0, dest[conn.index()], beta1, conn.value(), w1[i]);
		}
	}
	return true;
}


/// returns the number of non-zeroes (!= number of connections)
template<typename TSparseMatrix>
size_t GetNNZs(const TSparseMatrix &A)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
	size_t m=0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		iterator itEnd = A.end_row(i);
		for(iterator it = A.begin_row(i); it != itEnd; ++it)
			if(it.value() != 0.0) m++;
	}
	return m;
}

/// returns max number of non-zero connections in rows
template<typename TSparseMatrix>
size_t GetMaxConnections(const TSparseMatrix &A)
{
	typedef typename TSparseMatrix::const_row_iterator iterator;
	size_t m=0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		//if(m < A.num_connections(i)) m = A.num_connections(i);

		size_t n=0;
		iterator itEnd = A.end_row(i);
		for(iterator it = A.begin_row(i); it != itEnd; ++it)
			if(it.value() != 0.0) n++;
		if(m < n) m = n;
	}

	return m;
}


template<typename TSparseMatrix>
bool CheckRowIterators(const TSparseMatrix &A)
{
#ifndef NDEBUG
	bool iIter=0;

	for(size_t i=0; i<A.num_rows(); i++)
	{
		if (A.nrOfRowIterators[i] !=0) UG_LOG ("CheckRowIterators: Failed for row " << i << std::endl);
		iIter += A.nrOfRowIterators[i];
	}
	return iIter==0;
#else
	return true;
#endif
}

template<typename TSparseMatrix>
bool CheckDiagonalInvertible(const TSparseMatrix &A)
{
#ifndef NDEBUG
	typedef typename block_traits<typename TSparseMatrix::value_type>::inverse_type inverse_type;
	bool bsucc=true;
	inverse_type inv;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		bool b = GetInverse(inv, A(i,i));
		if(!b)
		{
			UG_LOG("WARNING: entry " << i << " = " << A(i,i) << " not invertible\n");
			bsucc = false;
		}
	}
	return bsucc;
#else
	return true;
#endif
}

template<typename TVector>
bool CheckVectorInvertible(const TVector &v)
{
#ifndef NDEBUG
	typedef typename block_traits<typename TVector::value_type>::inverse_type inverse_type;
	bool bsucc=true;
	inverse_type inv;

	for(size_t i=0; i<v.size(); i++)
	{
		bool b = GetInverse(inv, v[i]);
		if(!b)
		{
			UG_LOG("WARNING: entry " << i << " = " << v[i] << " not invertible\n");
			bsucc = false;
		}
	}
	return bsucc;
#else
	return true;
#endif
}

#ifdef UG_PARALLEL
template<typename TSparseMatrix>
bool CheckDiagonalInvertible(const ParallelMatrix<TSparseMatrix> &m)
{
	return AllProcsTrue(CheckDiagonalInvertible((TSparseMatrix&)m), m.layouts()->proc_comm());
}

template<typename TVector>
bool CheckVectorInvertible(const ParallelVector<TVector> &v)
{
	return AllProcsTrue(CheckVectorInvertible((TVector&)v), v.layouts()->proc_comm());
}
#endif


template<typename TSparseMatrix>
struct DenseMatrixFromSparseMatrix
{
	typedef DenseMatrix<VariableArray2<typename TSparseMatrix::value_type> > type;
};

template<typename TSparseMatrix>
typename DenseMatrixFromSparseMatrix<TSparseMatrix>::type &
GetDenseFromSparse(typename DenseMatrixFromSparseMatrix<TSparseMatrix>::type &A, const TSparseMatrix &S)
{

	typedef typename TSparseMatrix::const_row_iterator sparse_row_iterator;
	size_t numRows = S.num_rows();
	size_t numCols = S.num_cols();
//	PROGRESS_START(prog, n, "GetDenseFromSparse " << n);
	A.resize(numRows, numCols);
	for(size_t r=0; r<numRows; r++)
		for(size_t c=0; c<numCols; c++)
			A(r, c) = 0;
	for(size_t r=0; r<numRows; r++)
	{
//		PROGRESS_UPDATE(prog, r);
		sparse_row_iterator itEnd = S.end_row(r);
		for(sparse_row_iterator it = S.begin_row(r); it != itEnd; ++it)
			A(r, it.index()) = it.value();
	}
//	PROGRESS_FINISH(prog);
	return A;
}

template<typename TSparseMatrix>
size_t GetDoubleSize(const TSparseMatrix &S)
{
	const size_t nrOfRows = block_traits<typename TSparseMatrix::value_type>::static_num_rows;
	UG_COND_THROW(nrOfRows != block_traits<typename TSparseMatrix::value_type>::static_num_cols, "only square matrices supported");
	return S.num_rows() * nrOfRows;
}

template<typename TDoubleType, typename TSparseMatrix>
void GetDoubleFromSparseBlock(TDoubleType &A, const TSparseMatrix &S)
{
	const size_t nrOfRows = block_traits<typename TSparseMatrix::value_type>::static_num_rows;
	for(size_t r=0; r<S.num_rows(); r++)
		for(typename TSparseMatrix::const_row_iterator it = S.begin_row(r); it != S.end_row(r); ++it)
		{
			size_t rr = r*nrOfRows;
			size_t cc = it.index()*nrOfRows;
			for(size_t r2=0; r2<nrOfRows; r2++)
					for(size_t c2=0; c2<nrOfRows; c2++)
					  A(rr + r2, cc + c2) = BlockRef(it.value(), r2, c2);
		}
}

template<typename TDenseType, typename TSparseMatrix>
size_t GetDenseDoubleFromSparse(TDenseType &A, const TSparseMatrix &S)
{
	size_t N = GetDoubleSize(S);
	A.resize(0,0);
	A.resize(N, N);
	GetDoubleFromSparseBlock(A, S);
	return N;
}

template<typename TDoubleSparse, typename TSparseMatrix>
size_t GetDoubleSparseFromBlockSparse(TDoubleSparse &A, const TSparseMatrix &S)
{
	size_t N = GetDoubleSize(S);

	A.resize_and_clear(N, N);
	GetDoubleFromSparseBlock(A, S);
	A.defragment();
	return N;
}




/// @}
} // end namespace ug


#endif
