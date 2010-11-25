/**
 * \file amg_impl.h
 *
 * \author Martin Rupp
 *
 * \date 16.11.2010
 *
 * implementation file for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__

//#include "sparsematrix_util.h"

#include "amg_nodeinfo.h"
#include "stopwatch.h"


namespace ug{

//#define GRAPH_WITH_LOCAL_INVERSE


#define AMG_WRITE_MATRICES_PATH "/Users/mrupp/matrices/AMG_"
#define AMG_WRITE_MATRICES_MAX (200*200)

#define AMG_WRITE_GRAPH

#define LATE_COARSE_SOLVER // do coarsening down to 10 nodes.


#if 0

//#define AMG_WRITE_GRAPH


//#define AMG_PRINT_INDIRECT

#define AMG_PRINT_GRAPH

#define AMG_WRITE_COARSENING

#define AMG_PRINT_COARSENING
#define AMG_PRINT_P
#define AMG_PRINT_R
#define AMG_PRINT_AH

#define AMG_PRINT_COARSEN_RATINGS
#define AMG_PRINT_COARSEN
#endif

inline double amg_diag_value(const double &d) { return d; }
inline double amg_offdiag_value(const double &d) { return d; }

template<typename T> inline double amg_diag_value(const T &d) { return BlockNorm(d); }
template<typename T> inline double amg_offdiag_value(const T &d) { return -BlockNorm(d); }

void get_testvector_constant(DenseVector<VariableArray1<double> > &testvector, std::vector<size_t> &N2, int i_index)
{
	testvector.resize(N2.size());
	for(size_t j=0; j<testvector.size(); j++)
		testvector[j] = 1.0;
}



/*
template<typename value_type, typename pos_type>
void get_testvector_value_special(value_type &value, const pos_type &myPos, const pos_type &referencePos)
{
	BlockRef(value, TcomponentOut) = myPos[TcomponentIn] - referencePos[TcomponentIn];
}


inline void get_testvector_xx(DenseVector<VariableArray1<double> > &testvector, std::vector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 0>(testvector, N2, i_index);
}
inline void get_testvector_xy(DenseVector<VariableArray1<double> > &testvector, std::vector<size_t> &N2, int i_index)
{
	return get_testvector_special<0, 1>(testvector, N2, i_index);
}
...
2d :
(x 0), (y 0), (0 x), (0 y)
3d :
(x 0 0) (0 y 0) (0 0 z)  ..?

*/
template<typename value_type>
struct s_interpolation
{
	size_t i_from;
	value_type value;
};

struct neighborstruct2
{
	FixedArray1<s_interpolation, 2> parents;
	double F;
};

struct neighborstruct_var
{
	std::vector<s_interpolation> parents;
	double F;
};

void GetH(size_t N1size, DenseMatrix<VariableArray2<value_type> > &H, DenseMatrix<VariableArray2<value_type> > &S,
		vector<value_type> &Dinv, vector<value_type> &SF)
{
	size_t i_index = N1size;
	UG_ASSERT(S.num_cols() == S.num_rows());
	size_t N = S.num_rows();
	Dinv.resize(N);

	// get Dinv = 1/Aii
	for(size_t j=0; j < N; j++)
		GetInverse(Dinv[j], S(j,j));

	// get S = 1-w D^{-1} A
	for(size_t r=0; r < N; r++)
	{
		double omegadiaginv = omega/S(r, r);
		for(size_t c=0; c < N; c++)
			S(r, c) = 1 - omegadiaginv*S(r,c);
	}

	// get SF = 1-wDF^{-1} A  (F-smoothing)
	SF.resize(N);
	double omegadiaginv = omega/S(i_index, i_index);
	for(size_t j=0; j < N; j++)
		SF[j] = 1 - omegadiaginv * S(i_index, j);

	// get S = SF S
	for(size_t c=0; c < N; c++)
	{
		double s = 0;
		for(size_t j=0; j < N; j++)
			s += SF[j] * S(j, c);
		S(i_index, c) = s;
	}

	// get H = S S^T
	// todo: use symmetric H.
	for(size_t r=0; r < N1size; r++)
		for(size_t c=r; c < N1size; c++)
		{
			double s=0;
			for(size_t j=0; j<N; j++)
				s += S(r, j) * Dinv[j] * S(c, j);
			H(r, c) = s;
		}

	// get h0
	for(size_t r=0; r < N1size+1; r++)
	{
		double s=0;
		for(size_t j=0; j<N; j++)
			s += S(r, j) * Dinv[j] * S(i_index, j);
		Hi[r] = s;
	}

}

// neighbors including i, neighbors sorted
void GetAggressiveCoarseningInterpolation(const matrix_type &A, size_t i,
		DenseMatrix<VariableArray2<value_type> > &H, DenseMatrix<VariableArray2<value_type> > &S,
		vector<value_type> &Dinv, vector<value_type> &SF)
{
	DenseMatrix<VariableArray2<value_type> > H, S;
	vector<value_type> Dinv, SF;

	DenseMatrix<VariableArray2<double> > KKT;
	DenseVector<VariableArray1<double> > rhs;
	DenseVector<VariableArray1<double> > q;
	DenseVector<VariableArray1<double> > t;

	// aggressive coarsening
	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(nodes[node].is_uninterpolateable() == false) continue;

		// get Neighborhood N1 and N2

		GetNeighborhood(A, i, onlyN1, onlyN2, posInConnections);
		N2 = onlyN1;
		N2.push_back(i);
		N2.append(onlyN2);

		size_t i_index = onlyN1.size();

		// get submatrix A on N2

		S.resize(N2.size(), N2.size());
		A.get(S, neighbors);

		H.resize(onlyN1.size()+1, onlyN1.size()+1);
		GetH(onlyN1.size(), KKT, D, S, h0);

		// get testvector
		testvector.resize(onlyN1.size()+1);
		get_testvector(testvector, onlyN1, testvector[i_index], i);

		 /*  KKT = 	( H     t ) ( q_i,nm )    ( -H e_i )
		  *			( t^T   0 ) ( lambda )  = ( t[i]  ) */

		for(size_t j=0; j < onlyN1.size(); j++)
			KKT(i_index, j) = testvector[j];
		for(size_t j=0; j < onlyN1.size(); j++)
			rhs[j] = -h0[j];

		tmp.resize(N);

		rhs[onlyN1.size()] = - testvector[onlyN1.size()];
	}

}

template<typename heap_type>
void GetPossibleParentNodes(const matrix_type &A, std::vector<std::vector<neighborstruct2> > possible_neighbors,
		std::vector<famg_nodeinfo> &nodes)
{
	int *posInConnections = new int[A.num_rows()];
	for(size_t i=0; i<A.num_rows(); i++) posInConnections[i] = -1;

	std::vector<size_t> neighbors, N1, N2;
	DenseMatrix<VariableArray2<value_type> > S;
	DenseVector<VariableArray1<value_type> > v;
	std::vector<matrix_type::value_type> D;
	std::vector<neighborstruct2> i_neighborpairs;

	possible_neighbors.resize(A.num_rows());

	for(size_t i=0; i<A.num_rows(); i++)
	{
		// get Neighborhood N1 and N2

		GetNeighborhood(A, i, onlyN1, onlyN2, posInConnections);

		if(onlyN1.size() == 0)
		{
			nodes[i].set_fine();
			continue;
		}
		N2 = onlyN1;
		N2.push_back(i);
		N2.append(onlyN2);

		size_t i_index = onlyN1.size();

		// get submatrix A on N2

		S.resize(N2.size(), N2.size());
		A.get(S, neighbors);

		H.resize(onlyN1.size(), onlyN1.size());
		GetH(onlyN1.size(), H, D, S, h0);

		// get testvector
		testvector.resize(onlyN1.size()+1)
		get_testvector(testvector, onlyN1, testvector[i_index], i);

		DenseMatrix<FixedArray2<double, 3, 3> > KKT;
		DenseVector<FixedArray1<double, 3> > rhs;
		DenseVector<FixedArray1<double, 3> > q;
		DenseVector<FixedArray1<double, 3> > t;

		tmp.resize(N);
		rhs[3] = - testvector[i_index];

		i_neighborpairs.clear();
		i_neighborpairs.reserve(((N-1)*N)/2);

		neighborstruct2 n;

		int i_min = -1;
		double f_min = 1e12;
		for(size_t n=0; n < onlyN1.size(); n++)
		{
			for(size_t m=n+1; m < onlyN1.size(); m++)
			{
				// set KKT matrix
				/*
				 *  KKT = 	( H     t ) ( q_i,nm )    ( -H e_i )
				 * 			( t^T   0 ) ( lambda )  = ( t[i]  )
				 *
				 *  H *q = H(n,n) * q_n + H(n, m) * q_m + H(n, i) * q_i
				 *         H(m,n) * q_n + H(m, m) * q_m + H(m, i) * q_i
				 */

				KKT(0, 0) = H(n, n);
				KKT(0, 1) = H(n, m);
				KKT(1, 0) = H(m, n);
				KKT(1, 1) = H(m, m);

				KKT(2, 0) = KKT(0, 2) = testvector[n];
				KKT(2, 1) = KKT(1, 2) = testvector[m];

				rhs[0] = - Hi[n];
				rhs[1] = - Hi[m];

				InverseMatMult(q, 1.0, KKT, rhs);

				// calc q^T H q
				n.F = 	q[0] * (KKT(0,0) * q[0] + KKT(0,1) * q[1] + H(n,0)) +
						q[1] * (KKT(1,0) * q[0] + KKT(1,1) * q[1] + H(m,0)) +
						/*1 */ (Hi[n] * q[0] + Hi[m] * q[1] + Hi[i_index]);

				if(n.F > delta) continue;
				if(n.F < f_min)
				{
					f_min = n.F;
					i_min = i_neighborpairs.size();
				}
				n.parents[from].from = onlyN1[n];
				n.parents[from].value = -q[0];
				n.parents[1].from = onlyN1[m];
				n.parents[1].value = -q[1];
				i_neighborpairs.push_back(n);

				// F = q^T H q
			}
		}
		possible_neighbors[i].clear();
		if(i_min != -1)
		{
			// minimal element is always first in possible_neighbors[i] list.
			swap(i_neighborpairs[0], i_neighborpairs[i_min]);

			for(size_t j=0; j<i_neighborpairs.size(); j++)
				if(omega*i_neighborpairs[j].F <= f_min)
					possible_neighbors[i].push_back(i_neighborpairs[j]);
		}
		else
		{
			//UG_ASSERT(0, "node has no parents :'-(");
			nodes[i].set_uninterpolateable();
		}
	}

	delete[] posInConnections;
}




int GetRating(const FixedArray1<s_interpolation, 2> &M, const std::vector<famg_nodeinfo> &nodes)
{
	if(nodes[M[0].from].is_fine() || rating[M[1]].from].is_fine())
		return -1;
	return 2 - (nodes[M[0].from].is_coarse() ? 0 : 1) - (nodes[M[1].from].is_coarse() ? 0 : 1);
}

template<typename vec_type>
int GetRating(const vec_type &M, const std::vector<famg_nodeinfo> &nodes)
{
	size_t rating = 0;
	for(size_t i=0; i<M.size(); ++i)
	{
		famg_nodeinfo &ninfo = nodes[M[i].from];
		if(ninfo.is_fine()) return -1;
		else if(ninfo.is_uninterpolateable())
			rating += 1000;
		else if(!ninfo.is_coarse()) rating++;
	}
	return rating;
}

template<typename heap_type>
void GetRatings(std::vector<std::vector<neighborstruct2> > &possible_neighbors,
		std::vector<famg_nodeinfo> &nodes, heap_type &heap)
{
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i].rating == 0)
		{
			nodes[i].rating = UpdateRating(i, possible_neighbors[i], nodes);
			heap.insert_item(i);
		}
	}
}

// updates the rating of node <node>
// and writes it to rating
// possible neighbors are always stored as follows:
//	- possible_neighbors[i][0] is the best currently known pair for interpolating i
//	- this can change
//	a) if a parent node gets fine -> removal of this pair
//	b) if a parent node gets coarse
//	returns true if value has been updates, otherwise false
template<typename neighborstruct>
bool UpdateRating(size_t node, std::vector<neighborstruct2> &PN, const std::vector<famg_nodeinfo> &nodes)
{
	if(nodes[node].is_coarse()) return;

	int mini = -1;
	int minrating = 10000;
	double minF = 1e12;
	for(size_t i=0; i<PN.size(); )
	{
		int irating = GetRating(PN[i].parents, nodes);

		if(irating < 0)
		{
			swap(PN[i], PN.back()); // remove this pair
			P.resize(P.size()-1);
			continue;
		}

		// choose the pairs with best minimization F.
		if(irating <= minrating && PN[i].F < minF)
		{
			minrating = irating;
			minF = PN[i].F;
			mini = i;
		}
		i++;
	}

	if(mini != -1)
	{
		if(mini != 0) swap(PN[0], PN[mini]);
		if(nodes[node].rating != minrating)
		{
			nodes[node].rating = minrating;
			return true;
		}
		else return false;
	}
	else
	{
		if(nodes[node].is_uninterpolateable())
			return false;
		else
		{
			nodes[node].set_uninterpolateable();
			return true;
		}
	}
}

template<typename heap_type>
void UpdateNeighbors(size_t node, std::vector<std::vector<neighborstruct2> > &possible_neighbors,
		std::vector<famg_nodeinfo> &nodes)
{
	for(typename matrix_type::cRowIterator conn = A.beginRow(node); !conn.isEnd(); ++conn)
	{
		size_t neigh = conn.index();
		if(node == neigh) continue;
		if(!nodes[neigh].is_valid_rating())
			continue;
		if(GetUpdatedRating(neigh, possible_neighbors[neigh], nodes))
		{
			if(nodes[neigh].is_uninterpolateable())
				heap.remove(neigh);
			else
				heap.update(neigh);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// createFAMGLevel:
//-------------------------
/**
 * create FAMG matrix R, P, and AH = R A P
 * \param AH
 * \param R
 * \param A
 * \param P
 * \param level
 */
template<typename TAlgebra>
void famg<TAlgebra>::create_FAMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
		SparseMatrix<double> &P, int level)
{
	std::vector<std::vector<neighborstruct2> > possible_neighbors;
	std::vector<size_t> newIndex;
	std::vector<int> rating;
	std::vector<bool> fine;
	std::vector<bool> coarse;
	size_t N = A.num_rows();
	fine.resize(N);
	coarse.resize(N);

	new_index.resize(N);
	for(size_t i=0; i<N; i++)
		newIndex[i] = -1;
	size_t nrOfCoarse=0;

	rating.resize(N);
	maxheap heap(N, &rating[0]);

	// get possible parent nodes
	GetPossibleParentNodes(A, possible_neighbors, rating);

	// calculate ratings (not precalculateable because of coarse/uninterpolateable)
	GetRatings(possible_neighbors, rating, heap);

	// do coarsening
	size_t unassigned = N;
	while(heap.height() != 0)
	{
		size_t i = heap.remove_max();
		// node i gets fine. update neighbors.
		nodes[i].set_fine();
		UpdateNeighbors(i, fine, coarse, possible_neighbors, rating, heap);

		// get parent pair, set as coarse (if not already done), update neighbors.
		neighborstruct2 &n = possible_neighbors[i][0];
		for(size_t j=0; j < n.size(); j++)
		{
			size_t node = n.parents[j].from;
			if(!nodes[node].is_coarse())
			{
				nodes[node].set_coarse();
				heap.remove(node);
				newIndex[node] = nrOfCoarse++;
				UpdateNeighbors(node, possible_neighbors, nodes, heap);
			}
			P(node, newIndex[node]) = 1.0;
			P(i, newIndex[node]) = n.parents[j].value;
		}
	}

	// handle uninterpolateable nodes
	// todo: add aggressive coarsening here.
	for(size_t i=0; i<N; i++)
	{
		UG_ASSERT(nodes[i].is_valid_rating() == false, "?");
		if(nodes[i].is_uninterpolateable())
		{
			// set as coarse
			newIndex[i] = nrOfCoarse++;
			P(i, newIndex[i]);
			//nodes[i].set_coarse();
		}
	}
}



} // namespace ug

#endif //  __H__LIB_ALGEBRA__AMG__FAMG_IMPL_H__
