/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#include "famg_nodeinfo.h"

#ifndef __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
#define __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__

namespace ug {


template<typename matrix_type>
void GetNeighborhood(matrix_type &A, size_t node, vector<size_t> &onlyN1)
{
	onlyN1.clear();
	for(typename matrix_type::cRowIterator it = A.beginRow(node); !it.isEnd(); ++it)
	{
		if(it.value() != 0.0)
			onlyN1.push_back(it.index());
	}
}

// adds the indices of the neighbors of a node to the vector onlyN1, if they have visited flag 0.
// node: neither onlyN1 nor bvisited get reset
template<typename matrix_type>
void GetNeighborhoodRec(matrix_type &A, size_t node, vector<size_t> &onlyN1, vector<bool> &bvisited)
{
	bvisited[node] = true;

	for(typename matrix_type::cRowIterator it = A.beginRow(node); !it.isEnd(); ++it)
	{
		if(bvisited[it.index()] == false && it.value() != 0.0)
		{
			bvisited[it.index()] = true;
			onlyN1.push_back(it.index());
		}
	}
}
template<typename matrix_type>
void GetNeighborhoodRec(matrix_type &A, size_t node, vector<size_t> &onlyN1, vector<size_t> &onlyN2, vector<bool> &bvisited)
{
	bvisited[node] = true;

	GetNeighborhoodRec(A, node, onlyN1, bvisited);
	for(size_t i=0; i<onlyN1.size(); i++)
		GetNeighborhoodRec(A, onlyN1[i], onlyN2, bvisited);
}

template<typename matrix_type>
void GetNeighborhood(matrix_type &A, size_t node, vector<size_t> &onlyN1, vector<size_t> &onlyN2, vector<bool> &bvisited)
{
	onlyN1.clear();
	onlyN2.clear();
	GetNeighborhoodRec(A, node, onlyN1, onlyN2, bvisited);

	// cleanup
	bvisited[node] = false;
	for(size_t i=0; i<onlyN1.size(); i++)
		bvisited[onlyN1[i]] = false;
	for(size_t i=0; i<onlyN2.size(); i++)
		bvisited[onlyN2[i]] = false;
}


template<typename matrix_type>
class FAMGInterpolationCalculator
{
private:
	FAMGInterpolationCalculator(const FAMGInterpolationCalculator<matrix_type> &other);

public:
	FAMGInterpolationCalculator(const matrix_type &A_, double delta, double theta, double damping) : A(A_)
	{
		m_delta = delta;
		m_theta = theta;
		m_damping = damping;
		bvisited.resize(A.num_rows(), false);
	}


	void calculate_testvector(size_t i)
	{
		testvector.resize(onlyN1.size()+1);
		//get_testvector_constant(testvector, onlyN1, testvector[onlyN1.size()], i);
		for(size_t j=0; j<onlyN1.size(); j++)
			testvector[j] = big_testvector[onlyN1[j]];
		testvector[onlyN1.size()] = big_testvector[i];

	}




	// get_possible_parent_pairs:
	//---------------------------------------
	/** calculates of an index i a list of interpolating parent pairs.
	 * \param i						node index for which to calculate possible parent pairs
	 * \param possible_neighbors	list of possible parent pairs
	 * \param nodeinfo				info of the node. possibly set to coarse or fine.
	 */
	void get_possible_parent_pairs(size_t i, std::vector<neighborstruct2> &possible_neighbors, famg_nodeinfo &nodeinfo)
	{
		FAMG_LOG(2, "\n\n\n\n============================\n\n\n");
		FAMG_LOG(2, "node " << GetOriginalIndex(i) << "\n");


		if(get_H(i) == false)
		{
			FAMG_LOG(2, "node has no connections, set as fine\n");
			nodeinfo.set_fine();
			return;
		}

		size_t i_index = onlyN1.size();

		// calculate testvector
		calculate_testvector(i);

		DenseMatrix<FixedArray2<double, 3, 3> > KKT;
		DenseVector<FixedArray1<double, 3> > rhs;
		DenseVector<FixedArray1<double, 3> > q;
		DenseVector<FixedArray1<double, 3> > t;
		rhs[2] = - testvector[i_index];
		i_neighborpairs.clear();

		int i_min = -1;
		double f_min = 1e12;
		for(size_t n=0; n < onlyN1.size(); n++)
		{
			if(A.is_isolated(onlyN1[n]))
				continue;
			for(size_t m=n+1; m < onlyN1.size(); m++)
			{
				if(A.is_isolated(onlyN1[m]))
					continue;
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
				KKT(2, 2) = 0.0;

				rhs[0] = - Hi[n];
				rhs[1] = - Hi[m];

				FAMG_LOG(2, "checking parents " << n << " (" << GetOriginalIndex(onlyN1[n]) << ") and " << m << " (" << GetOriginalIndex(onlyN1[m]) << ")\n");
				IF_FAMG_LOG(3) KKT.maple_print("KKT");
				IF_FAMG_LOG(3) rhs.maple_print("rhs");

				InverseMatMult(q, 1.0, KKT, rhs);

				IF_FAMG_LOG(3) q.maple_print("q");

				neighborstruct2 s;

				IF_FAMG_LOG(3) KKT.maple_print("KKT");
				FAMG_LOG(3, "Hi[n] = " << Hi[n] << ", " << " Hi[m] = " << Hi[m] << " Hi[i_index] = " << Hi[i_index] << "\n");

				// calc q^T H q

				s.F = 	q[0] * (KKT(0,0) * q[0] + KKT(0,1) * q[1] + /* 1.0 */ Hi[n]) +
						q[1] * (KKT(1,0) * q[0] + KKT(1,1) * q[1] + /* 1.0 */ Hi[m]) +
						/*1 */ (Hi[n] * q[0] + Hi[m] * q[1] + Hi[i_index]);

				FAMG_LOG(2, "F: " << s.F << "\n");

				if(s.F > m_delta) continue;
				if(s.F < f_min)
				{
					f_min = s.F;
					i_min = i_neighborpairs.size();
				}
				s.parents[0].from = onlyN1[n];
				s.parents[0].value = -q[0];
				s.parents[1].from = onlyN1[m];
				s.parents[1].value = -q[1];
				i_neighborpairs.push_back(s);

				// F = q^T H q
			}
		}
		possible_neighbors.clear();
		if(i_min != -1)
		{
			// minimal element is always first in possible_neighbors[i] list.
			swap(i_neighborpairs[0], i_neighborpairs[i_min]);

			for(size_t j=0; j<i_neighborpairs.size(); j++)
				if(m_theta*i_neighborpairs[j].F <= f_min)
					possible_neighbors.push_back(i_neighborpairs[j]);

	//#ifdef FAMG_PRINT_POSSIBLE_PARENTS
			IF_FAMG_LOG(2)
			{
				UG_LOG(std::endl << i << ": onlyN1.size = " << onlyN1.size() << " f_min = " << f_min << std::endl);
				UG_LOG(possible_neighbors.size() << " accepted pairs:" << std::endl);
				for(size_t j=0; j<i_neighborpairs.size(); j++)
					if(m_theta*i_neighborpairs[j].F <= f_min)
						i_neighborpairs[j].print();
				UG_LOG("not accepted pairs: " << std::endl);
				for(size_t j=0; j<i_neighborpairs.size(); j++)
					if(m_theta*i_neighborpairs[j].F > f_min)
						i_neighborpairs[j].print();
				UG_LOG(std::endl);
			}
	//#endif
		}
		else
		{
			FAMG_LOG(2, std::endl << GetOriginalIndex(i) << ": UNINTERPOLATEABLE" << std::endl);
			//UG_ASSERT(0, "node has no parents :'-(");
			nodeinfo.set_uninterpolateable();
		}
	}

	// get_all_neighbors_interpolation:
	//---------------------------------------
	/** calculates an interpolation of the node i, when i is interpolating by all his neighbors
	 * if neighbors are fine, their interpolation form coarse nodes is used (indirect)
	 * \param i		node index for which to calculate possible parent pairs
	 * \param P		matrix for the interpolation
	 * \
	 */
	void get_all_neighbors_interpolation(size_t i, matrix_type &P,	famg_nodes &rating, std::vector<int> &newIndex,
					size_t &iNrOfCoarse, size_t &unassigned)
	{
		FAMG_LOG(3, "aggressive coarsening on node " << GetOriginalIndex(i) << "\n")

		get_H(i);

		size_t i_index = onlyN1.size();
		// get testvector
		calculate_testvector(i);

		vector<size_t> coarse_neighbors;
		for(size_t j=0; j<onlyN1.size(); j++)
			if(rating[j].is_coarse())
				coarse_neighbors.push_back(j);

		if(coarse_neighbors.size() >= 1)
		{
			size_t N = coarse_neighbors.size();
			vKKT.resize(N+1, N+1);

			for(size_t r=0; r<N; r++)
				for(size_t c=0; c<N; c++)
					vKKT(r, c) = H(coarse_neighbors[r], coarse_neighbors[c]);

			for(size_t j=0; j < N; j++)
				vKKT(j, N) = vKKT(N, j) = testvector[coarse_neighbors[j]];

			vKKT(N, N) = 0;

			rhs.resize(N+1);
			for(size_t j=0; j < N; j++)
				rhs[j] = -Hi[coarse_neighbors[j]];
			rhs[N] = -testvector[i_index];

			IF_FAMG_LOG(3) vKKT.maple_print("KKT");
			IF_FAMG_LOG(3) rhs.maple_print("rhs");

			InverseMatMult(q, 1.0, vKKT, rhs);

			IF_FAMG_LOG(3) q.maple_print("q");

			t.resize(q.size());
			// todo: calc F
			//MatMult(t, 1.0, H, q);
			double F = 0; //VecDot(q, t);

			if(F > m_delta)
			{
				// set as coarse
				newIndex[i] = iNrOfCoarse++;
				P(i, newIndex[i]) = 1.0;
				rating[i].set_coarse();
				unassigned --;
			}
			else
			{
				for(size_t j=0; j<N; j++)
				{
					size_t node = onlyN1[coarse_neighbors[j]];
					P(i, newIndex[node]) = -q[j];
				}
				rating[i].set_fine();
				unassigned--;
			}
		}
		else
		{

			vKKT.resize(onlyN1.size()+1, onlyN1.size()+1);
			for(size_t r=0; r<onlyN1.size(); r++)
				for(size_t c=0; c<onlyN1.size(); c++)
					vKKT(r, c) = H(r,c);

			 /*  KKT = 	( H     t ) ( q_i,nm )    ( -H e_i )
			  *			( t^T   0 ) ( lambda )  = ( t[i]  ) */

			for(size_t j=0; j < onlyN1.size(); j++)
				vKKT(j, i_index) = vKKT(i_index, j) = testvector[j];
			vKKT(i_index, i_index) = 0;

			rhs.resize(onlyN1.size()+1);
			for(size_t j=0; j < onlyN1.size(); j++)
				rhs[j] = -Hi[j];
			rhs[i_index] = -testvector[i_index];

			IF_FAMG_LOG(3) vKKT.maple_print("KKT");
			IF_FAMG_LOG(3) rhs.maple_print("rhs");


			InverseMatMult(q, 1.0, vKKT, rhs);

			IF_FAMG_LOG(3) q.maple_print("q");

			t.resize(q.size());
			// todo: calc F
			//MatMult(t, 1.0, H, q);
			double F = 0; //VecDot(q, t);

			if(F > m_delta)
			{
				// set as coarse
				newIndex[i] = iNrOfCoarse++;
				P(i, newIndex[i]) = 1.0;
				rating[i].set_coarse();
				unassigned --;
			}
			else
			{
				for(size_t j=0; j<onlyN1.size(); j++)
				{
					size_t node = onlyN1[j];
					for(typename matrix_type::rowIterator it=P.beginRow(node); !it.isEnd(); ++it)
						P(i, it.index()) += -q[j] * it.value();
				}
				rating[i].set_fine();
				unassigned--;


			}
		}






	}


private:
	void calculate_H_from_S()
	{
		size_t i_index = onlyN1.size();
		UG_ASSERT(S.num_cols() == S.num_rows(), "");
		UG_ASSERT(S.num_cols() == onlyN1.size()+1+onlyN2.size(), "");
		size_t N = S.num_rows();

		// get Dinv = 1/Aii
		Dinv.resize(N);
		for(size_t j=0; j < N; j++)
			GetInverse(Dinv[j], S(j,j));

		// get SF = 1-wDF^{-1} A  (F-smoothing)
		SF.resize(N);
		double diaginv = m_damping/S(i_index, i_index);
		for(size_t j=0; j < N; j++)
			SF[j] = - diaginv * S(i_index, j);
		SF[i_index] += 1.0;

		IF_FAMG_LOG(3) print_vector(SF, "SF");

		// get S = 1-w D^{-1} A
		for(size_t r=0; r < N; r++)
		{
			for(size_t c=0; c < N; c++)
			{
				S(r, c) = -m_damping*Dinv[r]*S(r,c);
				if(r==c) S(r, c) += 1.0;
			}
		}

		IF_FAMG_LOG(3) S.maple_print("Sjac");

		// get S = SF S
		// (possible without temporary since SF is mostly Id,
		//  and then SF S is addition of rows of S)
		// r<N1size, c<N1size

		// S'(i_index, i) = sum_t SF(i_index,t) * S(t, i)
		for(size_t i=0; i<N; i++)
			S(i_index, i) *= SF[i_index];
		for(size_t t = 0; t<N; t++)
		{
			if(t==i_index) continue;
			for(size_t i=0; i < N; i++)
				S(i_index, i) += SF[t] * S(t, i);
		}


		/* ?? for(size_t r=0; r<N; r++)
			for(size_t c=0; c < N; c++)
			{
				double s = 0;
				if(c != i_index)
					s += S(r, c); // SF[c,c] = 1.0
				//s += S(r, j) * SF(j, c);
				s += S(r, i_index) * SF[c];
				S(r, c) = s;
			}*/

		IF_FAMG_LOG(3)	S.maple_print("S_SF");

		H.resize(onlyN1.size(),onlyN1.size());
		// get H = S Y S^T
		// todo: use symmetric H.
		for(size_t r=0; r < onlyN1.size(); r++)
			for(size_t c=0; c < onlyN1.size(); c++)
			{
				double s=0;
				for(size_t j=0; j<N; j++)
					s += S(r, j) * Dinv[j] * S(c, j);
				H(r, c) = s;
			}

		IF_FAMG_LOG(2)	H.maple_print("H");

		// get Hi[.] = H(i_index, .)
		Hi.resize(onlyN1.size()+1);
		for(size_t r=0; r < onlyN1.size()+1; r++)
		{
			double s=0;
			for(size_t j=0; j<N; j++)
				s += S(r, j) * Dinv[j] * S(i_index, j);
			Hi[r] = s;
		}

		IF_FAMG_LOG(2) print_vector(Hi, "Hi");
	}
	bool get_H(size_t i)
	{

		GetNeighborhood(A, i, onlyN1, onlyN2, bvisited);

		//IF_FAMG_LOG(2) print_vector(onlyN1, "\nn1-neighbors");
		//IF_FAMG_LOG(2) print_vector(onlyN2, "n2-neighbors");

		if(onlyN1.size() == 0)
			return false;

		N2 = onlyN1;
		N2.push_back(i);
		N2.insert(N2.end(), onlyN2.begin(), onlyN2.end());

		IF_FAMG_LOG(2)
		{
			UG_LOG("\nN1 ");
			for(int i=0; i<onlyN1.size(); i++)
			{
				if(i>0) UG_LOG(", ");
				UG_LOG(GetOriginalIndex(onlyN1[i]));
			}
			UG_LOG("\n");
			//print_vector(onlyN2, "\nN2");
		}

		// get submatrix A on N2

		S.resize(N2.size(), N2.size());
		S = 0.0;
		A.get(S, &N2[0], &N2[0]);

		//IF_FAMG_LOG(3) S.maple_print("\nsubA");

		calculate_H_from_S();

		//IF_FAMG_LOG(3) H.maple_print("\nsubH");

		return true;
	}

private:
	// for speedup purposes, we don't want these arrays to be re-allocated all the time,
	// thats why they are stored here.
	std::vector<bool> bvisited;					// used for N2-neighborhood calculation

	DenseMatrix<VariableArray2<double> > S;		// local matrix S = 1 - wD^{-1}A on {N1, i_index, N2}
	std::vector<double> SF;						// SF = 1 -wD^{-1} A(i_index,.)
	std::vector<double> D;						// diagonal
	std::vector<double> Dinv;						// diagonal
	DenseMatrix<VariableArray2<double> > H;		// matrix H = S Y S^T  on {N1}
	std::vector<double> Hi;						// H(i_index, N1)

	// for the KKT system
	std::vector<double> testvector;
	DenseVector<std::vector<double> > rhs;
	DenseVector<std::vector<double> > q;
	DenseVector<std::vector<double> > t;

	std::vector<neighborstruct2> i_neighborpairs;

	std::vector<size_t> onlyN1;
	std::vector<size_t> onlyN2, N2;

	DenseMatrix<VariableArray2<double> > vKKT;

private:
	const matrix_type &A;
	double m_delta;
	double m_theta;
	double m_damping;

};


#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
