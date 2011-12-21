/**
 * \file amg_debug.h
 *
 * \author Martin Rupp
 *
 * \date 24.11.10
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
#define __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__

#include "famg.h"
#include "famg_nodeinfo.h"
#include "../rsamg/rsamg.h"


namespace ug {


template<typename matrix_type>
void GetNeighborhood(matrix_type &A, size_t node, std::vector<size_t> &onlyN1)
{
	onlyN1.clear();
	for(typename matrix_type::const_row_iterator it = A.begin_row(node); it != A.end_row(node); ++it)
	{
		if(it.value() != 0.0 && it.index() != node)
			onlyN1.push_back(it.index());
	}
}

// adds the indices of the neighbors of a node to the vector onlyN1, if they have visited flag 0.
// node: neither onlyN1 nor bvisited get reset
template<typename matrix_type>
void GetNeighborhoodRec(matrix_type &A, size_t node, std::vector<size_t> &onlyN1, std::vector<bool> &bvisited)
{
	bvisited[node] = true;

	for(typename matrix_type::const_row_iterator it = A.begin_row(node); it != A.end_row(node); ++it)
	{
		if(bvisited[it.index()] == false && it.value() != 0.0)
		{
			bvisited[it.index()] = true;
			onlyN1.push_back(it.index());
		}
	}
}
template<typename matrix_type>
void GetNeighborhoodRec(matrix_type &A, size_t node, std::vector<size_t> &onlyN1, std::vector<size_t> &onlyN2, std::vector<bool> &bvisited)
{
	bvisited[node] = true;

	GetNeighborhoodRec(A, node, onlyN1, bvisited);
	for(size_t i=0; i<onlyN1.size(); i++)
		GetNeighborhoodRec(A, onlyN1[i], onlyN2, bvisited);
}
void CleanupFlagArray(std::vector<bool> &bFlagged, const std::vector<size_t> &iFlagged)
{
	for(size_t i=0; i<iFlagged.size(); i++)
		bFlagged[iFlagged[i]] = false;
}

template<typename matrix_type>
void GetNeighborhood(matrix_type &A, size_t node, std::vector<size_t> &onlyN1, std::vector<size_t> &onlyN2, std::vector<bool> &bvisited)
{
	onlyN1.clear();
	onlyN2.clear();
	GetNeighborhoodRec(A, node, onlyN1, onlyN2, bvisited);

	// cleanup
	bvisited[node] = false;
	CleanupFlagArray(bvisited, onlyN1);
	CleanupFlagArray(bvisited, onlyN2);
}

/**
 * FAMGInterpolationCalculator
 *
 * This class calculates the interpolation weights. It should perhaps be merged with FAMGLevelCalculator
 * its only public member functions are
 * void get_possible_parent_pairs(size_t i, stdvector<neighborstruct2> &possible_neighbors, famg_nodes &rating)
 * void get_all_neighbors_interpolation(size_t i, prolongation_matrix_type &P,	famg_nodes &rating)
 */
template<typename matrix_type, typename vector_type>
class FAMGInterpolationCalculator
{
private:
	FAMGInterpolationCalculator(const FAMGInterpolationCalculator<matrix_type, vector_type> &other);

public:
	template<typename T>
	FAMGInterpolationCalculator(const matrix_type &A_, const matrix_type &A_OL2_,
			const FAMG<T> &famg,
			stdvector< vector_type > &testvectors, stdvector<double> &omega)
	: A(A_), A_OL2(A_OL2_), m_testvectors(testvectors), m_omega(omega)
	{
		m_delta = famg.get_delta();
		m_theta = famg.get_theta();
		m_damping = famg.get_damping_for_smoother_in_interpolation_calculation();
		m_dProlongationTruncation = famg.get_prolongation_truncation();
		m_dHReduceInterpolationNodesParameter = famg.get_H_reduce_interpolation_nodes_parameter();

		testvectorsExtern = (m_testvectors.size() > 0);
	}


	void init()
	{
		bvisited.resize(A_OL2.num_rows(), false);
	}

	template<typename prolongation_matrix_type>
	void check_weights(prolongation_matrix_type &P, size_t i)
	{
		IF_DEBUG(LIB_ALG_AMG, 2)
		{
			for(typename matrix_type::row_iterator it = P.begin_row(i); it != P.end_row(i); ++it)
			{
				if(it.value() < 0.01 || it.value() > 1)
					UG_LOG("P(" << i << ", " << it.index() << ") = " << it.value() << "\n");
			}
		}
	}

	// get_possible_parent_pairs:
	//---------------------------------------
	/** calculates of an index i a list of interpolating parent pairs.
	 * \param i						node index for which to calculate possible parent pairs
	 * \param possible_neighbors	here goes list of possible parent pairs
	 * \param rating				fine/coarse infos of the nodes
	 */
	void get_possible_parent_pairs(size_t i, stdvector<neighborstruct2> &possible_neighbors, FAMGNodes &rating)
	{
		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 2, "\n\n\n\n============================\n\n\n");
		UG_DLOG(LIB_ALG_AMG, 2, "node " << i << " [" << rating.get_original_index(i) << "]\n");

		UG_ASSERT(rating[i].is_valid_rating(), "node " << i << " has no valid rating " << rating[i].rating);

		if(testvectorsExtern && BlockNorm2(m_testvectors[0][i]) < 1e-12 )
		{
			UG_DLOG(LIB_ALG_AMG, 2, "testvector in " << i << " nearly 0, set as fine\n");
			rating.set_fine(i);
			return;
		}

		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);
		if(onlyN1.size() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 2, "node has no connections, set as fine\n");
			rating.set_fine(i);
			return;
		}

		size_t i_index = onlyN1.size();

		// get H on onlyN1,i_index,onlyN2.
		get_H(i, rating);

		// calculate testvector
		calculate_testvectors(i);

		DenseMatrix<FixedArray2<double, 3, 3> > KKT;
		DenseVector<FixedArray1<double, 3> > rhs;
		DenseVector<FixedArray1<double, 3> > q;
		DenseVector<FixedArray1<double, 3> > t;
		double t_i = localTestvector[0][i_index];
		rhs[2] = -t_i;
		i_neighborpairs.clear();


		AMG_PROFILE_NEXT(FAMG_getPPP_minCalc);
		int i_min = -1;
		double f_min = 1e12;

	
		IF_DEBUG(LIB_ALG_AMG, 5)
		{
			for(size_t n=0; n < onlyN1.size(); n++)
				UG_DLOG(LIB_ALG_AMG, 4, onlyN1[n] << " = " << rating.info(onlyN1[n]) << " ");
			UG_DLOG(LIB_ALG_AMG, 4, "\n");
		}

		neighborstruct2 s;

		const double &aii = A_OL2(i,i);
		for(size_t n=0; n < onlyN1.size(); n++)
		{
			if(!rating.i_can_set_coarse(onlyN1[n])) // || A.is_isolated(onlyN1[n]))
				continue;
			for(size_t m=n+1; m < onlyN1.size(); m++)
			{
				if(!rating.i_can_set_coarse(onlyN1[m])) // || A.is_isolated(onlyN1[m]))
					continue;
				// set KKT matrix
				/*
				 *  KKT = 	( H     t ) ( q_i,nm )    ( -H e_i )
				 * 			( t^T   0 ) ( lambda )  = ( t[i]  )
				 *
				 *  H *q = H(n,n) * q_n + H(n, m) * q_m + H(n, i) * q_i
				 *         H(m,n) * q_n + H(m, m) * q_m + H(m, i) * q_i
				 */

#if 1
				KKT(0, 0) = H(n, n);
				KKT(0, 1) = H(n, m);
				KKT(1, 0) = H(m, n);
				KKT(1, 1) = H(m, m);

				KKT(2, 0) = KKT(0, 2) = localTestvector[0][n];
				KKT(2, 1) = KKT(1, 2) = localTestvector[0][m];
				KKT(2, 2) = 0.0;

				rhs[0] = - H(i_index, n);
				rhs[1] = - H(i_index, m);

				UG_DLOG(LIB_ALG_AMG, 2, "checking parents " << n << " (" << rating.get_original_index(onlyN1[n]) << ") and " << m << " (" << rating.get_original_index(onlyN1[m]) << ")\n");
				IF_DEBUG(LIB_ALG_AMG, 5) KKT.maple_print("KKT");
				IF_DEBUG(LIB_ALG_AMG, 5) rhs.maple_print("rhs");

				if(InverseMatMult(q, 1.0, KKT, rhs) == false)
				{
					UG_DLOG(LIB_ALG_AMG, 3, "get_possible_parent_pairs: could not invert KKT system.\n");
					continue;
				}

				IF_DEBUG(LIB_ALG_AMG, 5) q.maple_print("q");

				IF_DEBUG(LIB_ALG_AMG, 5) KKT.maple_print("KKT");
				UG_DLOG(LIB_ALG_AMG, 5, "H(i, n) = " << H(i_index, n) << ", " << " H(i, m) = " << H(i_index, m) << " H(i, i) = " << H(i_index, i_index) << "\n");

				if(q[0] > 0 || q[1] > 0) continue;

				// calc q^T H q

				s.F = 	q[0] * (KKT(0,0) * q[0] + KKT(0,1) * q[1] + /* 1.0 */ H(i_index, n)) +
						q[1] * (KKT(1,0) * q[0] + KKT(1,1) * q[1] + /* 1.0 */ H(i_index, m)) +
						/*1 */ (H(i_index, n) * q[0] + H(i_index, m) * q[1] + H(i_index, i_index));
				// diagonal scaling is made here:
				s.F *= aii;
				UG_DLOG(LIB_ALG_AMG, 2, "F: " << s.F << "\n");

				if(s.F > m_delta) continue;


				s.parents[0].from = onlyN1[n];
				s.parents[0].value = -q[0];
				s.parents[1].from = onlyN1[m];
				s.parents[1].value = -q[1];
				if(s.F < f_min)
				{
					f_min =s.F;
					i_min = i_neighborpairs.size();
				}

				//std::cout << s.parents[0].value << "   " << s.parents[1].value << " : " << s.F << "\n";

#else
				double qn, qm;
				UG_DLOG(LIB_ALG_AMG, 5, "H(i, n) = " << H(i_index, n) << ", " << " H(i, m) = " << H(i_index, m) << " H(i, i) = " << H(i_index, i_index) << "\n");

				double tn = localTestvector[0][n];
				double tm = localTestvector[0][m];
				double nenner1 = tn*tn*H(m,m) - 2.0*tn*tm*H(n,m) + tm*tm*H(n,n);
				if(nenner1 < 1e-12) continue;
				 // construct prolongation
				if(abs(tm) > abs(tn))
				{
					//std::cout << "a\n";
					// Hnn = Hnn
					// Hm0 Hnn
					double zaehler1 = -(tm*tm * H(i_index, n) - tm*tn * H(i_index, m) + (tn*H(m,m) - tm*H(n,m)) * t_i);
					double laenge1 = H(i_index, i_index)*tm*tm - 2.0*t_i*tm*H(i_index, m) + t_i*t_i*H(m,m);

					s.F = (laenge1 - zaehler1*zaehler1/nenner1)/(tm*tm);

					s.F *= aii; // diagonal scaling is made here
					if(s.F > m_delta) continue;

					qn = zaehler1/nenner1;
					qm = (-t_i - qn*tn)/tm;
				}

				else
				{
					//std::cout << "b\n";
					double zaehler1 =  -( - tn*tm * H(i_index, n) + tn*tn * H(i_index, m) + (tm*H(n,n) - tn*H(n,m)) * t_i);
					double laenge1 = H(i_index, i_index)*tn*tn - 2.0*t_i*tn*H(i_index, n) + t_i*t_i*H(n,n);
					s.F = (laenge1 - zaehler1*zaehler1/nenner1)/(tn*tn);

					s.F *= aii; // diagonal scaling is made here
					if(s.F > m_delta) continue;

					qm = zaehler1/nenner1;
					qn = (-t_i - qm*tm)/tn;
			   }


				if(qm > 0 || qn > 0) continue;

				IF_DEBUG(LIB_ALG_AMG, 5) q.maple_print("q");
				s.parents[0].from = onlyN1[n];
				s.parents[0].value = -qn;
				s.parents[1].from = onlyN1[m];
				s.parents[1].value = -qm;
				//std::cout << s.parents[0].value << "   " << s.parents[1].value << " : " << s.F << "\n\n";
				// calc q^T H q

				UG_DLOG(LIB_ALG_AMG, 2, "F: " << s.F << "\n");

				if(s.F < f_min)
				{
					f_min =s.F;
					i_min = i_neighborpairs.size();
				}
#endif
				i_neighborpairs.push_back(s);

				// F = q^T H q
			}
		}
		AMG_PROFILE_NEXT(FAMG_getPPP_cleanup);
		possible_neighbors.clear();
		if(i_min != -1)
		{
			// minimal element is always first in possible_neighbors[i] list.
			std::swap(i_neighborpairs[0], i_neighborpairs[i_min]);

			for(size_t j=0; j<i_neighborpairs.size(); j++)
				if(m_theta*i_neighborpairs[j].F <= f_min)
					possible_neighbors.push_back(i_neighborpairs[j]);

	//#ifdef FAMG_PRINT_POSSIBLE_PARENTS
			IF_DEBUG(LIB_ALG_AMG, 2)
			{
				UG_LOG(std::endl << i << ": onlyN1.size = " << onlyN1.size() << " f_min = " << f_min << std::endl);
				UG_LOG(possible_neighbors.size() << " accepted pairs:" << std::endl);
				for(size_t j=0; j<i_neighborpairs.size(); j++)
					if(m_theta*i_neighborpairs[j].F <= f_min)
						i_neighborpairs[j].print(rating);
				UG_LOG("not accepted pairs: " << std::endl);
				for(size_t j=0; j<i_neighborpairs.size(); j++)
					if(m_theta*i_neighborpairs[j].F > f_min)
						i_neighborpairs[j].print(rating);
				UG_LOG(std::endl);
			}
	//#endif
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 2, std::endl << rating.get_original_index(i) << ": UNINTERPOLATEABLE" << std::endl);
			//UG_ASSERT(0, "node has no parents :'-(");
			UG_ASSERT(rating.i_must_assign(i), i);
			rating.set_uninterpolateable(i);
		}
	}

	inline bool solve_KKT()
	{
		AMG_PROFILE_FUNC();
		size_t N = onlyN1.size();
		size_t i_index = N;
		vKKT.resize(N+1, N+1);

		for(size_t r=0; r<N; r++)
			for(size_t c=0; c<N; c++)
				vKKT(r, c) = H(r, c);

		for(size_t j=0; j < N; j++)
			vKKT(j, N) = vKKT(N, j) = localTestvector[0][j];

		vKKT(N, N) = 0;

		rhs.resize(N+1);
		for(size_t j=0; j < N; j++)
			rhs[j] = -H(i_index, j);
		rhs[N] = -localTestvector[0][i_index];

		IF_DEBUG(LIB_ALG_AMG, 5) vKKT.maple_print("KKT");
		IF_DEBUG(LIB_ALG_AMG, 5) rhs.maple_print("rhs");

		q.resize(N+1);
		if(InverseMatMult(q, 1.0, vKKT, rhs))
		{
			F = H(i_index, i_index) -localTestvector[0][i_index];
			for(size_t j=0; j<N; j++)
				F += q[j] * rhs[j];
			return true;
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 1, "solve_KKT: could not invert KKT :" << vKKT);
			return false;

		}
	}

	inline bool solve_KKT_on_subset(const std::vector<size_t> &indices, size_t i)
	{
		AMG_PROFILE_FUNC();
		size_t i_index = onlyN1.size();
		size_t N = indices.size();
		vKKT.resize(N+1, N+1);

		for(size_t r=0; r<N; r++)
			for(size_t c=0; c<N; c++)
				vKKT(r, c) = H(indices[r], indices[c]);

		for(size_t j=0; j < N; j++)
			vKKT(j, N) = vKKT(N, j) = localTestvector[0][indices[j]];


		vKKT(N, N) = 0;

		rhs.resize(N+1);
		for(size_t j=0; j < N; j++)
			rhs[j] = -H(i_index, indices[j]);
		rhs[N] = -localTestvector[0][i_index];

		IF_DEBUG(LIB_ALG_AMG, 5) vKKT.maple_print("KKT");
		IF_DEBUG(LIB_ALG_AMG, 5) rhs.maple_print("rhs");

		q.resize(N+1);

		if(InverseMatMult(q, 1.0, vKKT, rhs))
		{
			F = H(i_index, i_index) -localTestvector[0][i_index];
			for(size_t j=0; j<N; j++)
				F += q[j] * rhs[j];
			return true;
		}
		else
		{
			UG_DLOG(LIB_ALG_AMG, 1, "solve_KKT_on_subset: could not invert KKT in "
				 << i << ":\n" << vKKT);
			return false;
		}
	}

	void printH(size_t i)
	{
		std::vector<size_t> N;
		for(size_t j=0; j < onlyN1.size(); j++)
			N.push_back(onlyN1[j]);
		N.push_back(i);

		for(size_t r=0; r < H.num_rows(); r++)
		{
			UG_LOG(N[r] <<": ");
			for(size_t c=0; c < H.num_cols();c ++)
				UG_LOG("(" << N[c] << "): " << H(r, c) << " ");
			UG_LOG("\n");
		}
	}

	bool get_N2(size_t i, std::vector<size_t> &myN1)
	{
		AMG_PROFILE_FUNC();
		swap(onlyN1, myN1);

		for(size_t j=0; j<bvisited.size(); j++)
			bvisited[j] = false;
		bvisited[i] = true;
		for(size_t j=0; j<onlyN1.size(); j++)
			bvisited[onlyN1[j]] = true;

		onlyN2.clear();
		for(size_t j=0; j<onlyN1.size(); j++)
		{
			UG_LOG("Neighbhors of " << onlyN1[j] << ": ");
			size_t k=onlyN2.size();
			GetNeighborhoodRec(A_OL2, onlyN1[j], onlyN2, bvisited);
			for(; k<onlyN2.size(); k++)
				UG_LOG(onlyN2[k] << " ");
		}

		CleanupFlagArray(bvisited, onlyN1);
		CleanupFlagArray(bvisited, onlyN2);
		bvisited[i] = false;
		return true;
	}

	void reduceToStrong(const std::vector<size_t> &coarseNeighbors, std::vector<size_t> &strongNeighbors, double eps)
	{
		size_t i_index = onlyN1.size();
		double maxVal = H(i_index, coarseNeighbors[0]);
		for(size_t j=1; j<coarseNeighbors.size(); j++)
			maxVal = std::max(H(i_index, coarseNeighbors[j]), maxVal);
		for(size_t j=0; j<coarseNeighbors.size(); j++)
			if(H(i_index, coarseNeighbors[j]) > eps*maxVal)
				strongNeighbors.push_back(coarseNeighbors[j]);
	}

	template<typename prolongation_matrix_type>
	bool indirect_interpolation(size_t i, prolongation_matrix_type &P,	FAMGNodes &rating,
			std::vector<size_t> &coarseNeighbors)
	{
		AMG_PROFILE_FUNC();
		if(coarseNeighbors.size() == 1)
		{
			get_H(i, rating);
			calculate_testvectors(i);
			P(i, onlyN1[coarseNeighbors[0]]) = localTestvector[0][onlyN1.size()] / localTestvector[0][coarseNeighbors[0]];
			rating.set_fine(i);
			return true;
		}

		get_H(i, rating);
		calculate_testvectors(i);

		std::vector<size_t> strongNeighbors;
		//reduceToStrong(coarseNeighbors, strongNeighbors, 0.5);
		reduceToStrong(coarseNeighbors, strongNeighbors, 0.1);

		if(solve_KKT_on_subset(strongNeighbors, i))
		{
			IF_DEBUG(LIB_ALG_AMG, 5) q.maple_print("q");

			t.resize(q.size());

			if(F > m_delta)
			{
				UG_DLOG(LIB_ALG_AMG, 3, "coarse neighbors, had to set node " << i << " coarse!");
				rating.external_set_coarse(i);
			}
			else
			{
				UG_DLOG(LIB_ALG_AMG, 3, "coarse neighbors, Interpolating from ");
				size_t N = strongNeighbors.size();
				if(m_dProlongationTruncation > 0.0)
				{
					double dmax = -q[0];
					for(size_t j=1; j<N; j++)
						if(dmax < -q[j]) dmax = -q[j];
					for(size_t j=0; j<N; j++)
					{
						if(-q[j] < dmax*m_dProlongationTruncation) continue;
						size_t node = onlyN1[strongNeighbors[j]];
						UG_DLOG(LIB_ALG_AMG, 3, node << ": " << q[j] << ", ");
						P(i, node) = -q[j];
					}
				}
				else
				{
					for(size_t j=0; j<N; j++)
					{
						size_t node = onlyN1[strongNeighbors[j]];
						P(i, node) = -q[j];
						UG_DLOG(LIB_ALG_AMG, 3, node << ": " << q[j] << ", ");
					}
				}
				check_weights(P, i);
				rating.set_fine(i);
			}
		}
		else
		{
			rating.external_set_coarse(i);
			UG_DLOG(LIB_ALG_AMG, 3, "indirect_interpolation of node " << i << ": could not invert KKT system (coarse neighbors).\n");
		}
		return true;
	}

	template<typename prolongation_matrix_type>
	bool direct_interpolation(size_t i, prolongation_matrix_type &P, FAMGNodes &rating,
			std::vector<size_t> &interpolateNeighbors)
	{
		AMG_PROFILE_FUNC();
		get_H(i, rating);
		calculate_testvectors(i);
		F=0;
		if(solve_KKT_on_subset(interpolateNeighbors, i) && F < m_delta)
		{
			IF_DEBUG(LIB_ALG_AMG, 5) q.maple_print("q");
			std::map<size_t, double> localP;
			for(size_t j=0; j<interpolateNeighbors.size(); j++)
			{
				size_t node = onlyN1[interpolateNeighbors[j]];
				for(typename matrix_type::row_iterator it=P.begin_row(node); it != P.end_row(node); ++it)
					localP[it.index()] += -q[j] * it.value();
			}
			double dmax = -1000;
			for(std::map<size_t, double>::iterator it = localP.begin(); it != localP.end(); ++it)
				if(dmax < (*it).second) dmax = (*it).second;

			for(std::map<size_t, double>::iterator it = localP.begin(); it != localP.end(); ++it)
			{
				if((*it).second > dmax*m_dProlongationTruncation)
					P(i, (*it).first) = (*it).second;
			}

			check_weights(P, i);
			rating.set_aggressive_fine(i);
		}
		else
		{
			rating.external_set_coarse(i);
			UG_DLOG(LIB_ALG_AMG, 3, "had to set node " << i << " coarse: ");
			if(F==0)
			{	UG_DLOG(LIB_ALG_AMG, 3, "could not invert KKT system.\n"); }
			else
			{	UG_DLOG(LIB_ALG_AMG, 3, "F = " << F << " < delta = " << m_delta << "\n"); }
		}
		return true;
	}



	void get_strong_connected_neighborhood(const matrix_type &A, const FAMGNodes &rating, size_t i,
			stdvector<size_t> &onlyN1, stdvector<size_t> &onlyN2, stdvector<size_t> &coarseNeighbors,
			stdvector<bool> &bvisited)
	{
		double minConnValue, maxConnValue, diag;
		GetNeighborValues(A, i, minConnValue, maxConnValue, diag);
		onlyN1.clear();
		onlyN2.clear();

		bvisited[i] = true;
		const double theta = 0.3;
		double barrier = theta*std::max(maxConnValue, dabs(minConnValue));
		for(typename matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
		{
			if(dabs(amg_offdiag_value(conn.value())) < barrier) // only strong
				continue;
			if(conn.index() == i) continue; // skip diagonal
			onlyN1.push_back(conn.index());
			bvisited[conn.index()] = true;
			if(!rating[conn.index()].is_coarse()) continue; // only coarse

			coarseNeighbors.push_back(onlyN1.size()-1);
		}

		for(size_t i=0; i<onlyN1.size(); i++)
		{
			GetNeighborValues(A, onlyN1[i], minConnValue, maxConnValue, diag);
			barrier = theta*std::max(maxConnValue, dabs(minConnValue));
			for(typename matrix_type::const_row_iterator conn = A.begin_row(onlyN1[i]);
					conn != A.end_row(onlyN1[i]); ++conn)
			{
				if(dabs(amg_offdiag_value(conn.value())) < barrier) // only strong
					continue;
				if(bvisited[conn.index()] == false)
				{
					onlyN2.push_back(conn.index());
					bvisited[conn.index()] = true;
				}
			}
			// GetNeighborhoodRec(A, onlyN1[i], onlyN2, bvisited);
		}

		CleanupFlagArray(bvisited, onlyN1);
		CleanupFlagArray(bvisited, onlyN2);
	}

	void get_strong_connected_neighborhood2(const matrix_type &A, const FAMGNodes &rating, size_t i, stdvector<size_t> &onlyN1,
				stdvector<size_t> &onlyN2, stdvector<size_t> &coarseNeighbors, stdvector<bool> &bvisited)
	{
		GetNeighborhood(A_OL2, i, onlyN1, onlyN2, bvisited);

		for(size_t j=0; j<onlyN1.size(); j++)
		{
			if(rating[onlyN1[j]].is_coarse())
				coarseNeighbors.push_back(j);
		}
	}

	// get_all_neighbors_interpolation:
	//---------------------------------------
	/** calculates an interpolation of the node i, when i is interpolating by all his neighbors
	 * if neighbors are fine, their interpolation form coarse nodes is used (indirect)
	 * \param i					node index for which to calculate possible parent pairs
	 * \param P					matrix for the interpolation
	 * \param rating			fine/coarse infos of the nodes
	 * \param indirectInterpolationLevel	if true, only use coarse neighbors
	 */
	template<typename prolongation_matrix_type>
	bool get_all_neighbors_interpolation(size_t i, prolongation_matrix_type &P,	FAMGNodes &rating,
			int indirectInterpolationLevel=-1)
	{
		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 3, "aggressive coarsening on node " << rating.get_original_index(i) << "\n")

		stdvector<size_t> coarseNeighbors;

		get_strong_connected_neighborhood2(A_OL2, rating, i, onlyN1, onlyN2, coarseNeighbors, bvisited);
		if(onlyN1.size() == 0)
		{
			UG_DLOG(LIB_ALG_AMG, 3, "no neighbors -> gets fine");
			rating.set_fine(i);
			return true;
		}

		//const double &aii = A_OL2(i,i);


		if(coarseNeighbors.size() >= 1)
			return indirect_interpolation(i, P, rating, coarseNeighbors);
		else if(indirectInterpolationLevel >= 0)
		{
			UG_DLOG(LIB_ALG_AMG, 3, "indirectInterpolationLevel = " << indirectInterpolationLevel << ", but no coarse neighbors of this level");
			return false;
		}
		else
		{
			std::vector<size_t> innerNodes;
			for(size_t j=0; j<onlyN1.size(); j++)
				if(rating.is_inner_node(onlyN1[j]) || rating.is_master(onlyN1[j]))
					innerNodes.push_back(j);
			return direct_interpolation(i, P, rating, innerNodes);
		}
	}


private:
	/**
	 * get_H
	 * \param i 		central node i
	 * \param rating	only used for debug output (get_orinal_index)
	 */
	bool get_H(size_t i, FAMGNodes &rating)
	{
		AMG_PROFILE_FUNC();
		// replace this with
		// stdvector<stdvector<size_t> > neighbors(3);
		/// stdvector<size_t> &onlyN1 = neighbors[1];
		/// stdvector<size_t> &onlyN2 = neighbors[2];
		// GetNeighborhoodHierachy(A, i, neighbors, bvisited);

		// 1. Get Neighborhood N1 and N2 of i.
		//AMG_PROFILE_BEGIN(AMG_H_GetNeighborhood);


		IF_DEBUG(LIB_ALG_AMG, 5) print_vector(onlyN1, "\nn1-neighbors");
		IF_DEBUG(LIB_ALG_AMG, 5) print_vector(onlyN2, "n2-neighbors");

		//AMG_PROFILE_NEXT(AMG_H_CreateN2);
		N2 = onlyN1;
		N2.push_back(i);
		N2.insert(N2.end(), onlyN2.begin(), onlyN2.end());

		IF_DEBUG(LIB_ALG_AMG, 2)
		{
			print_vector(onlyN1, "\nonlyN1");
			print_vector(onlyN2, "\nonlyN2");
			print_vector(N2, "\nN2");
		}

		// 2. get submatrix in A_OL2 on N2
		//AMG_PROFILE_NEXT(AMG_H_GetLocalMatrix);

		S.resize(N2.size(), N2.size());
		S = 0.0;

		GetLocalMatrix(A_OL2, S, &N2[0], &N2[0]);


		// make S symmetric
		/*for(size_t r=0; r<S.num_rows(); r++)
			for(size_t c=0; c<S.num_cols(); c++)
			{
				S(r, c) = (S(r,c)+S(c,r))/2
			}*/

		IF_DEBUG(LIB_ALG_AMG, 5) S.maple_print("\nsubA");

		//AMG_PROFILE_END();
		// 3. calculate H from submatrix A
		calculate_H_from_local_A(i);

		IF_DEBUG(LIB_ALG_AMG, 5) H.maple_print("\nsubH");

		return true;
	}

	/**
	 * calculate_H_from_local_A
	 * calculates locally
	 * - F-Smoothing SF = 1-w DF^{-1} A  (on only1, i, only2)
	 * - S = SF (1-w D^{-1} A)			 (on only1, i, only2)
	 * - H = H = S Y S^T				 (on only1)
	 * - Hi[.] = H(., i)				 (on only1)
	 */
	void calculate_H_from_local_A(size_t i)
	{
		AMG_PROFILE_FUNC();
		size_t i_index = onlyN1.size();
		UG_ASSERT(S.num_cols() == S.num_rows(), "");
		UG_ASSERT(S.num_cols() == onlyN1.size()+1+onlyN2.size(), "");
		size_t N = S.num_rows();

		if(Dinv.size() != N) Dinv.resize(N);
		if(SF.size() != N) SF.resize(N);
		//AMG_PROFILE_NEXT(AMG_HA_calculate_H);
		if(H.num_rows() != onlyN1.size()+1)
			H.resize(onlyN1.size()+1, onlyN1.size()+1);

		//AMG_PROFILE_NEXT(AMG_HA_Dinv);
		// get Dinv = 1/Aii
		
		IF_DEBUG(LIB_ALG_AMG, 1)
		{
			for(size_t j=0; j<N; j++)
				if(dabs(S(j,j)) < 1e-12)
					UG_LOG(i << " j=" << j << " S(j,j)=" << S(j,j) << "\n");
		}

		for(size_t j=0; j < N; j++)
			GetInverse(Dinv[j], S(j,j));
		
		// get SF = 1-wDF^{-1} A  (F-smoothing)
		//AMG_PROFILE_BEGIN(AMG_HA_calculate_SF);

		// bei f-smoothing nie und nimmer damping (arne 3.juni)
		if(dabs(S(i_index, i_index)) < 1e-12)
			UG_LOG("i=" << i << "S(i,i) = " << S(i,i) << "\n");
		for(size_t j=0; j < N; j++)
			SF[j] = - Dinv[i_index] * S(i_index, j);
		SF[i_index] += 1.0;


		IF_DEBUG(LIB_ALG_AMG, 5) print_vector(SF, "SF");

		// get S = 1-w D^{-1} A
		for(size_t r=0; r < N; r++)
		{
			for(size_t c=0; c < N; c++)
			{
				S(r, c) = -m_damping*Dinv[r]*S(r,c);
				if(r==c) S(r, c) += 1.0;
			}
		}

		IF_DEBUG(LIB_ALG_AMG, 5) S.maple_print("Sjac");

		// get S' = SF S
		//AMG_PROFILE_NEXT(AMG_HA_calculate_SFS);
		// (possible without temporary since SF is mostly Id,
		//  and then SF S is addition of rows of S)
		// r<N1size, c<N1size

		// S'(i_index, i) = sum_t SF(i_index,t) * S(t, i)
		for(size_t j=0; j<N; j++)
			S(i_index, j) *= SF[i_index];
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

		IF_DEBUG(LIB_ALG_AMG, 5)	S.maple_print("S_SF");


		// get H = S' Y S'^T
		// todo: use symmetric H.
		for(size_t r=0; r < onlyN1.size()+1; r++)
			for(size_t c=r; c < onlyN1.size()+1; c++)
			{
				double s=0;
				for(size_t j=0; j<N; j++)
					s += S(r, j) * Dinv[j] * S(c, j);
				H(r, c) = s;
				H(c, r) = s;
			}

		IF_DEBUG(LIB_ALG_AMG, 2)	H.maple_print("H");
	}

	// global_to_local_testvectors
	/** calculates the local testvectors
	 * \param node the fine node to be interpolated
	 *	note: das mit den testvektoren ist noch nicht so sicher:
	 *	sie mŸssen ja theoretisch 2 mal geglŠttet werden. 1x normal jacobi,
	 *	und dann ein 2. mal nur auf den feinen knoten. und danach muss noch
	 *	$\frac 1 {\abs{t}_A}$ geteilt werden. Leider wei§ man aber ja vorher noch nicht,
	 *	welche Knoten fein sind. Dh. im Moment ist das so: Man macht einmal Jacobi,
	 *	berechnet dann mal $\frac 1 {\abs{t}_A}$, und macht dann lokal nur noch so eine
	 *	NachglŠttung, wenn die globalen Testvektoren in lokale Testvektoren umgerechnet werden.
	 */
	void global_to_local_testvectors(size_t node)
	{
		AMG_PROFILE_FUNC();
		localTestvector.resize(m_testvectors.size());
		for(size_t k=0; k<m_testvectors.size(); k++)
		{
			localTestvector[k].resize(onlyN1.size()+1);

			for(size_t j=0; j<onlyN1.size(); j++)
				localTestvector[k][j] = m_testvectors[k][onlyN1[j]];

			localTestvector[k][onlyN1.size()] = m_testvectors[k][node];
		}

		// smooth localTestvectors:
		for(size_t k=0; k<m_testvectors.size(); k++)
		{
			double s=0;
			for(size_t j=0; j < onlyN1.size()+1; j++)
				s += localTestvector[k][j] * SF[j];

			localTestvector[k][onlyN1.size()] = s;
			IF_DEBUG(LIB_ALG_AMG, 5) print_vector(localTestvector[k], "\nlocalTestvector");
		}

	}

	void calculate_EV_testvectors(size_t node)
	{
		std::vector<std::vector<size_t> > neighbors(4);
		std::vector<size_t> &onlyN1 = neighbors[1];
		std::vector<size_t> &onlyN2 = neighbors[2];
		std::vector<size_t> &onlyN3 = neighbors[3];
		GetNeighborhoodHierachy(A_OL2, node, 3, neighbors);

		// print_vector(onlyN1, "\nn1-neighbors");
		// print_vector(onlyN2, "\nn2-neighbors");
		// print_vector(onlyN3, "\nn3-neighbors");

		stdvector<size_t> myNeigh;
		myNeigh = onlyN1;
		myNeigh.push_back(node);
		myNeigh.insert(myNeigh.end(), onlyN2.begin(), onlyN2.end());
		myNeigh.insert(myNeigh.end(), onlyN3.begin(), onlyN3.end());

		DenseMatrix<VariableArray2<double> > localA;
		localA.resize(myNeigh.size(), myNeigh.size());

		GetLocalMatrix(A_OL2, localA, &myNeigh[0], &myNeigh[0]);

		localA.maple_print("A");

		size_t N1 = onlyN1.size();
		size_t N2 = onlyN2.size();
		DenseMatrix<VariableArray2<double> > localS;
		localS.resize(N1+1+N2, N1+1+N2);


		for(size_t r=0; r<N1+1+N2; r++)
		{
			for(size_t c=0; c < N1+1+N2; c++)
				localS(r, c) = localA(r, c);
			for(size_t c=N1+1+N2; c < localA.num_cols(); c++)
				localS(r, r) += localA(r, c);
		}

		 localS.maple_print("S1");


		// get S = 1-D^{-1} A
		for(size_t r=0; r < N1+1+N2; r++)
		{
			double diagInv = 1/localS(r,r);
			for(size_t c=0; c < N1+1+N2; c++)
			{
				localS(r, c) = -diagInv*localS(r,c);
				if(r==c) localS(r, c) += 1.0;
				else if(localS(r, c) < 1e-12)
					localS(r, c) = 0.0;
			}
		}

		localS.maple_print("S2");

		size_t N = N1+1+N2;
		DenseMatrix<VariableArray2<double> > X;
		X.resize(N, N);
		DenseMatrix<VariableArray2<double> > B;
		B.resize(N, N);
		B = 1.0;
		DenseVector<VariableArray1<double> > lambda;
		lambda.resize(N);

		/*int res =*/
		GeneralizedEigenvalueProblem(localS, X, lambda, B, true);


		// X.maple_print("X");
		// lambda.maple_print("lambda");

		localTestvector.resize(1);
		localTestvector[0].resize(onlyN1.size()+1);
		for(size_t i=0; i<onlyN1.size()+1; i++)
			localTestvector[0][i] = X(i, N-1);
	}

	// calculate_testvectors
	//---------------------------------------
	void calculate_testvectors(size_t node)
	{
		AMG_PROFILE_FUNC();

		if(testvectorsExtern == true)
		{
			global_to_local_testvectors(node);
			add_additional_testvectors_to_H();
		}
		else
			calculate_EV_testvectors(node);
	}

	/** add_additional_testvectors_to_H
	 * adds additional factors from the testvectors to H, namely
	 * H += \omega_k t^{(k)} t^{(k)}^T
	 * for the testvectors 1 .. m_testvectors.size()-1
	 * (testvector 0 is approximated exactly)
	 */
	void add_additional_testvectors_to_H()
	{
		AMG_PROFILE_FUNC();
		// skip first vector (it is approximated exactly)
		// /test if other way round is faster/
		for(size_t r = 0; r < onlyN1.size(); r++)
			for(size_t c = 0; c < onlyN1.size(); c++)
			{
				double s = 0;
				for(size_t k=1; k<m_testvectors.size(); k++)
					s+= m_omega[k] * localTestvector[k][r] * localTestvector[k][c];
				H(r, c) += s;
			}
	}
private:
	// for speedup purposes, we don't want these arrays to be re-allocated all the time,
	// thats why they are stored here.
	stdvector<bool> bvisited;					// used for N2-neighborhood calculation

	// todo: instead of VariableArray2, use ReserveableArray2
	// onlyN1: 1-neighborhood without i
	// onlyN2: 2-neighborhood without i and N1.
	DenseMatrix<VariableArray2<double> > S;		//< local matrix S = 1 - wD^{-1}A on {onlyN1, i, onlyN2}
	stdvector<double> SF;						//< SF = 1 -wD^{-1} A(ix,.)  on {onlyN1, i, onlyN2}
	stdvector<double> D;						//< diagonal
	stdvector<double> Dinv;						//< diagonal on {onlyN1, i, onlyN2}
	DenseMatrix<VariableArray2<double> > H;		//< matrix H = S Y S^T  on {onlyN1}

	// for the KKT system
	stdvector<stdvector<double> > localTestvector;			//< on {onlyN1, i}
	DenseVector<stdvector<double> > rhs;
	DenseVector<stdvector<double> > q;
	DenseVector<stdvector<double> > t;

	stdvector<neighborstruct2> i_neighborpairs;

	stdvector<size_t> onlyN1;
	stdvector<size_t> onlyN2, N2;

	DenseMatrix<VariableArray2<double> > vKKT;

private:
	double F;
	const matrix_type &A;
	const matrix_type &A_OL2;
	double m_delta;
	double m_theta;
	double m_damping;
	double m_dProlongationTruncation;
	double m_dHReduceInterpolationNodesParameter;
	stdvector< vector_type > &m_testvectors;
	stdvector<double> &m_omega;

	bool testvectorsExtern;
};

} // namespace ug
#endif // __H__LIB_ALGEBRA__FAMG_SOLVER__FAMG_INTERPOLATION_CALCULATOR_H__
