/**
 * \file amg_rs_prolongation.h
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_RS_PROLONGATION_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_RS_PROLONGATION_H__

#include <vector>
#include "amg_nodeinfo.h"

namespace ug {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateProlongation:
//-------------------------
/**
 * Calculates Prolongation P with Matrix_type A and coarse/fine markers in
 * nodes[i].is_/is_coarse by direct interpolation
 * \param 	P				Matrix P: here goes the calculated prolongation
 * \param	A				Matrix A: matrix for which to calculate prolongation on next level
 * \param	nodes			fine/coarse marks of the nodes.
 * \param	theta			\f$\epsilon_{str}\f$.
 */
template<typename Matrix_type>
void CreateRugeStuebenProlongation(SparseMatrix<double> &P, const Matrix_type &A, AMGNodes &nodes,
		stdvector<int> &newIndex, double theta)
{
	UG_ASSERT(newIndex.size() == nodes.size(), "");
	P.resize(A.num_rows(), nodes.get_nr_of_coarse());

	std::vector<SparseMatrix<double>::connection> con(255);
	SparseMatrix<double>::connection c;

	// DIRECT INTERPOLATION

	for(size_t i=0; i < A.num_rows(); i++)
	{
		if(nodes[i].is_coarse())
		{
			// a coarse node
			UG_ASSERT(newIndex[i] != -1, "coarse node but no new index?");
			P(i, newIndex[i]) = 1.0;
		}
		else if(A.is_isolated(i))
		{
			//P[i].initWithoutDiag(); // boundary values need not to be prolongated
		}
		else if(nodes[i].is_fine_direct())
		{
			// a non-interpolated fine node. calculate interpolation weights

			// calc min off-diag-entry, and sum of Neighbors
			double dmax = 0, connValue, maxConnValue = 0;
			double sumNeighbors =0, sumInterpolatory=0;

			double diag = amg_diag_value(A(i, i));

			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				if(conn.index() == i) continue; // skip diag
				connValue = amg_offdiag_value(conn.value());

				if(connValue > 0)
				{
					diag += connValue;
					continue;
				}

				sumNeighbors += connValue;

				if(dmax > connValue)
					dmax = connValue;
				if(nodes[conn.index()].is_coarse() && maxConnValue > connValue)
					maxConnValue = connValue;

			}

			double barrier;
			// todo: check if it is ok to do it THIS way:
			/*if(eps_truncation_of_interpolation > 0)  // [AMGKS99] 7.2.4 truncation of interpolation
				barrier = min(theta*dmax, eps_truncation_of_interpolation*maxConnValue);
			else*/
				barrier = theta*dmax;

			con.clear();
			// step 1: set w'_ij = a_ij/a_jj for suitable j
			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				if(conn.index() == i) continue; // skip diagonal
				if(!nodes[conn.index()].is_coarse()) continue;

				connValue = amg_offdiag_value(conn.value());
				if(connValue > barrier)
					continue;
				c.iIndex = newIndex[conn.index()];
				c.dValue = connValue;

				UG_ASSERT(c.iIndex >= 0, "not coarse?");

				con.push_back(c);
				sumInterpolatory += connValue;
			}

			if(con.size() > 0)
			{
				// step 2: calculate alpha_i
				double alpha = - (sumNeighbors / sumInterpolatory) / diag;
				// step 3: set w_ij = alpha * w'_ij = alpha * a_ii/a_jj.
				for(size_t j=0; j<con.size(); j++)
					con[j].dValue *= alpha;

				//UG_ASSERT(con.size() > 0, "0 connections in point i = " << i << " ?");
				// set w_ij in matrix row P[i] forall j.
				P.set_matrix_row(i, &con[0], con.size());
			}
			else
			{
				// no suitable interpolating nodes for node i,
				// so this node has to be treated by indirect interpolation
				nodes.set_unassigned_fine_indirect(i);
			}
		}
		else
		{
			UG_ASSERT(0, "?");
			//unassigned++;
			//UG_ASSERT(aggressiveCoarsening != 0, "no aggressive Coarsening but node " << i << " is fine and indirect??");
		}
	}

	if(nodes.get_unassigned())
		UG_DLOG(LIB_ALG_AMG, 1, "Pass 1: " << nodes.get_unassigned() << " left. ")
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateIndirectProlongation:
//-------------------------
/**
 * Assume Prolongation of all normal fine nodes is already computed, it calculates the Interpolation of
 * fineIndirect nodes with Matrix_type A and Coarse/FineIndirect markers in nodes[i].is_coarse/is_fine_indirect
 *
 * [AMGKS99] 7.2.2 Multi-pass-interpolation
 * Using A1 or A2 aggressive coarsening, it can be shown that we need at most 4 passes.
 * \note CreateIndirectProlongation is starts with pass 2, since pass 1 is CreateRugeStuebenProlongation.
 *
 * \param 	P					Matrix P: here goes the calculated prolongation
 * \param	A					Matrix A: matrix for which to calculate prolongation on next level
 * \param	nodes
 * \param 	posInConnections	array of size A.num_rows() for speedup of neighbor-neighbor-calculation inited with -1.
 * \param	theta
 */
template<typename Matrix_type>
void CreateIndirectProlongation(SparseMatrix<double> &P, const Matrix_type &A,
		AMGNodes &nodes, int *posInConnections, double theta)
{
	std::vector<SparseMatrix<double>::connection > con, con2;
	std::vector<int> nrOfPaths;
	con.reserve(255); con2.reserve(255); nrOfPaths.reserve(255);
	SparseMatrix<double>::connection c;
	//P.print();
	// INDIRECT INTERPOLATION

	size_t oldUnassigned = -1;
	int pass=2;
	while(nodes.get_unassigned())
	{
#ifdef AMG_PRINT_INDIRECT
		UG_DLOG(LIB_ALG_AMG, 1, std::endl);
#endif
		UG_DLOG(LIB_ALG_AMG, 1, "Pass " << pass << ": ");
		for(size_t i=0; i<A.num_rows() && nodes.get_unassigned() > 0; i++)
		{
			if(!nodes[i].is_unassigned_fine_indirect() || A.is_isolated(i))
				continue;

			double diag = amg_diag_value(A(i, i));
			// calculate min offdiag-entry
			double dmax = 0;

			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				if(conn.index() == i) continue; // skip diagonal
				double connValue = amg_offdiag_value(conn.value());
				if(connValue > 0)
				{
					diag += connValue;
					continue;
				}
				if(dmax > connValue)
					dmax = connValue;
			}

			con.clear();
			double sumInterpolatory=0, sumNeighbors=0;

			//cout << "indirect interpolating node " << i << endl;

			// look at neighbors of node i, try to interpolate indirectly through them
			for(typename Matrix_type::const_row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn)
			{
				size_t indexN = conn.index();
				if(indexN == i) continue; // skip diagonal

				double connValue = amg_offdiag_value(conn.value());
				sumNeighbors += connValue;
				if(connValue > theta * dmax)
					continue;

				// we dont want fine nodes which were indirectly interpolated in THIS pass
				// (no gauss-seidel-style)
				if(nodes[indexN].is_fine_indirect_level(pass))
					continue;
				// all interpolate neighbors are now from pass (pass-1)


				UG_ASSERT(!nodes[indexN].is_coarse(), "Node " << i << " indirect, but neighbor " <<  indexN << " coarse?");

				// now we look at the interpolation values of this neighbor (in P !!!)
				// and set w'(i,indexNN) += A(i,indexN) * P(indexN,indexNN)
				for(typename SparseMatrix<double>::row_iterator conn2 = P.begin_row(indexN); conn2 != P.end_row(indexN); ++conn2)
				{
					size_t indexNN = conn2.index();
					int pos = posInConnections[indexNN];

					if(pos == -1)
					{
						pos = posInConnections[indexNN] = con.size();
						c.iIndex = indexNN; assert(c.iIndex >= 0);

						AssignMult(c.dValue, connValue, conn2.value());
						con.push_back(c);
					}
					else
						AddMult(con[pos].dValue, connValue, conn2.value());
				}
			}

			for(size_t j=0; j<con.size(); j++)
			{
				sumInterpolatory += con[j].dValue;		// calc sumInterpolatory
				posInConnections[con[j].iIndex] = -1; 	// reset posInConnections
			}

			if(con.size() == 0)
				continue;

			nodes.set_fine_indirect_level(i, pass);

#ifdef AMG_PRINT_INDIRECT
			cout << i << " ";
#endif
			//cout << endl;

			UG_ASSERT(sumInterpolatory != 0.0, " numerical unstable?");

			double alpha = - (sumNeighbors / sumInterpolatory)/diag;

			for(size_t j=0; j<con.size(); j++)
			{
				// set w(i,j) = alpha * w'(i,j)
				//cout << con[j].dValue << " - N:" << sumNeighbors << " I: " << sumInterpolatory << " alpha: " << alpha << ". " << con[j].dValue*alpha << " : " << A.get_diag(i) << endl;
				con[j].dValue *= alpha;
			}

			// add connections
			P.set_matrix_row(i, &con[0], con.size());
		}


		if(nodes.get_unassigned() == oldUnassigned)
		{
			UG_LOG(std::endl << "unassigned nodes left: " << std::endl);
			for(size_t i=0; i<A.num_rows(); i++)
			{
				if(nodes[i].is_unassigned_fine_indirect())
				 UG_LOG(i << " ");
			}
			UG_LOG("\n");
		}
		UG_ASSERT(nodes.get_unassigned() != oldUnassigned, "Pass " << pass <<
				": Indirect Interpolation hangs at " << nodes.get_unassigned() << " unassigned nodes.");

#ifdef AMG_PRINT_INDIRECT
		UG_LOG("calculated, ");
#endif
		UG_LOG(nodes.get_unassigned() << " left. ");
		pass++;
		oldUnassigned = nodes.get_unassigned();
		break;
	}

	UG_ASSERT(nodes.get_unassigned() == 0, "number of unassigned nodes is still " << nodes.get_unassigned());

	P.defragment();
	//P.print();
}

} // namespace ug

#endif /* __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_RS_PROLONGATION_H__ */
