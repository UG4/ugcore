/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PILUT__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PILUT__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "common/progress.h"
#include "common/util/ostream_util.h"
namespace ug{

template <typename TAlgebra>
class PILUTPreconditioner : public IPreconditioner<TAlgebra>
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

	private:
		typedef typename matrix_type::value_type block_type;
		using IPreconditioner<TAlgebra>::debug_writer;
		using IPreconditioner<TAlgebra>::set_debug;

	public:
	//	Constructor
		PILUTPreconditioner(double eps=0) :
			m_eps(eps), m_info(false), m_show_progress(true)
		{};

	// 	Clone

		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<PILUTPreconditioner<algebra_type> > newInst(new PILUTPreconditioner<algebra_type>(m_eps));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_info(m_info);
			return newInst;
		}

	// 	Destructor
		virtual ~PILUTPreconditioner()
		{
		};

	///	sets threshold for incomplete LU factorisation (added 01122010ih)
		void set_threshold(number thresh)
		{
			m_eps = thresh;
		}

	///	sets storage information output to true or false
		void set_info(bool info)
		{
			m_info = info;
		}
	
	///	set whether the progress meter should be shown
		void set_show_progress(bool s)
		{
			m_show_progress = s;
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "PILUTPreconditioner";}

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type mat;
			mat.set_as_copy_of(*pOp);
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

		//	Prepare Inverse Matrix
			matrix_type* A = &mat;
			typedef typename matrix_type::connection connection;
			m_L.resize_and_clear(A->num_rows(), A->num_cols());
			m_U.resize_and_clear(A->num_rows(), A->num_cols());

			// con is the current line of L/U
			// i also tried using std::list here or a special custom vector-based linked list
			// but vector is fastest, even with the insert operation.
			std::vector<typename matrix_type::connection> con;
			con.reserve(300);
			con.resize(0);

			// init row 0 of U
			for(typename matrix_type::row_iterator i_it = A->begin_row(0); i_it != A->end_row(0); ++i_it)
				con.push_back(connection(i_it.index(), i_it.value()));
			m_U.set_matrix_row(0, &con[0], con.size());

			size_t totalentries=0;
			size_t maxentries=0;
			Progress prog;
			if(m_show_progress)
				PROGRESS_START_WITH(prog, A->num_rows(),
					"Using ILUT(" << m_eps << ") on " << A->num_rows() << " x " << A->num_rows() << " matrix...");

			PROFILE_BEGIN(PILUT1);
			for(size_t i=1; i<A->num_rows()/2; i++)
			{
				if(m_show_progress) {PROGRESS_UPDATE(prog, i);}
				con.resize(0);
				size_t u_part=0;

				// get the row A(i, .) into con
				double dmax=0;
				for(typename matrix_type::row_iterator i_it = A->begin_row(i); i_it != A->end_row(i); ++i_it)
				{
					con.push_back(connection(i_it.index(), i_it.value()));
					if(dmax < BlockNorm(i_it.value()))
						dmax = BlockNorm(i_it.value());
				}


				u_part = eliminate_row(con, 0, i, dmax);
				totalentries+=con.size();
				if(maxentries < con.size()) maxentries = con.size();
				// safe L and U
				for(size_t i_it=0; i_it<con.size(); i_it++)
				{
					if(con[i_it].iIndex < i) m_L(i, con[i_it].iIndex) = con[i_it].dValue;
					else m_U(i, con[i_it].iIndex) = con[i_it].dValue;
				}
			}
			PROFILE_END();
			PROFILE_BEGIN(PILUT2);
			for(size_t i=A->num_rows()/2; i<A->num_rows(); i++)
			{
				con.resize(0);
				size_t u_part=0;

				// get the row A(i, .) into con
				double dmax=0;
				for(typename matrix_type::row_iterator i_it = A->begin_row(i); i_it != A->end_row(i); ++i_it)
				{
					con.push_back(connection(i_it.index(), i_it.value()));
					if(dmax < BlockNorm(i_it.value()))
						dmax = BlockNorm(i_it.value());
				}


				u_part = eliminate_row(con, 0, A->num_rows()/2, dmax);
				totalentries+=con.size();
				if(maxentries < con.size()) maxentries = con.size();
				// safe L and U
				A->set_matrix_row(i, &con[0], con.size());
			}
			con.clear();
			size_t u_part;
			for(typename matrix_type::row_iterator i_it = A->begin_row(A->num_rows()/2); i_it != A->end_row(A->num_rows()/2); ++i_it)
			{
				if(i_it.index() <= A->num_rows()/2)
					u_part=con.size();
				con.push_back(connection(i_it.index(), i_it.value()));
			}
			m_L.set_matrix_row(A->num_rows()/2, &con[0], u_part);
			m_U.set_matrix_row(A->num_rows()/2, &con[u_part], con.size()-u_part);
			PROFILE_END();
			PROFILE_BEGIN(PILUT3);
			for(size_t i=A->num_rows()/2+1; i<A->num_rows(); i++)
			{
				if(m_show_progress) {PROGRESS_UPDATE(prog, i);}
				con.resize(0);
				size_t u_part=0;

				// get the row A(i, .) into con
				double dmax=0;
				for(typename matrix_type::row_iterator i_it = A->begin_row(i); i_it != A->end_row(i); ++i_it)
				{
					con.push_back(connection(i_it.index(), i_it.value()));
					if(dmax < BlockNorm(i_it.value()))
						dmax = BlockNorm(i_it.value());
				}


				u_part = eliminate_row(con, A->num_rows()/2, i, dmax);
				totalentries+=con.size();
				if(maxentries < con.size()) maxentries = con.size();
				// safe L and U
				for(size_t i_it=0; i_it<con.size(); i_it++)
				{
					if(con[i_it].iIndex < i) m_L(i, con[i_it].iIndex) = con[i_it].dValue;
					else m_U(i, con[i_it].iIndex) = con[i_it].dValue;
				}
			}
			PROFILE_END();
			//m_L.print();
			//m_U.print();

			if(m_show_progress) {PROGRESS_FINISH(prog);}

			m_L.defragment();
			m_U.defragment();

			if (m_info==true){
				UG_LOG("\n	ILUT storage information:\n");
				UG_LOG("	A nr of connections: " << A->total_num_connections()  << "\n");
				UG_LOG("	L+U nr of connections: " << m_L.total_num_connections()+m_U.total_num_connections() << "\n");
				UG_LOG("	Increase factor: " << (float)(m_L.total_num_connections() + m_U.total_num_connections() )/A->total_num_connections() << "\n");
				UG_LOG(reset_floats << "Total entries: " << totalentries << " (" << ((double)totalentries) / (A->num_rows()*A->num_rows()) << "% of dense)");
			}

			return true;
		}

		size_t eliminate_row(std::vector<typename matrix_type::connection> &con, size_t start, size_t stop, double dmax)
		{
			size_t u_part;
			// eliminate all entries A(i, k) with k<i with rows U(k, .) and k<i
			for(size_t i_it = 0; i_it < con.size(); ++i_it)
			{
				size_t k = con[i_it].iIndex;
				if(k < start) continue;
				if(k >= stop)
				{
					// safe where U begins / L ends in con
					u_part = i_it;
					break;
				}
				if(con[i_it].dValue == 0.0) continue;
				UG_ASSERT(m_U.num_connections(k) != 0 && m_U.begin_row(k).index() == k, "");
				block_type &ukk = m_U.begin_row(k).value();

				// add row k to row i by A(i, .) -= U(k,.)  A(i,k) / U(k,k)
				// so that A(i,k) is zero.
				// safe A(i,k)/U(k,k) in con, (later L(i,k) )
				block_type &d = con[i_it].dValue = con[i_it].dValue / ukk;

				typename matrix_type::row_iterator k_it = m_U.begin_row(k); // upper row iterator
				++k_it; // skip diag
				size_t j = i_it+1;
				while(k_it != m_U.end_row(k) && j < con.size())
				{
					// (since con and U[k] is sorted, we can do sth like a merge on the two lists)
					if(k_it.index() == con[j].iIndex)
					{
						// match
						con[j].dValue -= k_it.value() * d;
						++k_it;	++j;
					}
					else if(k_it.index() < con[j].iIndex)
					{
						// we have a value in U(k, (*k_it).iIndex), but not in A.
						// check tolerance criteria

						typename matrix_type::connection
							c(k_it.index(), k_it.value() * d * -1.0);

						if(BlockNorm(c.dValue) > dmax * m_eps)
						{
							// insert sorted
							con.insert(con.begin()+j, c);
							++j; // don't do this when using iterators
						}
						// else do some lumping
						++k_it;
					}
					else
					{
						// we have a value in A(k, con[j].iIndex), but not in U.
						++j;
					}
				}
				// insert new connections after last connection of row i
				if (k_it!=m_U.end_row(k)){
					for (;k_it!=m_U.end_row(k);++k_it){
						typename matrix_type::connection c(k_it.index(),-k_it.value()*d);
						if(BlockNorm(c.dValue) > dmax * m_eps)
						{
							con.push_back(c);
						}
					}
				};
			}
			return u_part;

		}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			// apply iterator: c = LU^{-1}*d (damp is not used)
			// L
			for(size_t i=0; i < m_L.num_rows(); i++)
			{
				// c[i] = d[i] - m_L[i]*c;
				c[i] = d[i];
				for(typename matrix_type::row_iterator it = m_L.begin_row(i); it != m_L.end_row(i); ++it)
					MatMultAdd(c[i], 1.0, c[i], -1.0, it.value(), c[it.index()] );
				// lii = 1.0.
			}
			// U
			//
			// last row diagonal U entry might be close to zero with corresponding zero rhs
			// when solving Navier Stokes system, therefore handle separately
			{
				size_t i=m_U.num_rows()-1;
				typename matrix_type::row_iterator it = m_U.begin_row(i);
				UG_ASSERT(it.index() == i, "");
				block_type &uii = it.value();
				typename vector_type::value_type s = c[i];
				// check if close to zero
				if (BlockNorm(uii)<m_small){
					// set correction to zero
					c[i] = 0;
					if (BlockNorm(s)>m_small){
						UG_LOG("Warning: near-zero diagonal entry in last row of U with corresponding non-near-zero rhs entry (" << BlockNorm(s) << ")\n");
					}
				} else {
					// c[i] = s/uii;
					InverseMatMult(c[i], 1.0, uii, s);
				}
			}
			// handle all other rows
			for(size_t i=m_U.num_rows()-2; ; i--)
			{
				typename matrix_type::row_iterator it = m_U.begin_row(i);
				UG_ASSERT(it.index() == i, "");
				block_type &uii = it.value();

				typename vector_type::value_type s = c[i];
				++it; // skip diag
				for(; it != m_U.end_row(i); ++it)
					// s -= it.value() * c[it.index()];
					MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()] );

					// c[i] = s/uii;
					InverseMatMult(c[i], 1.0, uii, s);

					if(i==0) break;
			}

#ifdef 	UG_PARALLEL
		//	Correction is always consistent
		//	todo: We set here correction to consistent, but it is not. Think about how to use ilu in parallel.
			c.set_storage_type(PST_CONSISTENT);
#endif
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		matrix_type m_L;
		matrix_type m_U;
		double m_eps;
		bool m_info;
		bool m_show_progress;
		static const number m_small;
};

// define constant
template <typename TAlgebra>
const number PILUTPreconditioner<TAlgebra>::m_small = 1e-8;

} // end namespace ug

#endif
