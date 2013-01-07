/*
 * \file ilut.h
 *
 * \date 04.07.2010
 * \author Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT__

#include "lib_algebra/operator/interface/operator.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TAlgebra>
class ILUTPreconditioner : public IPreconditioner<TAlgebra>
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
		ILUTPreconditioner(double eps=1e-6) :
			m_eps(eps), m_info(false)
		{};

	// 	Clone

		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ILUTPreconditioner<algebra_type> > newInst(new ILUTPreconditioner<algebra_type>(m_eps));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_info(m_info);
			return newInst;
		}

	// 	Destructor
		virtual ~ILUTPreconditioner()
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
		
	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUTPreconditioner";}

	//	Preprocess routine
		virtual bool preprocess(matrix_operator_type& mat)
		{
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

		//	Prepare Inverse Matrix
			matrix_type* A = &mat;
			typedef typename matrix_type::connection connection;
			m_L.resize(A->num_rows(), A->num_cols());
			m_U.resize(A->num_rows(), A->num_cols());

			// con is the current line of L/U
			std::vector<typename matrix_type::connection> con;
			con.reserve(300);
			con.resize(0);

			// init row 0 of U
			for(typename matrix_type::row_iterator i_it = A->begin_row(0); i_it != A->end_row(0); ++i_it)
				con.push_back(connection(i_it.index(), i_it.value()));
			m_U.set_matrix_row(0, &con[0], con.size());

			size_t totalentries=0;
			size_t maxentries=0;
			for(size_t i=1; i<A->num_rows(); i++)
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

				// eliminate all entries A(i, k) with k<i with rows U(k, .) and k<i
				for(size_t i_it = 0; i_it < con.size(); ++i_it)
				{
					size_t k = con[i_it].iIndex;
					if(k >= i)
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
								++j;
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

				totalentries+=con.size();
				if(maxentries < con.size()) maxentries = con.size();

				// safe L and U
				m_L.set_matrix_row(i, &con[0], u_part);
				m_U.set_matrix_row(i, &con[u_part], con.size()-u_part);

			}

			m_L.defragment();
			m_U.defragment();

			if (m_info==true){
				UG_LOG("\n	ILUT storage information:\n");
				UG_LOG("	A nr of connections: " << A->total_num_connections()  << "\n");
				UG_LOG("	L+U nr of connections: " << m_L.total_num_connections()+m_U.total_num_connections() << "\n");
				UG_LOG("	Increase factor: " << (float)(m_L.total_num_connections() + m_U.total_num_connections() )/A->total_num_connections() << "\n");
			}

			return true;
		}

	//	Stepping routine
		virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
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
				if (BlockNorm(uii)<m_small_lower){
					// set correction to zero
					c[i] = 0;
					if (BlockNorm(s)>m_small_upper){
						UG_LOG("Warning: zero entry in last row of U with corresponding non-zero rhs entry (" << BlockNorm(s) << ")\n");
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
		static const number m_small_lower;
		static const number m_small_upper;
};

template <typename TAlgebra>
const number ILUTPreconditioner<TAlgebra>::m_small_lower = 1e-9;

template <typename TAlgebra>
const number ILUTPreconditioner<TAlgebra>::m_small_upper = 1e-6;

} // end namespace ug

#endif
