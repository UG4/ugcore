/*
 * ilu.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILUT__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILUT__

#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/operator/operator_interface.h"

namespace ug{

template <typename TDiscreteFunction>
class AssembledILUTOperator : public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// domain function type
		typedef TDiscreteFunction domain_function_type;

		// codomain function type
		typedef TDiscreteFunction codomain_function_type;

	private:
		typedef typename TDiscreteFunction::algebra_type algebra_type;

	public:
		AssembledILUTOperator(double eps=1e-6) : m_bOpChanged(true), m_eps(eps)
		{
			//UG_LOG("INIT");
		};

		bool init(ILinearizedOperator<domain_function_type,codomain_function_type>& Op)
		{
			AssembledLinearizedOperator<TDiscreteFunction>* A
				= dynamic_cast<AssembledLinearizedOperator<TDiscreteFunction>*>(&Op);
			UG_ASSERT(A != NULL, 	"Operator used does not use a matrix."
									" Currently only matrix based operators can be inverted by this Jacobi.\n");

			// remember operator
			m_pOperator = A;
			m_pMatrix = &m_pOperator->get_matrix();
			m_bOpChanged = true;
			return true;
		}

		// prepare Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c)
		{
			if(m_bOpChanged)
			{
				m_L.destroy();	m_L.create(m_pMatrix->num_rows(), m_pMatrix->num_cols());
				m_U.destroy();	m_U.create(m_pMatrix->num_rows(), m_pMatrix->num_cols());

				typedef typename algebra_type::matrix_type Matrix_type;

				// con is the current line of L/U
				std::vector<typename Matrix_type::connection> con;
				con.reserve(300);	
				con.resize(0);
				
				// init row 0 of U
				for(typename Matrix_type::rowIterator i_it = m_pMatrix->beginRow(0); !i_it.isEnd(); ++i_it)
					con.push_back(*i_it);
				m_U.set_matrix_row(0, &con[0], con.size());
				
				size_t totalentries=0;
				size_t maxentries=0;
				for(size_t i=1; i<m_pMatrix->num_rows(); i++)
				{
					con.resize(0);
					size_t u_part=0;

					// get the row A(i, .) into con
					double dmax=0;
					for(typename Matrix_type::rowIterator i_it = m_pMatrix->beginRow(i); !i_it.isEnd(); ++i_it)
					{
						con.push_back(*i_it);
						if(dmax < dabs((*i_it).dValue))
							dmax = dabs((*i_it).dValue);
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
						UG_ASSERT(!m_U.beginRow(k).isEnd() && (*m_U.beginRow(k)).iIndex == k, "");
						double ukk = (*m_U.beginRow(k)).dValue;
						double d = con[i_it].dValue / ukk;
						
						// add row k to row i by A(i, .) -= U(k,.)  A(i,k) / U(k,k)
						// so that A(i,k) is zero.
						// safe A(i,k)/U(k,k) in con, (later L(i,k) )
						con[i_it].dValue = d;
						
						typename Matrix_type::rowIterator k_it = m_U.beginRow(k); // upper row iterator
						++k_it; // skip diag
						size_t j = i_it+1;
						while(!k_it.isEnd() && j < con.size())
						{
							// (since con and U[k] is sorted, we can do sth like a merge on the two lists)
							if((*k_it).iIndex == con[j].iIndex)
							{
								// match
								con[j].dValue -= (*k_it).dValue * d;
								++k_it;	++j;
							}
							else if((*k_it).iIndex < con[j].iIndex)
							{
								// we have a value in U(k, (*k_it).iIndex), but not in A.
								// check tolerance criteria

								typename Matrix_type::connection c;
								c.iIndex = (*k_it).iIndex;
								c.dValue = - (*k_it).dValue * d;
								if(dabs(c.dValue) > dmax * m_eps)
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
					}
					
					totalentries+=con.size();
					if(maxentries < con.size()) maxentries = con.size();
					
					// safe L and U
					m_L.set_matrix_row(i, &con[0], u_part);
					m_U.set_matrix_row(i, &con[u_part], con.size()-u_part);
				}	
				
				m_L.finalize();
				m_U.finalize();

				m_bOpChanged = false;
			}

			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c)
		{
			if(!d.has_storage_type(PST_ADDITIVE)) return false;

			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();


			// Apply ILU on d:   c := ILU^{-1}*d

			// make correction consistent
			c.set_storage_type(PST_ADDITIVE);
			if(c.change_storage_type(PST_CONSISTENT) != true)
				return false;

			// apply iterator: c = LU^{-1}*d (damp is not used)
			for(size_t i=0; i < m_L.num_rows(); i++)
				c_vec[i] = d_vec[i] - m_L[i]*c_vec;
			for(size_t i=m_U.num_rows()-1; ; i--)
			{
				typename algebra_type::matrix_type::rowIterator it = m_U.beginRow(i);
				UG_ASSERT((*it).iIndex == i, "");
				double uii = (*it).dValue;
				double s = c_vec[i];
				++it;
				for(; !it.isEnd(); ++it)
					s -= (*it).dValue * c_vec[(*it).iIndex];
				c_vec[i] = s/uii;
				if(i==0) break;
			}		

			// update defect
			// TODO: Check that matrix has correct type (additive)
			m_pMatrix->matmul_minus(d_vec, c_vec);

			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
			d.set_storage_type(PST_ADDITIVE);
			return true;
		}

		// clone
		ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
		{
			return new AssembledILUTOperator<TDiscreteFunction>(m_eps);
		}


		// destructor
		virtual ~AssembledILUTOperator() {};

	protected:
		AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		typename algebra_type::matrix_type m_L;
		typename algebra_type::matrix_type m_U;
	
		bool m_bOpChanged;
		double m_eps;
};


} // end namespace ug

#endif
