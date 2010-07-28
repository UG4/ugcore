/*
 * ilu.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILU__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILU__

#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/operator/operator_interface.h"

namespace ug{

template<typename Matrix_type>
bool FactorizeILU(Matrix_type &A)
{


	for(size_t i=1; i < A.num_rows(); i++)
	{
		for(typename Matrix_type::rowIterator it_k = A.beginRow(i); !it_k.isEnd() && ((*it_k).iIndex < i); ++it_k)
		{
			const size_t k = (*it_k).iIndex;
			typename Matrix_type::entry_type a_ik = (*it_k).dValue;
			if(BlockNorm(a_ik) < 1e-7)	continue;
			typename Matrix_type::entry_type a_kk = A(k,k);

			a_ik /= a_kk;

			 (*it_k).dValue = a_ik;

			typename Matrix_type::rowIterator it_j = it_k;
			for(++it_j; !it_j.isEnd(); ++it_j)
			{
				const size_t j = (*it_j).iIndex;
				typename Matrix_type::entry_type& a_ij = (*it_j).dValue;
				bool bFound;
				typename Matrix_type::rowIterator p = A.get_connection(k,j, bFound);
				if(bFound)
				{
					const typename Matrix_type::entry_type a_kj = (*p).dValue;
					a_ij -= a_ik*a_kj;
				}
			}
		}
	}

	return true;
}

// solve x = L^-1 b
template<typename Matrix_type, typename Vector_type>
bool invert_L(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	typename Vector_type::entry_type s;
	for(size_t i=0; i < x.size(); i++)
	{
		s = b[i];
		for(typename Matrix_type::cRowIterator it = A.beginRow(i); !it.isEnd(); ++it)
		{
			if((*it).iIndex >= i) continue;
			s -= (*it).dValue * x[(*it).iIndex];
		}
		x[i] = s;
	}

	return true;
}

// solve x = U^-1 * b
template<typename Matrix_type, typename Vector_type>
bool invert_U(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	typename Vector_type::entry_type s;
	for(size_t i = x.size()-1; ; --i)
	{
		s = b[i];
		for(typename Matrix_type::cRowIterator it = A.beginRow(i); !it.isEnd(); ++it)
		{
			if((*it).iIndex <= i) continue;
			s -= (*it).dValue * x[(*it).iIndex];
		}
		x[i] = s/A(i,i);
		if(i == 0) break;
	}

	return true;
}


template <typename TDiscreteFunction>
class AssembledILUOperator : public ILinearizedIteratorOperator<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// domain function type
		typedef TDiscreteFunction domain_function_type;

		// codomain function type
		typedef TDiscreteFunction codomain_function_type;

	private:
		typedef typename TDiscreteFunction::algebra_type algebra_type;

	public:
		AssembledILUOperator() : m_bOpChanged(true)
		{};

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
				m_ILU.destroy();
				m_h.destroy();

				// copy matrix
				typename algebra_type::matrix_type& A = *m_pMatrix;
				m_ILU.create(A.num_rows(), A.num_cols());
				for(size_t i=0; i < A.num_rows(); i++)
				{
					for(typename algebra_type::matrix_type::rowIterator it_k = A.beginRow(i); !it_k.isEnd(); ++it_k)
					{
						const size_t k = (*it_k).iIndex;
						m_ILU(i,k) = (*it_k).dValue;
					}
				}

				// Compute ILU Factorization
				FactorizeILU(m_ILU);

				m_h.create(A.num_cols());
				m_bOpChanged = false;
			}

			return true;
		}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c, bool updateDefect)
		{
#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE)) return false;
#endif
			typename domain_function_type::vector_type& d_vec = d.get_vector();
			typename codomain_function_type::vector_type& c_vec = c.get_vector();

			// Apply ILU on d:   c := ILU^{-1}*d

			// make correction consistent
#ifdef UG_PARALLEL
			c.set_storage_type(PST_ADDITIVE);
			if(c.change_storage_type(PST_CONSISTENT) != true)
				return false;
#endif
			// apply iterator: c = LU^{-1}*d (damp is not used)
			invert_L(m_ILU, m_h, d_vec); // h := L^-1 d
			invert_U(m_ILU, c_vec, m_h); // c := U^-1 h = (LU)^-1 d

			// update defect
			// TODO: Check that matrix has correct type (additive)
			if(updateDefect) m_pMatrix->matmul_minus(d_vec, c_vec);

			// defect is now no longer unique (maybe we should handle this in matrix multiplication)
#ifdef UG_PARALLEL
			d.set_storage_type(PST_ADDITIVE);
#endif
			return true;
		}

		// clone
		ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* clone()
		{
			return new AssembledILUOperator<TDiscreteFunction>();
		}


		// destructor
		virtual ~AssembledILUOperator() {};

	protected:
		AssembledLinearizedOperator<TDiscreteFunction>* m_pOperator;

		typename algebra_type::matrix_type* m_pMatrix;

		typename algebra_type::matrix_type m_ILU;
		typename algebra_type::vector_type m_h;

		bool m_bOpChanged;
};


} // end namespace ug

#endif
