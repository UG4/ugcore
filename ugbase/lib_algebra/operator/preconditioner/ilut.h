/*
 * \file ilut.h
 *
 * \date 04.07.2010
 * \author Martin Rupp
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILUT__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ILUT__

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/operator_interface.h"
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

	public:
	//	Constructor
		ILUTPreconditioner(double eps=1e-6) :
			m_eps(eps)
		{};

	// 	Clone
		ILinearIterator<vector_type,vector_type>* clone()
		{
			return new ILUTPreconditioner<algebra_type>(m_eps);
		}

	// 	Destructor
		virtual ~ILUTPreconditioner()
		{
			m_L.destroy();
			m_U.destroy();
		};

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUTPreconditioner";}

	//	Preprocess routine
		virtual bool preprocess()
		{
		//	Prepare Inverse Matrix
			matrix_type* A = this->m_pMatrix;
			m_L.destroy();	m_L.create(A->num_rows(), A->num_cols());
			m_U.destroy();	m_U.create(A->num_rows(), A->num_cols());

			// con is the current line of L/U
			std::vector<typename matrix_type::connection> con;
			con.reserve(300);
			con.resize(0);

			// init row 0 of U
			for(typename matrix_type::rowIterator i_it = A->beginRow(0); !i_it.isEnd(); ++i_it)
				con.push_back(*i_it);
			m_U.set_matrix_row(0, &con[0], con.size());

			size_t totalentries=0;
			size_t maxentries=0;
			for(size_t i=1; i<A->num_rows(); i++)
			{
				con.resize(0);
				size_t u_part=0;

				// get the row A(i, .) into con
				double dmax=0;
				for(typename matrix_type::rowIterator i_it = A->beginRow(i); !i_it.isEnd(); ++i_it)
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

					typename matrix_type::rowIterator k_it = m_U.beginRow(k); // upper row iterator
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

							typename matrix_type::connection c;
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

			return true;
		}

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
			// apply iterator: c = LU^{-1}*d (damp is not used)
			for(size_t i=0; i < m_L.num_rows(); i++)
				c[i] = d[i] - m_L[i]*c;
			for(size_t i=m_U.num_rows()-1; ; i--)
			{
				typename matrix_type::rowIterator it = m_U.beginRow(i);
				UG_ASSERT((*it).iIndex == i, "");
				double uii = (*it).dValue;
				double s = c[i];
				++it;
				for(; !it.isEnd(); ++it)
					s -= (*it).dValue * c[(*it).iIndex];
				c[i] = s/uii;
				if(i==0) break;
			}
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		matrix_type m_L;
		matrix_type m_U;
		double m_eps;
};


} // end namespace ug

#endif
