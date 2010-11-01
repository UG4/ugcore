/*
 * jacobi.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__JACOBI__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__JACOBI__

#include "lib_algebra/operator/operator.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/* This Operator type behaves different on application. It not only computes v = L*u, but also changes u. */
/* It is used in iterative schemes. */
template <typename TAlgebra>
class JacobiPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
		JacobiPreconditioner() :
			m_damp(1.0), m_bOpChanged(true)
		{};

		JacobiPreconditioner(number damp) :
			m_damp(damp), m_bOpChanged(true)
		{};

		void set_damp(number damp) {m_damp = damp;}

	//	Clone
		virtual ILinearIterator<vector_type,vector_type>* clone()
		{
			JacobiPreconditioner<TAlgebra>* clone = new JacobiPreconditioner<TAlgebra>(m_damp);

			return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
		}

	// 	Destructor
		virtual ~JacobiPreconditioner()
		{};

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Jacobi";}

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat)
		{
		//	Currently remember that Operator has changed
			m_bOpChanged = true;

			return true;
		}

		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			// todo: 	this should be done in 'init'. It is currently placed here, since a
			//			ParallelMatrix is not yet implemented and we have no possibility to
			//			access the communicator needed below
			if(m_bOpChanged)
			{
				// create help vector to apply diagonal
				size_t size = mat.num_rows();
				if(size != mat.num_cols())
				{
					UG_LOG("Square Matrix needed for Jacobi Iteration.\n");
					return false;
				}

				m_diagInv.resize(size);
				m_diag.create(size);

				m_diag.set_slave_layout(c.get_slave_layout());
				m_diag.set_master_layout(c.get_master_layout());
				m_diag.set_communicator(c.get_communicator());

				// copy diagonal
				for(size_t i = 0; i < m_diag.size(); ++i){
					m_diag[i] = mat.get_diag(i);
				}

				//	make diagonal consistent
				m_diag.set_storage_type(PST_ADDITIVE);
				m_diag.change_storage_type(PST_CONSISTENT);

				// invert diagonal and multiply by damping
				for(size_t i = 0; i < m_diag.size(); ++i)
				{
					m_diag[i] *= 1/m_damp;
					GetInverse(m_diagInv[i], m_diag[i]);
				}
				m_diag.destroy();

				m_bOpChanged = false;
			}
#endif

#ifdef UG_PARALLEL
			// multiply defect with diagonal, c = damp * D^{-1} * d
			for(size_t i = 0; i < m_diagInv.size(); ++i)
			{
				// c[i] = m_diagInv[i] * d[i];
				MatMult(c[i], 1.0, m_diagInv[i], d[i]);
			}

		// 	the computed correction is additive
			c.set_storage_type(PST_ADDITIVE);

		//	we make it consistent
			if(!c.change_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in 'JacobiPreconditioner::apply': "
						"Cannot change parallel status of correction to consistent.\n");
				return false;
			}
#else
			// apply iterator: c = B*d (damp is not used)
			for(size_t i = 0; i < c.size(); ++i)
			{
				// c[i] = m_damp * d[i] / mat(i,i)
				InverseMatMult(c[i], m_damp, mat(i,i), d[i]);
			}
#endif
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		ParallelVector<Vector< typename matrix_type::value_type > > m_diag;
		std::vector< typename block_matrix_traits<typename matrix_type::value_type>::inverse_type > m_diagInv;
#endif
		number m_damp;

		bool m_bOpChanged;
};


} // end namespace ug

#endif
