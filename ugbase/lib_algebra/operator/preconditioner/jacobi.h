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
		virtual const char* name() const {return "GSPreconditioner";}

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

				if(m_diagInv.size() == size) m_diagInv.set(0.0);
				else {
					m_diagInv.destroy();
					m_diagInv.create(size);
				}

				m_diagInv.set_slave_layout(c.get_slave_layout());
				m_diagInv.set_master_layout(c.get_master_layout());
				m_diagInv.set_communicator(c.get_communicator());

				// copy diagonal
				for(size_t i = 0; i < m_diagInv.size(); ++i){
					m_diagInv[i] = mat.get_diag(i);
				}

				//	make diagonal consistent
				m_diagInv.set_storage_type(PST_ADDITIVE);
				m_diagInv.change_storage_type(PST_CONSISTENT);

				// invert diagonal and multiply by damping
				for(size_t i = 0; i < m_diagInv.size(); ++i)
				{
					m_diagInv[i] *= 1/m_damp;
					Invert(m_diagInv[i]);
				}
				m_bOpChanged = false;
			}
#endif

#ifdef UG_PARALLEL
			// multiply defect with diagonal, c = damp * D^{-1} * d
			// TODO: We should handle this by a VecEntrywiseMultiply
			for(size_t i = 0; i < m_diagInv.size(); ++i)
			{
				c[i] = m_diagInv[i] * d[i];
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
				c[i] = d[i];
				c[i] /= mat.get_diag(i);
				// damp correction
				c[i] *= m_damp;
			}
#endif
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		ParallelVector<Vector<typename matrix_type::entry_type> > m_diagInv;
#endif
		number m_damp;

		bool m_bOpChanged;
};


} // end namespace ug

#endif
