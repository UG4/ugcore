/*
 * jacobi.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__JACOBI__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__JACOBI__

#include "lib_algebra/operator/interface/operator_iterator.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

///	Jacobi Preconditioner
template <typename TAlgebra>
class Jacobi : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		Jacobi() : m_damp(1.0) {};

	///	constructor setting the damping parameter
		Jacobi(number damp) :	m_damp(damp){};

	///	sets the damping parameter
		void set_damp(number damp) {m_damp = damp;}

	///	Clone
		virtual ILinearIterator<vector_type,vector_type>* clone()
		{
			Jacobi<algebra_type>* newInst = new Jacobi<algebra_type>(m_damp);
			newInst->set_debug(debug_writer());
			return newInst;
		}

	///	Destructor
		virtual ~Jacobi()
		{};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "Jacobi";}

	///	Preprocess routine
		virtual bool preprocess(matrix_operator_type& mat)
		{
#ifdef UG_PARALLEL
		// 	create help vector to apply diagonal
			size_t size = mat.num_rows();
			if(size != mat.num_cols())
			{
				UG_LOG("Square Matrix needed for Jacobi Iteration.\n");
				return false;
			}

		//	temporary vector for the diagonal
			ParallelVector<Vector< typename matrix_type::value_type > > diag;

		//	resize
			m_diagInv.resize(size);
			diag.resize(size);

		//	copy the layouts+communicator into the vector
			diag.set_layouts(mat.master_layout(), mat.slave_layout());
			diag.set_communicator(mat.communicator());

		// 	copy diagonal
			for(size_t i = 0; i < diag.size(); ++i){
				diag[i] = mat(i, i);
			}

		//	make diagonal consistent
			diag.set_storage_type(PST_ADDITIVE);
			diag.change_storage_type(PST_CONSISTENT);

		// 	invert diagonal and multiply by damping
			for(size_t i = 0; i < diag.size(); ++i)
			{
				diag[i] *= 1/m_damp;
				GetInverse(m_diagInv[i], diag[i]);
			}
#endif
		//	done
			return true;
		}

		virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
		// 	multiply defect with diagonal, c = damp * D^{-1} * d
		//	note, that the damping is already included in the inverse diagonal
			for(size_t i = 0; i < m_diagInv.size(); ++i)
			{
			// 	c[i] = m_diagInv[i] * d[i];
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
		// 	apply iterator: c = B*d (damp is not used)
			for(size_t i = 0; i < c.size(); ++i)
			{
			// 	c[i] = m_damp * d[i] / mat(i,i)
				InverseMatMult(c[i], m_damp, mat(i,i), d[i]);
			}
#endif
		//	done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
	///	damping parameter
		number m_damp;

#ifdef UG_PARALLEL
	///	type of block-inverse
		typedef typename block_traits<typename matrix_type::value_type>::inverse_type inverse_type;

	///	storage of the inverse diagonal in parallel
		std::vector<inverse_type> m_diagInv;
#endif

};


} // end namespace ug

#endif
