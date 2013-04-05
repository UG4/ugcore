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
		using base_type::m_spOperator;
		using base_type::damping;

	public:
	///	default constructor
		Jacobi() {this->set_damp(1.0);};

	///	constructor setting the damping parameter
		Jacobi(number damp) {this->set_damp(damp);};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<Jacobi<algebra_type> > newInst(new Jacobi<algebra_type>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(damping());
			return newInst;
		}

	///	Destructor
		virtual ~Jacobi()
		{};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "Jacobi";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{

			PROFILE_BEGIN_GROUP(Jacobi_preprocess, "algebra Jacobi");
#ifdef UG_PARALLEL
			matrix_type &mat = *pOp;
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
			diag.set_layouts(mat.layouts());

		// 	copy diagonal
			for(size_t i = 0; i < diag.size(); ++i){
				diag[i] = mat(i, i);
			}

		//	make diagonal consistent
			diag.set_storage_type(PST_ADDITIVE);
			diag.change_storage_type(PST_CONSISTENT);

		//	get damping in constant case to damp at once
			number damp = 1.0;
			if(damping()->constant_damping())
				damp = damping()->damping();

		// 	invert diagonal and multiply by damping
			for(size_t i = 0; i < diag.size(); ++i)
			{
				diag[i] *= 1./damp;
				GetInverse(m_diagInv[i], diag[i]);
			}
#endif
		//	done
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(Jacobi_step, "algebra Jacobi");
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
			matrix_type &mat = *pOp;
		//	get damping in constant case to damp at once
			number damp = 1.0;
			if(damping()->constant_damping())
				damp = damping()->damping();

		// 	apply iterator: c = B*d (damp is not used)
			for(size_t i = 0; i < c.size(); ++i)
			{
			// 	c[i] = damp * d[i] / mat(i,i)
				InverseMatMult(c[i], damp, mat(i,i), d[i]);
			}
#endif
		//	done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}


	//	overwrite function in order to specially treat constant damping
		virtual bool apply(vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(Jacobi_apply, "algebra Jacobi");
		//	Check that operator is initialized
			if(!this->m_bInit)
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Iterator not initialized.\n");
				return false;
			}

		//	Check parallel status
			#ifdef UG_PARALLEL
			if(!d.has_storage_type(PST_ADDITIVE))
				UG_THROW(name() << "::apply: Wrong parallel "
				               "storage format. Defect must be additive.");
			#endif

		//	Check sizes
			if(d.size() != m_spOperator->num_rows())
				UG_THROW("Vector [size= "<<d.size()<<"] and Row [size= "
				               <<this->m_spOperator->num_rows()<<"] sizes have to match!");
			if(c.size() != m_spOperator->num_cols())
				UG_THROW("Vector [size= "<<c.size()<<"] and Column [size= "
				               <<m_spOperator->num_cols()<<"] sizes have to match!");
			if(d.size() != c.size())
				UG_THROW("Vector [d size= "<<d.size()<<", c size = "
				               << c.size() << "] sizes have to match!");

		// 	apply iterator: c = B*d
			if(!step(m_spOperator, c, d))
			{
				UG_LOG("ERROR in '"<<name()<<"::apply': Step Routine failed.\n");
				return false;
			}

		//	apply scaling
			if(!damping()->constant_damping()){
				const number kappa = damping()->damping(c, d, m_spOperator);
				if(kappa != 1.0){
					c *= kappa;
				}
			}

		//	Correction is always consistent
			#ifdef 	UG_PARALLEL
			if(!c.change_storage_type(PST_CONSISTENT))
				UG_THROW(name() << "::apply': Cannot change "
						"parallel storage type of correction to consistent.");
			#endif

		//	we're done
			return true;
		}

	protected:
#ifdef UG_PARALLEL
	///	type of block-inverse
		typedef typename block_traits<typename matrix_type::value_type>::inverse_type inverse_type;

	///	storage of the inverse diagonal in parallel
		std::vector<inverse_type> m_diagInv;
#endif

};


} // end namespace ug

#endif
