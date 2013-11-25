/*
 * jacobi.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__JACOBI__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__JACOBI__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/small_algebra/additional_math.h"
#include "lib_algebra/cpu_algebra/vector.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/////////////////////////////////////////////////////////////////////////////////////////////
///		Jacobi-Iteration
/**
 * Here, the Jacobi-iteration is described for solving the linear equation
 *
 * 		\f$ A * x = b.			A \in R^{nxn}, x \in R^n, b \in R^n \f$.
 *
 * 	Most of the common linear iteration-methods base on the decomposition of A into
 * 	its diagonal (D) and strict-upper(-U) and strict-lower part (-L),
 *
 * 		\f$ A = D - L - U \f$.
 *
 * 	Among others, W. Hackbusch ('Iterative Loesung grosser Gleichungssysteme'),
 * 	distinguishes three different forms for describing a linear iteration scheme.
 * 	The general 'first normal-form' of a linear iteration scheme takes the form
 *
 * 		\f$ x^{m+1} = M * x^m + N * b \f$,
 *
 * 	with some Matrices \f$ M \f$ and \f$ N \in R^{nxn} \f$. m denotes the iteration index.
 * 	The general 'second normal-form' of a linear iteration scheme takes the form
 *
 * 		\f$ x^{m+1} = x^m - N * (A * x^m - b) \f$.
 *
 * 	Those linear iteration schemes, which can be represented by the second normal-form
 * 	are the linear, consistent iteration schemes.
 *
 * 	Introducing the correction \f$ c{m+1} := x^{m+1} - x^m \f$ and the defect
 * 	\f$ d^m := b - A * x^m \f$ the second normal-form can be rewritten as
 *
 * 		\f$ c = N * d \f$.
 *
 *	The matrix of the second normal-form for the Jacobi-method takes the simple form
 *
 *		\f$ N = D^{-1} \f$. 			.
 *
 *	References:
 * <ul>
 * <li> W. Hackbusch. Iterative Loesung grosser Gleichungssysteme
 * </ul>
 */

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
			newInst->set_block(m_bBlock);
			return newInst;
		}

	///	Destructor
		virtual ~Jacobi()
		{};

		void set_block(bool b)
		{
			m_bBlock = b;
		}

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "Jacobi";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{

			PROFILE_BEGIN_GROUP(Jacobi_preprocess, "algebra Jacobi");

			matrix_type &mat = *pOp;
		// 	create help vector to apply diagonal
			size_t size = mat.num_rows();
			if(size != mat.num_cols())
			{
				UG_LOG("Square Matrix needed for Jacobi Iteration.\n");
				return false;
			}

			//	resize
			m_diagInv.resize(size);
#ifdef UG_PARALLEL
					//	temporary vector for the diagonal
			ParallelVector<Vector< typename matrix_type::value_type > > diag;
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


			if(diag.size() > 0)
				if(CheckVectorInvertible(diag) == false)
					return false;
//			UG_ASSERT(CheckVectorInvertible(diag), "Jacobi: A has noninvertible diagonal");

#endif

//	get damping in constant case to damp at once
			number damp = 1.0;
			if(damping()->constant_damping())
				damp = damping()->damping();

			typename matrix_type::value_type m;
		// 	invert diagonal and multiply by damping
			for(size_t i = 0; i < mat.num_rows(); ++i)
			{
#ifdef UG_PARALLEL
				typename matrix_type::value_type &d = diag[i];
#else
				typename matrix_type::value_type &d = mat(i,i);
#endif
				if(m_bBlock)
					GetDiag(m, d);
				else
					m = d;
				m *= 1./damp;
				GetInverse(m_diagInv[i], m);
			}

		//	done
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(Jacobi_step, "algebra Jacobi");

		// 	multiply defect with diagonal, c = damp * D^{-1} * d
		//	note, that the damping is already included in the inverse diagonal
			for(size_t i = 0; i < m_diagInv.size(); ++i)
			{
			// 	c[i] = m_diagInv[i] * d[i];
				MatMult(c[i], 1.0, m_diagInv[i], d[i]);
			}

#ifdef UG_PARALLEL

		// 	the computed correction is additive
			c.set_storage_type(PST_ADDITIVE);

		//	we make it consistent
			if(!c.change_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in 'JacobiPreconditioner::apply': "
						"Cannot change parallel status of correction to consistent.\n");
				return false;
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
			THROW_IF_NOT_EQUAL(d.size(), m_spOperator->num_rows());
			THROW_IF_NOT_EQUAL(c.size(), m_spOperator->num_cols());
			THROW_IF_NOT_EQUAL(c.size(), m_spOperator->num_cols());
			THROW_IF_NOT_EQUAL(c.size(), d.size());

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
	///	type of block-inverse
		typedef typename block_traits<typename matrix_type::value_type>::inverse_type inverse_type;

	///	storage of the inverse diagonal in parallel
		std::vector<inverse_type> m_diagInv;
		bool m_bBlock;


};

} // end namespace ug

#include "gpujacobi.h"

#endif
