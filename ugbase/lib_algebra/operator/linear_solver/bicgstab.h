/*
 * bicgstab.h
 *
 *  Created on: 05.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__

#include <iostream>
#include <sstream>

#include "lib_algebra/operator/operator_interface.h"
#include "common/profiler/profiler.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

///	the BiCGStab method as a solver for linear operators
/**
 * This class implements the BiCGStab - method for the solution of linear
 * operator problems like A*x = b, where the solution x = A^{-1} b is computed.
 *
 * For detailed description of the algorithm, please refer to:
 *
 * - Barrett, Berry, Chan, Demmel, Donatom Dongarra, Eijkhout, Pozo, Romine,
 * 	 Van der Vorst, "Templates for the Solution of Linear Systems: Building
 * 	 Blocks for Iterative Methods", p.24, Fig, 2.10
 *
 * - Saad, "Iterative Methods For Sparse Linear Systems", p246, Alg. 7.7
 *
 * \tparam 	TAlgebra		algebra type
 */
template <typename TAlgebra>
class BiCGStab :
	public ILinearOperatorInverse< 	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	public:
	///	default constructor
		BiCGStab() : m_pPrecond(NULL), m_pConvCheck(NULL){};

	///	constructor setting the preconditioner and the convergence check
		BiCGStab	( 	ILinearIterator<vector_type,vector_type>* Precond,
							IConvergenceCheck& ConvCheck) :
							m_pPrecond(Precond), m_pConvCheck(&ConvCheck)
		{};

	///	name of solver
		virtual const char* name() const {return "BiCGStab";}

	///	set the convergence check
		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
		}

	///	returns the convergence check
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}

	///	sets the preconditioner
		void set_preconditioner(ILinearIterator<vector_type, vector_type>& precond)
		{
			m_pPrecond = &precond;
		}

	///	initializes the solver
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
		//	remember operator
			m_A = &J;

		// 	init Preconditioner for operator J
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(J, u))
				{
					UG_LOG("ERROR in 'BiCGStabSolver::prepare': Cannot init "
							"Iterator Operator for Operator J.\n");
					return false;
				}

		//	done
			return true;
		}

	///	initializes the solver
		virtual bool init(ILinearOperator<vector_type, vector_type>& L)
		{
		//	remember operator
			m_A = &L;

		// 	init Preconditioner for operator L
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(L))
				{
					UG_LOG("ERROR in 'BiCGStabSolver::prepare': "
							"Cannot init Iterator Operator for Operator L.\n");
					return false;
				}

		//	done
			return true;
		}

	// 	Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
		//	convergence check is required
			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'BiCGStabSolver::apply_return_defect': "
						"Convergence check not set.\n");
				return false;
			}

		//	check correct storage type in parallel
			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("ERROR: In 'BiCGStabSolver::apply_return_defect':"
							"Inadequate storage format of Vectors.\n");
					return false;
				}
			#endif

		// 	build defect:  b := b - A*x
			m_A->apply_sub(b, x);

			// create start r_0^* vector
		//	todo: 	It would be sufficient to copy only the pattern and
		//			without initializing
			vector_type r; r.create(b.size()); r = b;
			vector_type p; p.create(b.size()); p = b;
			vector_type v; v.create(b.size()); v = b;
			vector_type t; t.create(b.size()); t = b;
			vector_type s; s.create(b.size()); s = b;
			vector_type q; q.create(x.size()); q = x;

		//	prepare convergence check
			prepare_conv_check();

		//	compute start defect norm
			m_pConvCheck->start(b);

		//	convert b to unique (should already be unique due to norm calculation)
			#ifdef UG_PARALLEL
			if(!b.change_storage_type(PST_UNIQUE))
			{
				UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
						"Cannot convert b to unique vector.\n");
				return false;
			}
			#endif

		//	needed variables
			number rhoOld=1,rho=1,alpha=1, omega=1;

		// 	Iteration loop
			while(!m_pConvCheck->iteration_ended())
			{
			//	check if start values have to be set
				if(m_pConvCheck->step() == 0 /*or restart*/)
				{
				//	if at restart recompute start defect
					if(m_pConvCheck->step() != 0)
					{
						m_pConvCheck->update(b);
						if(m_pConvCheck->iteration_ended()) break;
					}

				// 	reset vectors
					r = b;

				// 	make r additive unique
					#ifdef UG_PARALLEL
					if(!r.change_storage_type(PST_UNIQUE))
					{
						UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
								"Cannot convert r to unique vector.\n");
						return false;
					}
					#endif

				//	set start vectors
					p.set(0.0);
					v.set(0.0);

				//	set start values
					rhoOld = alpha = omega = 1.0;
				}

			// 	Compute rho new
				rho = VecProd(b, r);

			//	check values
				if(rhoOld == 0.0 || omega == 0.0)
				{
					UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
							"rho= "<<rhoOld<<" and omega= "<<omega<<" are invalid "
							"values. Aborting iteration.\n");
					return false;
				}

			// 	Compute new beta
				const number beta = (rho/rhoOld) * (alpha/omega);

			// 	scale p := beta * p
				p *= beta;

			// 	add b to p (p:= p + b)
				p += b;

			// 	subtract: p := p - beta * omega * v
				VecScaleAppend(p, v, (-1)*beta*omega);

			// 	Precondition
				if(m_pPrecond)
				{
				// 	apply q = M^-1 * p
					if(!m_pPrecond->apply(q, p))
					{
						UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
								"Cannot apply preconditioner. Aborting.\n");
						return false;
					}
				}
				else
				{
				// 	copy q = p
					q = p;

				// 	make q consistent
					#ifdef UG_PARALLEL
					if(!q.change_storage_type(PST_CONSISTENT))
					{
						UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
								"Cannot convert q to consistent vector.\n");
						return false;
					}
					#endif
				}

			// 	compute v := A*q
				m_A->apply(v, q);

			// 	make v unique
				#ifdef UG_PARALLEL
				if(!v.change_storage_type(PST_UNIQUE))
				{
					UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
							"Cannot convert v to unique vector.\n");
					return false;
				}
				#endif

			//	alpha = (v,r)
				alpha = VecProd(v, r);

			//	check validity of alpha
				if(alpha == 0.0)
				{
					UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
							"alpha= "<<alpha<<" is an invalid value."
							" Aborting iteration.\n");
					return false;
				}

			//	alpha = rho/(v,r)
				alpha = rho/alpha;

			// 	add: x := x + alpha * q
				VecScaleAppend(x, q, alpha);

			// 	set s := b
				s = b;

			// 	update s := s - alpha*v
				VecScaleAppend(s, v, (-1)*alpha);

			// 	check convergence
				m_pConvCheck->update(s);
				if(m_pConvCheck->iteration_ended())
				{
				//	set output to last defect
					b = s; break;
				}

			//	apply preconditioner
				if(m_pPrecond)
				{
				// 	apply q = M^-1 * t
					if(!m_pPrecond->apply(q, s))
					{
						UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
								"Cannot apply preconditioner. Aborting.\n");
						return false;
					}
				}
				else
				{
				// 	set q:=s
					q = s;

				// 	make q consistent
					#ifdef UG_PARALLEL
					if(!q.change_storage_type(PST_CONSISTENT))
					{
						UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
								"Cannot convert q to consistent vector.\n");
						return false;
					}
					#endif
				}

			// 	compute t := A*q
				m_A->apply(t, q);

			// 	make t unique
				#ifdef UG_PARALLEL
				if(!t.change_storage_type(PST_UNIQUE))
				{
					UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
							"Cannot convert t to unique vector.\n");
					return false;
				}
				#endif

			// 	tt = (t,t)
				const number tt = VecProd(t, t);

			// 	omega = (s,t)
				omega = VecProd(s, t);

			//	check tt
				if(tt == 0.0)
				{
					UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
							"tt= "<<tt<<" is an invalid value. "
							"Aborting iteration.\n");
					return false;
				}

			// 	omega = (s,t)/(t,t)
				omega = omega/tt;

			// 	add: x := x + omega * q
				VecScaleAppend(x, q, omega);

			// 	set b := s
				b = s;

			// 	update of b:  b:= b - omega*t
				VecScaleAppend(b, t, (-1)*omega);

			// 	check convergence
				m_pConvCheck->update(b);

			// 	remember current rho
				rhoOld = rho;
			}

		//	print ending output
			return m_pConvCheck->post();
		}

	///	apply the solver
		virtual bool apply(vector_type& x, const vector_type& b)
		{
		//	copy defect
			vector_type bTmp; bTmp.resize(b.size()); bTmp = b;

		//	solve on copy of defect
			return apply_return_defect(x, bTmp);
		}


	/// destructor
		virtual ~BiCGStab() {};

	protected:
	///	prepares the output of the convergence check
		void prepare_conv_check()
		{
		//	set iteration symbol an name
			m_pConvCheck->set_name(name());
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());

		//	set preconditioner string
			std::stringstream ss;
			if(m_pPrecond) ss<<" (Precond: "<<m_pPrecond->name()<<")";
			else ss << " (No Preconditioner) ";
			m_pConvCheck->set_info(ss.str());
		}

	protected:
	///	adds a scaled vector to a second one
		bool VecScaleAppend(vector_type& a, vector_type& b, number s)
		{
			#ifdef UG_PARALLEL
			if(a.has_storage_type(PST_UNIQUE) && b.has_storage_type(PST_UNIQUE));
			else if(a.has_storage_type(PST_CONSISTENT) && b.has_storage_type(PST_CONSISTENT));
			else if (a.has_storage_type(PST_ADDITIVE) && b.has_storage_type(PST_ADDITIVE));
			else
			{
				a.change_storage_type(PST_ADDITIVE);
				b.change_storage_type(PST_ADDITIVE);
			}
			#endif

            for(size_t i = 0; i < a.size(); ++i)
            {
            	// todo: move VecScaleAppend to ParallelVector
            	VecScaleAdd(a[i], 1.0, a[i], s, b[i]);
            }
            return true;
		}

	///	computes the vector product
		number VecProd(vector_type& a, vector_type& b)
		{
			return a.dotprod(b);
		}

	protected:
	/// Operator that is inverted by this Inverse Operator
		ILinearOperator<vector_type,vector_type>* m_A;

	/// Preconditioner
		ILinearIterator<vector_type,vector_type>* m_pPrecond;

	/// Convergence Check
		IConvergenceCheck* m_pConvCheck;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__ */
