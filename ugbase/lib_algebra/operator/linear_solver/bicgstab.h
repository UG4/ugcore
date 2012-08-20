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

#include "lib_algebra/operator/interface/operator.h"
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
 * \tparam 	TVector		vector type
 */
template <typename TVector>
class BiCGStab
	: public IPreconditionedLinearOperatorInverse<TVector>
{
	public:
	///	Vector type
		typedef TVector vector_type;

	///	Base type
		typedef IPreconditionedLinearOperatorInverse<vector_type> base_type;

	protected:
		using base_type::convergence_check;
		using base_type::linear_operator;
		using base_type::preconditioner;
		using base_type::write_debug;

	public:
	///	default constructor
		BiCGStab() {};

	///	constructor setting the preconditioner and the convergence check
		BiCGStab( SmartPtr<ILinearIterator<vector_type> > spPrecond,
		          SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
			: base_type(spPrecond, spConvCheck)
		{};

	///	name of solver
		virtual const char* name() const {return "BiCGStab";}

	// 	Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
		//	check correct storage type in parallel
			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW("BiCGStabSolver::apply_return_defect:"
								"Inadequate storage format of Vectors.");
			#endif

		// 	build defect:  b := b - A*x
			linear_operator()->apply_sub(b, x);

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
			convergence_check()->start(b);

		//	convert b to unique (should already be unique due to norm calculation)
			#ifdef UG_PARALLEL
			if(!b.change_storage_type(PST_UNIQUE))
				UG_THROW("BiCGStab::apply_return_defect: "
								"Cannot convert b to unique vector.");
			#endif

		//	needed variables
			number rhoOld=1,rho=1,alpha=1, omega=1;

		// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
			//	check if start values have to be set
				if(convergence_check()->step() == 0 /*or restart*/)
				{
				//	if at restart recompute start defect
					if(convergence_check()->step() != 0)
					{
						convergence_check()->update(b);
						if(convergence_check()->iteration_ended()) break;
					}

				// 	reset vectors
					r = b;

				// 	make r additive unique
					#ifdef UG_PARALLEL
					if(!r.change_storage_type(PST_UNIQUE))
						UG_THROW("BiCGStab::apply_return_defect: "
										"Cannot convert r to unique vector.");
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
					UG_LOG("BiCGStab::apply_return_defect: "
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
				if(preconditioner().valid())
				{
				// 	apply q = M^-1 * p
					if(!preconditioner()->apply(q, p))
					{
						UG_LOG("ERROR in 'BiCGStab::apply_return_defect': "
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
						UG_THROW("BiCGStab::apply_return_defect: "
										"Cannot convert q to consistent vector.");
					#endif
				}

			// 	compute v := A*q
				linear_operator()->apply(v, q);

			// 	make v unique
				#ifdef UG_PARALLEL
				if(!v.change_storage_type(PST_UNIQUE))
					UG_THROW("BiCGStab::apply_return_defect: "
									"Cannot convert v to unique vector.");
				#endif

			//	alpha = (v,r)
				alpha = VecProd(v, r);

			//	check validity of alpha
				if(alpha == 0.0)
				{
					UG_LOG("ERROR in 'BiCGStab::apply_return_defect': "
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
				convergence_check()->update(s);
				if(convergence_check()->iteration_ended())
				{
				//	set output to last defect
					b = s; break;
				}

			//	apply preconditioner
				if(preconditioner().valid())
				{
				// 	apply q = M^-1 * t
					if(!preconditioner()->apply(q, s))
					{
						UG_LOG("ERROR in 'BiCGStab::apply_return_defect': "
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
						UG_THROW("BiCGStab::apply_return_defect: "
										"Cannot convert q to consistent vector.");
					#endif
				}

			// 	compute t := A*q
				linear_operator()->apply(t, q);

			// 	make t unique
				#ifdef UG_PARALLEL
				if(!t.change_storage_type(PST_UNIQUE))
					UG_THROW("BiCGStab::apply_return_defect: "
									"Cannot convert t to unique vector.");
				#endif

			// 	tt = (t,t)
				const number tt = VecProd(t, t);

			// 	omega = (s,t)
				omega = VecProd(s, t);

			//	check tt
				if(tt == 0.0)
				{
					UG_LOG("ERROR in 'BiCGStab::apply_return_defect': "
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
				convergence_check()->update(b);

			// 	remember current rho
				rhoOld = rho;
			}

		//	print ending output
			return convergence_check()->post();
		}

	protected:
	///	prepares the output of the convergence check
		void prepare_conv_check()
		{
		//	set iteration symbol and name
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');

		//	set preconditioner string
			std::stringstream ss;
			if(preconditioner().valid())
				ss <<" (Precond: "<<preconditioner()->name()<<")";
			else
				ss << " (No Preconditioner) ";
			convergence_check()->set_info(ss.str());
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
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__ */
