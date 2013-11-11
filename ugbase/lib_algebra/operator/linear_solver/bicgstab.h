/*
 * bicgstab.h
 *
 *  Created on: 05.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__

#include <iostream>
#include <string>
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
		BiCGStab() :
			m_numRestarts(0), m_minOrtho(0.0)
		{};

	///	constructor setting the preconditioner and the convergence check
		BiCGStab( SmartPtr<ILinearIterator<vector_type> > spPrecond,
		          SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
			: base_type(spPrecond, spConvCheck),
			  m_numRestarts(0), m_minOrtho(0.0)
		{};

	///	name of solver
		virtual const char* name() const {return "BiCGStab";}

	// 	Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
		//	check correct storage type in parallel
			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW("BiCGStab: Inadequate storage format of Vectors.");
			#endif

		// 	build defect:  r := b - A*x
			linear_operator()->apply_sub(b, x);
			vector_type& r = b;

		// 	create vectors
			SmartPtr<vector_type> spR = r.clone_without_values(); vector_type& r0 = *spR;
			SmartPtr<vector_type> spP = r.clone_without_values(); vector_type& p = *spP;
			SmartPtr<vector_type> spV = r.clone_without_values(); vector_type& v = *spV;
			SmartPtr<vector_type> spT = r.clone_without_values(); vector_type& t = *spT;
			SmartPtr<vector_type> spS = r.clone_without_values(); vector_type& s = *spS;
			SmartPtr<vector_type> spQ = x.clone_without_values(); vector_type& q = *spQ;

		//	prepare convergence check
			prepare_conv_check();

		//	compute start defect norm
			convergence_check()->start(r);

		//	convert b to unique (should already be unique due to norm calculation)
			#ifdef UG_PARALLEL
			if(!r.change_storage_type(PST_UNIQUE))
				UG_THROW("BiCGStab: Cannot convert b to unique vector.");
			#endif

		//	needed variables
			number rho = 1, alpha = 1, omega = 1, norm_r0 = 0.0;

		//	restart flag (set to true at first run)
			bool bRestart = true;

		// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
			//	check for restart based on fixed step number restart
				if(m_numRestarts > 0 &&
					(convergence_check()->step() % m_numRestarts == 0)){
					std::stringstream ss; ss <<
					"Restarting: at every "<<m_numRestarts<<" Iterations";
					convergence_check()->print_line(ss.str());
					bRestart = true;
				}

			//	check if start values have to be set
				if(bRestart)
				{
				// 	reset arbitrary vectors
					r0 = r;

				// 	make r additive unique
					#ifdef UG_PARALLEL
					if(!r0.change_storage_type(PST_UNIQUE))
						UG_THROW("BiCGStab: Cannot convert r to unique vector.");
					#endif

				//	set start vectors and alpha, omega:
				//	This will lead to p = r, since beta = 0.0
					p = 0.0; alpha = 0.0;
					v = 0.0; omega = 1.0;

				//	set rhoOld to 1 (rhoOld = rho, see below)
					rho = 1.0;

				//	remember start norm
					norm_r0 = convergence_check()->defect();

				//	remove restart flag
					bRestart = false;
				}

			// 	remember current rho
				const number rhoOld = rho;

			// 	Compute rho new
				rho = VecProd(r0, r);

			//	check for restart compare (r, r0) > m_minOrtho * ||r|| ||r0||
				const number norm_r = convergence_check()->defect();
				if(fabs(rho)/(norm_r * norm_r0) <= m_minOrtho){
					std::stringstream ss; ss <<
					"Restarting: Min Orthogonality "<<m_minOrtho<<" missed: "
					<<"(r,r0)="<<fabs(rho)<<", ||r||="<<norm_r<<", ||r0||= "
					<<norm_r0;
					convergence_check()->print_line(ss.str());
					bRestart = true;
				}

			//	check that rhoOld valid
				if(rhoOld == 0.0){
					UG_LOG("BiCGStab: Method breakdown with rhoOld = "<<rhoOld<<
						   ". Aborting iteration.\n");
					return false;
				}

			// 	Compute new beta
				const number beta = (rho/rhoOld) * (alpha/omega);

			//	update p = r + beta * p - beta * omega * v
				VecScaleAdd(p, 1.0, r, beta, p, -beta*omega, v);

			// 	apply q = M^-1 * p ...
				if(preconditioner().valid()){
					if(!preconditioner()->apply(q, p)){
						UG_LOG("BiCGStab: Cannot apply preconditioner. Aborting.\n");
						return false;
					}
			// 	... or copy q = p
				} else{
					q = p;

				// 	make q consistent
					#ifdef UG_PARALLEL
					if(!q.change_storage_type(PST_CONSISTENT))
						UG_THROW("BiCGStab: Cannot convert q to consistent vector.");
					#endif
				}

			// 	compute v := A*q
				linear_operator()->apply(v, q);

			// 	make v unique
				#ifdef UG_PARALLEL
				if(!v.change_storage_type(PST_UNIQUE))
					UG_THROW("BiCGStab: Cannot convert v to unique vector.");
				#endif

			//	alpha = (v,r)
				alpha = VecProd(v, r0);

			//	check validity of alpha
				if(alpha == 0.0){
					UG_LOG("BiCGStab: Method breakdown: alpha = "<<alpha<<
					       " is an invalid value. Aborting iteration.\n");
					return false;
				}

			//	alpha = rho/(v,r)
				alpha = rho/alpha;

			// 	add: x := x + alpha * q
				VecScaleAdd(x, 1.0, x, alpha, q);

			//  compute s = r - alpha*v
				VecScaleAdd(s, 1.0, r, -alpha, v);

			// 	check convergence
				convergence_check()->update(s);

			//	if finished: set output to last defect and exist loop
				if(convergence_check()->iteration_ended()){
					r = s; break;
				}

			// 	apply q = M^-1 * t ...
				if(preconditioner().valid()){
					if(!preconditioner()->apply(q, s)){
						UG_LOG("BiCGStab: Cannot apply preconditioner. Aborting.\n");
						return false;
					}
			// 	... or set q:=s
				}else{
					q = s;

				// 	make q consistent
					#ifdef UG_PARALLEL
					if(!q.change_storage_type(PST_CONSISTENT))
						UG_THROW("BiCGStab: Cannot convert q to consistent vector.");
					#endif
				}

			// 	compute t := A*q
				linear_operator()->apply(t, q);

			// 	make t unique
				#ifdef UG_PARALLEL
				if(!t.change_storage_type(PST_UNIQUE))
					UG_THROW("BiCGStab: Cannot convert t to unique vector.");
				#endif

			// 	tt = (t,t)
				const number tt = VecProd(t, t);

			// 	omega = (s,t)
				omega = VecProd(s, t);

			//	check tt
				if(tt == 0.0){
					UG_LOG("BiCGStab: Method breakdown tt = "<<tt<<" is an "
							"invalid value. Aborting iteration.\n");
					return false;
				}

			// 	omega = (s,t)/(t,t)
				omega = omega/tt;

			// 	add: x := x + omega * q
				VecScaleAdd(x, 1.0, x, omega, q);

			//  compute r = s - omega*t
				VecScaleAdd(r, 1.0, s, -omega, t);

			// 	check convergence
				convergence_check()->update(r);

			//	check values
				if(omega == 0.0){
					UG_LOG("BiCGStab: Method breakdown with omega = "<<omega<<
					       ". Aborting iteration.\n");
					return false;
				}
			}

		//	print ending output
			return convergence_check()->post();
		}


	///	sets to restart at given number of iteration steps
		void set_restart(int numRestarts) {m_numRestarts = numRestarts;}

	///	sets to restart if given orthogonality missed
		void set_min_orthogonality(number minOrtho) {m_minOrtho = minOrtho;}

	protected:
	///	prepares the output of the convergence check
		void prepare_conv_check()
		{
		//	set iteration symbol and name
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');

		//	set preconditioner string
			std::string s;
			if(preconditioner().valid())
			  s = std::string(" (Precond: ") + preconditioner()->name() + ")";
			else
				s = " (No Preconditioner) ";
			convergence_check()->set_info(s);
		}

	public:
		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "BiCGStab( restart = " << m_numRestarts << ", min_orthogonality = " << m_minOrtho << ")\n";
			ss << base_type::config_string_preconditioner_convergence_check();
			return ss.str();
		}
	protected:
	/// restarts at every numRestarts steps (numRestarts <= 0 --> never)
		int m_numRestarts;

	///	minimal value in (0,1) accepted for Orthoginality before restart
		number m_minOrtho;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__BICGSTAB__ */
