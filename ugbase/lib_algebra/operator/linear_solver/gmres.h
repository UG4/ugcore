/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GMRES__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GMRES__

#include <iostream>
#include <string>

#include "lib_algebra/operator/interface/operator.h"
#include "common/profiler/profiler.h"
#include "lib_algebra/operator/interface/pprocess.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

///	the GMREs method as a solver for linear operators
/**
 * This class implements the GMRES - method for the solution of linear
 * operator problems like A*x = b, where the solution x = A^{-1} b is computed.
 *
 * For detailed description of the algorithm, please refer to:
 *
 * - Barrett, Berry, Chan, Demmel, Donatom Dongarra, Eijkhout, Pozo, Romine,
 * 	 Van der Vorst, "Templates for the Solution of Linear Systems: Building
 * 	 Blocks for Iterative Methods"
 *
 * - Saad, "Iterative Methods For Sparse Linear Systems"
 *
 * \tparam 	TVector		vector type
 */
template <typename TVector>
class GMRES
	: public IPreconditionedLinearOperatorInverse<TVector>
{
	public:
	///	Vector type
		using vector_type = TVector;

	///	Base type
		using base_type = IPreconditionedLinearOperatorInverse<vector_type>;

	protected:
		using base_type::convergence_check;
		using base_type::linear_operator;
		using base_type::preconditioner;
		using base_type::write_debug;

	public:
	///	default constructor
		GMRES(size_t restart) : m_restart(restart) {};

	///	constructor setting the preconditioner and the convergence check
		GMRES( size_t restart,
		       SmartPtr<ILinearIterator<vector_type> > spPrecond,
		       SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
			: base_type(spPrecond, spConvCheck), m_restart(restart)
		{};
		~GMRES() override = default;

	///	name of solver
		const char* name() const override {return "GMRES";}

	///	returns if parallel solving is supported
		bool supports_parallel() const override {
			if(preconditioner().valid())
				return preconditioner()->supports_parallel();
			return true;
		}

	// 	Solve J(u)*x = b, such that x = J(u)^{-1} b
		bool apply_return_defect(vector_type& x, vector_type& b) override {
		//	check correct storage type in parallel
			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW("GMRES: Inadequate storage format of Vectors.");
			#endif

		//	copy rhs
			SmartPtr<vector_type> spR = b.clone();

		// 	build defect:  b := b - A*x
			linear_operator()->apply_sub(*spR, x);

		//	prepare convergence check
			prepare_conv_check();

		//	compute start defect norm
			convergence_check()->start(*spR);

		//	storage for v, h, gamma
			std::vector<SmartPtr<vector_type> > v(m_restart+1);
			std::vector<std::vector<number> > h(m_restart+1);
			for(size_t i = 0; i < h.size(); ++i) h[i].resize(m_restart+1);
			std::vector<number> gamma(m_restart+1);
			std::vector<number> c(m_restart+1);
			std::vector<number> s(m_restart+1);

		//	old norm
			number oldNorm;

		// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
			//	get storage for first vector v[0]
				if(v[0].invalid()) v[0] = x.clone_without_values();

			// 	apply v[0] = M^-1 * (b-A*x)
				if(preconditioner().valid()){
					if(!preconditioner()->apply(*v[0], *spR)){
						UG_LOG("GMRES: Cannot apply preconditioner to b-A*x0.\n");
						return false;
					}
				}
			// 	... or reuse v[0] = (b-A*x)
				else{
					SmartPtr<vector_type> tmp = v[0]; v[0] = spR; spR = tmp;
				}

			// 	make v[0] unique
				#ifdef UG_PARALLEL
				if(!v[0]->change_storage_type(PST_UNIQUE))
					UG_THROW("GMRES: Cannot convert v0 to consistent vector.");
				#endif

			//	post-process the correction
				m_corr_post_process.apply (*v[0]);

			// 	Compute norm of inital residuum:
				oldNorm = gamma[0] = v[0]->norm();

			//	normalize v[0] := v[0] / ||v[0]||
				*v[0] *= 1./gamma[0];

			//	loop gmres iterations
				size_t numIter = 0;
				for(size_t j = 0; j < m_restart; ++j)
				{
					numIter = j;

				//	get storage for v[j+1]
					if(v[j+1].invalid()) v[j+1] = x.clone_without_values();

#ifdef UG_PARALLEL
					if(!v[j]->change_storage_type(PST_CONSISTENT))
						UG_THROW("GMRES: Cannot convert v["<<j+1<<"] to consistent vector.");
#endif

				//	compute r = A*v[j]
					linear_operator()->apply(*spR, *v[j]);

				// 	apply v[j+1] = M^-1 * A * v[j]
					if(preconditioner().valid()){
						if(!preconditioner()->apply(*v[j+1], *spR)){
							UG_LOG("GMRES: Cannot apply preconditioner to A*v["<<j<<"].\n");
							return false;
						}
					}
				// 	... or reuse v[j+1] = A * v[j]
					else{
						SmartPtr<vector_type> tmp = v[j+1]; v[j+1] = spR; spR = tmp;
					}

				// 	make v[j], v[j+1] unique
					#ifdef UG_PARALLEL
					if(!v[j]->change_storage_type(PST_UNIQUE))
						UG_THROW("GMRES: Cannot convert v0 to consistent vector.");
					if(!v[j+1]->change_storage_type(PST_UNIQUE))
						UG_THROW("GMRES: Cannot convert v["<<j<<"] to consistent vector.");
					#endif

				//	post-process the correction
					m_corr_post_process.apply (*v[j+1]);

				//	loop previous steps
					for(size_t i = 0; i <= j; ++i)
					{
					//	h_ij := (r, v[j])
						h[i][j] = VecProd(*v[j+1], *v[i]);

					//	v[j+1] -= h_ij * v[i]
						VecScaleAppend(*v[j+1], *v[i], (-1)*h[i][j]);
					}

				//	compute h_{j+1,j}
					h[j+1][j] = v[j+1]->norm();

				//	update h
					for(size_t i = 0; i < j; ++i)
					{
						const number hij = h[i][j];
						const number hi1j = h[i+1][j];

						h[i][j]   =  c[i+1]*hij + s[i+1]*hi1j;
						h[i+1][j] =  s[i+1]*hij - c[i+1]*hi1j;
					}

				//	alpha := sqrt(h_jj ^2 + h_{j+1,j}^2)
					const number alpha = sqrt(h[j][j]*h[j][j] + h[j+1][j]*h[j+1][j]);

				//	update s, c
					s[j+1] = h[j+1][j] / alpha;
					c[j+1] = h[j][j]   / alpha;
					h[j][j] = alpha;

				//	compute new norm
					gamma[j+1] = s[j+1]*gamma[j];
					gamma[j] = c[j+1]*gamma[j];

					if(preconditioner().valid()) {
						UG_LOG(std::string(convergence_check()->get_offset(),' '));
						UG_LOG("% GMRES "<<std::setw(4) <<j+1<<": "
							   << gamma[j+1] << "    " << gamma[j+1] / oldNorm);
						UG_LOG(" (in Precond-Norm) \n");
						oldNorm = gamma[j+1];
					}
					else{
						convergence_check()->update_defect(gamma[j+1]);
					}

				//	normalize v[j+1]
					*v[j+1] *= 1./(h[j+1][j]);
				}

			//	compute current x
				for(size_t i = numIter; ; --i){
					for(size_t j = i+1; j <= numIter; ++j)
						gamma[i] -= h[i][j] * gamma[j];

					gamma[i] /= h[i][i];

				//	x = x + gamma[i] * v[i]
					VecScaleAppend(x, *v[i], gamma[i]);

					if(i == 0) break;
				}

			//	compute fresh defect: b := b - A*x
				*spR = b;
				linear_operator()->apply_sub(*spR, x);

				if(preconditioner().valid())
					convergence_check()->update(*spR);
			}

		//	print ending output
			return convergence_check()->post();
		}

	public:
		std::string config_string() const override {
			std::stringstream ss;
			ss << "GMRes ( restart = " << m_restart << ")\n";
			ss << base_type::config_string_preconditioner_convergence_check();
			return ss.str();
		}
		
	///	adds a post-process for the iterates
		void add_postprocess_corr (SmartPtr<IPProcessVector<vector_type> > p)
		{
			m_corr_post_process.add (p);
		}

	///	removes a post-process for the iterates
		void remove_postprocess_corr (SmartPtr<IPProcessVector<vector_type> > p)
		{
			m_corr_post_process.remove (p);
		}

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

	protected:
	///	restart parameter
		size_t m_restart;

	///	postprocessor for the correction in the iterations
		/**
		 * These postprocess operations are applied to the preconditioned
		 * defect before the orthogonalization. The goal is to prevent the
		 * useless kernel parts to prevail in the (floating point) arithmetics.
		 */
		PProcessChain<vector_type> m_corr_post_process;

	///	adds a scaled vector to a second one
		static bool VecScaleAppend(vector_type& a, vector_type& b, number s)
		{
			#ifdef UG_PARALLEL // ø todo check this logical
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

#endif