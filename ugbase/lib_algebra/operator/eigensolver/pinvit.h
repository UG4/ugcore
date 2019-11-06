/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__LIB_ALGEBRA__PINVIT_H__
#define __H__UG__LIB_ALGEBRA__PINVIT_H__

#include <complex>
#include <vector>
#include <string>
#include "common/util/string_util.h"
#include "common/util/sort_util.h"
#include "additional_math.h"

#include "lib_algebra/small_algebra/lapack/eigenvalue2.h"

#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/debug_writer.h"
#include "common/progress.h"

#include "lib_algebra/algebra_common/vector_util.h"

#include "smart_ptr_vector.h"

// constructors
namespace ug{


/**
 * PINVIT Eigensolver
 *
 * This Eigensolver solves problems of the form
 * Ax = lambda B x
 * For sparse matrices A and B, and we want to find the smallest (in terms of abs(lambda) ) e.g. 1-100 solutions (eigenvalues) of the problem.
 * For this we need a preconditioner, calculating c = Pd (e.g. Jacobi, Gauss-Seidel, Geometric/Algebraic Multigrid, ILU).
 *
 * This implements the PINVIT(s) methods, with s=2 = LOBPCG
 * see Martin Rupp - Berechnung der Resonanzschwingungen einer Gitarrendecke (Diploma thesis)
 * and Andrew Knyazew, Toward the optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned
 * 	Conjugate Gradient Method. http://epubs.siam.org/doi/pdf/10.1137/S1064827500366124
 *
 * iPINVIT=1  -> Preconditioned Inverse Block Iteration [Neymeyr]
 * iPINVIT=2  -> Preconditioned Block Gradient Method
 * iPINVIT=3  -> LOBPCG (locally optimal block preconditioned gradient) [Knyazew]
 * iPINVIT>=4 -> gerneralized methods.
 *
 * example:
 * \code
 * EigenSolver eig;
 * eig:set_linear_operator_A(A)
 * eig:set_linear_operator_B(B)
 * eig:set_max_iterations(50)
 * eig:set_precision(evPrec)
 * eig:set_preconditioner(gmg)
 * eig:set_max_iterations(30)
 * eig:set_pinvit(2)
 * eig:set_debug(dbgWriter)
 * ev = {}
 * for i=1,nev do
 *   print("adding ev "..i)
 * 	 ev[i] = GridFunction(approxSpace)
 * 	 ev[i]:set_random(-1.0, 1.0)
 * 	 domainDisc:adjust_solution(ev[i])
 * 	 eig:add_vector(ev[i])
 * end
 * \endcode
 *
 * If you want to use a solver as preconditioner, use
 * precond = OperatorInverseIterator(linSolver)
 */
template<typename TAlgebra>
class PINVIT
		: public DebugWritingObject<TAlgebra>
{
public:
// 	Algebra type
	typedef TAlgebra algebra_type;

// 	Vector type
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename vector_type::value_type vector_value_type;

///	Base type
	typedef DebugWritingObject<TAlgebra> base_type;


private:
	using base_type::set_debug;
	using base_type::debug_writer;
	using base_type::write_debug;

	//ILinearOperator<vector_type,vector_type>* m_pA;
	//ILinearOperator<vector_type,vector_type>* m_pB;
	SmartPtr<ILinearOperator<vector_type> > m_pA;

	typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;
	matrix_operator_type *m_pB;

	double m_dMinimumDefectToCalcCorrection;

	SmartPtrVector<vector_type> px;
	SmartPtr<ILinearIterator<vector_type> > m_spPrecond;

	size_t m_maxIterations;
	double m_dPrecision;
	bool m_bRelativePrecision;
	size_t m_iPINVIT;
	std::vector<double> lambda;
	std::vector<bool> vbDirichlet;

	double m_linearDependentEps;

	bool m_bLaplacian;

	std::vector<double> freq_change;
	double m_dDensity_linear_elasticity;
	bool m_bPrintFrequencies_linear_elasticity;
	bool m_bAbortOnFreqConverged_linear_elasticity;
	double m_dFreqPrecision;

	bool m_bPrintEigenvaluesAndDefect;
	bool m_bPrintProjectedEigenvalues;
	bool m_bPrintProjectedEigenvectors;
	bool m_bPrintProjectedEigenproblem;
	bool m_bPrintUsedTestvectors;
	bool m_bUseAdditionalCorrections;
	bool m_bDebugCalcProjectedEigenvalues;
	size_t m_additionalEigenvectorsToKeep;
	size_t m_currentAdditionalEigenvectors;
	size_t m_currentAdditionalCorrections;

public:
	std::string tostring()
	{
		return config_string();
	}
	virtual std::string config_string() const
	{
		std::stringstream ss;
		ss << "PINVIT Eigensolver by Martin Rupp / G-CSC 2013-2015." <<
				"\n MaxIterations = " << m_maxIterations <<
				"\n Precision = " << m_dPrecision << (m_bRelativePrecision ? " (relative) " : " (absolute) ") <<
				"\n MinimumDefectToCalcCorrection = " << m_dMinimumDefectToCalcCorrection <<
				"\n Number of EV = " << px.size() <<
				"\n PINVIT = " << m_iPINVIT << " = ";
		if(m_iPINVIT == 1)
			 ss << "Preconditioned Inverse Block Iteration [Neymeyr]";
		else if(m_iPINVIT==2) ss << "Preconditioned Block Gradient Method";
		else if(m_iPINVIT==3) ss << "LOBPCG (locally optimal block preconditioned gradient) [Knyazew]";
		else if(m_iPINVIT>=4) ss << "gerneralized method PINVIT-method looking back " << m_iPINVIT-1 << " ev approximations.";
		ss << "\n Preconditioner: ";
		if(m_spPrecond.valid()) ss << ConfigShift(m_spPrecond->config_string()) << "\n";
		else ss << "  NOT SET!\n";

		if(m_additionalEigenvectorsToKeep>0)
		{
			ss << "\n\tAdditionaly storing " << m_additionalEigenvectorsToKeep << " eigenvectors";
		}
		if(m_bUseAdditionalCorrections)
			ss << "\n\tUsing additional Corrections.";
		if(m_bLaplacian)
			ss << "\n\tUsing LAPLACIAN!";

		ss << "\n";
		return ss.str();

	}

	PINVIT()
	{
		m_pA = SPNULL;
		m_pB = NULL;
		m_iPINVIT = 3;
		m_dMinimumDefectToCalcCorrection = 1e-8;
		m_dPrecision = 1e-8;

		m_bPrintEigenvaluesAndDefect = true;
		m_bPrintProjectedEigenvalues = false;
		m_bPrintProjectedEigenvectors = false;
		m_bPrintProjectedEigenproblem = false;
		m_bPrintUsedTestvectors = false;
		m_bUseAdditionalCorrections = false;
		m_bDebugCalcProjectedEigenvalues = false;

		m_additionalEigenvectorsToKeep = 0;
		m_currentAdditionalEigenvectors = 0;
		m_linearDependentEps = 1e-6;
		m_bLaplacian = false;
		m_bStoreDefects = false;
		m_bStoreLambdas = false;
		m_bRelativePrecision = false;

		m_dDensity_linear_elasticity = 0.0;
		m_bPrintFrequencies_linear_elasticity = false;
		m_bAbortOnFreqConverged_linear_elasticity = false;
		m_dFreqPrecision = 0.0;
	}

	void set_abort_on_frequencies_converged_linear_elasticity(double freq_precision, double density)
	{
		m_dDensity_linear_elasticity = density;
		m_bPrintFrequencies_linear_elasticity = true;
		m_bAbortOnFreqConverged_linear_elasticity = true;
		m_dFreqPrecision = freq_precision;
	}

	void set_linear_dependent_eps(double eps)
	{
		m_linearDependentEps = eps;
	}
	/**
	 * adds a vector which should be used for eigenvalue computation
	 * @param vec
	 */
	void add_vector(SmartPtr<vector_type> vec)
	{
		px.push_back(vec);
	}

	/**
	 * set the preconditioner (or Linear Iterator)
	 * e.g. Gauss-Seidel, AMG/GMG, Jacobi, ILU, ...
	 * @param precond
	 */
	void set_preconditioner(SmartPtr<ILinearIterator<vector_type> > precond)
	{
		m_spPrecond = precond;
	}

	void set_linear_operator_A(SmartPtr<ILinearOperator<vector_type> > loA)
	{
		m_pA = loA;
		// get dirichlet nodes
		SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
				m_pA.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();
		matrix_type &A = pOp->get_matrix();
		vbDirichlet.resize(A.num_rows());

		if(m_bLaplacian)
		{
			bool bFirst=false;
			for(size_t i=0; i<A.num_rows(); i++)
			{
				//vbDirichlet[i] = A.is_isolated(i);
				if(bFirst && A.is_isolated(i))
				{
					bFirst = false;
					continue;
				}

				for(typename matrix_type::row_iterator it = A.begin_row(i); it != A.end_row(i); ++it)
				{
					it.value() = -1.0;
				}
				A(i, i) = A.num_connections(i)*1.0;

			}

			for(size_t i=0; i<A.num_rows(); i++)
				vbDirichlet[i] = false;
		}
		else
		{
			for(size_t i=0; i<A.num_rows(); i++)
				vbDirichlet[i] = A.is_isolated(i);
		}

	}

	void set_linear_operator_B(matrix_operator_type &B)
	{
		m_pB = &B;
	}

	void set_max_iterations(size_t maxIterations)
	{
		m_maxIterations = maxIterations;
	}

	void set_precision(double precision)
	{
		m_dPrecision = precision;
		m_dMinimumDefectToCalcCorrection = precision;
		m_bRelativePrecision = false;
	}

	void set_additional_eigenvectors_to_keep(size_t i)
	{
		m_additionalEigenvectorsToKeep = i;
	}

	void set_debug_calc_projected_eigenvalues(bool b)
	{
		m_bDebugCalcProjectedEigenvalues = b;
	}

	void set_use_additional_corrections(bool b)
	{
		m_bUseAdditionalCorrections = b;
	}

	/**
	 * iPINVIT=1  -> Preconditioned Inverse Block Iteration [Neymeyr]
	 * iPINVIT=2  -> Preconditioned Block Gradient Method
	 * iPINVIT=3  -> LOBPCG (locally optimal block preconditioned gradient) [Knyazew]
	 * iPINVIT>=4 -> gerneralized methods.
	 * @param iPINVIT
	 */
	void set_pinvit(size_t iPINVIT)
	{
		m_iPINVIT = iPINVIT;
		UG_ASSERT(iPINVIT >=1 && iPINVIT <= 3, "i has to be >= 1 and <= 3, but is " << iPINVIT);
	}

	size_t num_eigenvalues()
	{
		return px.size();
	}

	double get_eigenvalue(size_t i)
	{
		if(lambda.size() != px.size()) return 0.0;
		return lambda[i];
	}

	SmartPtr<vector_type> get_eigenvector(size_t i)
	{
		return px[i];
	}

	void set_print_eigenvalues_and_defect(bool b)
	{
		m_bPrintEigenvaluesAndDefect = b;
	}

	void set_print_projected_eigenvectors(bool b)
	{
		m_bPrintProjectedEigenvectors = b;
	}
	void set_print_projected_eigenvalues(bool b)
	{
		m_bPrintProjectedEigenvalues = b;
	}

	void set_print_projected_eigenproblem(bool b)
	{
		m_bPrintProjectedEigenproblem = b;
	}

	void set_print_used_testvectors(bool b)
	{
		m_bPrintUsedTestvectors = b;
	}

	void print_projected_eigenvectors(DenseMatrix<VariableArray2<double> > &r_ev, std::vector<std::complex<double> > &r_lambda, std::vector<std::string> &vTestVectorDescription)
	{
		UG_LOG("evs: \n");
		for(size_t c=0; c < r_ev.num_cols(); c++)
		{
			UG_LOG("eigenvalue [" << c << "] =  " << r_lambda[c] << ", vector:\n");
			std::vector<double> tmpEV(r_ev.num_rows());
			for(size_t r=0; r<r_ev.num_rows(); r++)
				tmpEV[r] = r_ev(r, c);
			std::vector<size_t> s = GetSortedIndices(tmpEV, absCompare);


			for(size_t i=0; i<r_ev.num_rows(); i++)
			{
				//if(r_ev.num_rows() > 3 && abs(r_ev(s[i], c))/abs(r_ev(s[r_ev.num_rows()-1], c)) < 0.01) continue;
				size_t j = s[i];
				UG_LOG(std::setw(14) << r_ev(j, c) << "   " << vTestVectorDescription[j] << "\n");
			}
			UG_LOG("\n\n");
		}


	}

	SmartPtr<vector_type> create_approximation_vector()
	{
		SmartPtr<vector_type> t = px(0).clone_without_values();
		t->set(0.0);
		return t;
	}

	void set_store_defects(bool b)
	{
		m_bStoreDefects = b;
	}

	void set_store_lambdas(bool b)
	{
		m_bStoreLambdas = b;
	}

	bool m_bStoreDefects;
	bool m_bStoreLambdas;
	size_t m_iteration;
	std::vector<std::vector<double> > m_defects;
	std::vector<std::vector<double> > m_lambdas;

	size_t get_iterations()
	{
		return m_iteration;
	}

	double get_defect(size_t nev, size_t it)
	{
		UG_COND_THROW(m_bStoreDefects == false, "Not storing defects. Use set_store_defects(true).");
		UG_COND_THROW(nev > m_defects.size(), nev << " > #EV" << m_defects.size());
		UG_COND_THROW(it > m_defects[nev].size(), "Iteration " << it << " not available.");
		return m_defects[nev][it];
	}

	double get_lambda(size_t nev, size_t it)
	{
		UG_COND_THROW(m_bStoreLambdas == false, "Not storing lambdas. Use set_store_lambdas(true).");
		UG_COND_THROW(nev > m_lambdas.size(), nev << " > #EV" << m_lambdas.size());
		UG_COND_THROW(it > m_lambdas[nev].size(), "Iteration " << it << " not available.");

		return m_lambdas[nev][it];
	}

	/**
	 * perform the calculation
	 * @return
	 */
	int apply()
	{
		PINVIT_PROFILE_FUNC();
		UG_LOG(config_string() << "\nInitializing... ");


		size_t nEigenvalues = px.size();
		DenseMatrix<VariableArray2<double> > rA;
		DenseMatrix<VariableArray2<double> > rB;
		DenseMatrix<VariableArray2<double> > r_ev;
		std::vector<std::complex<double> > r_lambda;
		std::vector<std::string> vCorrectionName;

		SmartPtrVector<vector_type> pTestVectors;
		SmartPtrVector<vector_type> vAdditional;
		SmartPtrVector<vector_type> vCorr;
		SmartPtrVector<vector_type> vOldX;

		std::vector<double> vDefectNorm(nEigenvalues, m_dPrecision*10);
		std::vector<double> oldXnorm(nEigenvalues);
		std::vector<bool> bConverged(nEigenvalues, false);

		std::vector<std::string> vTestVectorDescription;

		SmartPtr<vector_type> spDefect;

		try{
			spDefect = create_approximation_vector();

			m_currentAdditionalEigenvectors = 0;

			//size_t size = px[0]->size();
			/*
			ParallelMatrix<SparseMatrix<double> > B;
			B.resize(size, size);
			for(size_t i=0; i<size; i++)
				B(i,i) = 0.00390625;
			B.set_storage_type(PST_ADDITIVE);*/

			for(size_t i=0; i<m_additionalEigenvectorsToKeep; i++)
				vAdditional.push_back(create_approximation_vector());

			lambda.resize(nEigenvalues);
			freq_change.resize(nEigenvalues);
			for(size_t i=0; i<nEigenvalues; i++)
			{
				vCorr.push_back(create_approximation_vector() );
				if(m_iPINVIT >= 3)
					vOldX.push_back(create_approximation_vector() );
			}

			vCorrectionName.resize(nEigenvalues);
		}UG_CATCH_THROW("PINVIT::apply init vectors failed.");

		UG_LOG("done.\n");
		UG_LOG("PINVIT: Initializing preconditioner... ");

		try{
			m_spPrecond->init(m_pA);
		}UG_CATCH_THROW("PINVIT::apply init preconditioner.");

		UG_LOG("done.\n");

#ifdef UG_PARALLEL
		pcl::ProcessCommunicator pc;
#endif

		m_iteration=0;
		m_defects.clear();
		if(m_bStoreDefects)
			m_defects.resize(nEigenvalues);
		m_lambdas.clear();
		if(m_bStoreLambdas)
			m_lambdas.resize(nEigenvalues);

		for(m_iteration=0; m_iteration<m_maxIterations; m_iteration++)
		{
			try{

	//			UG_LOG("normalize...\n");
				// 0. normalize
				normalize_approximations();

				//  2. compute rayleigh quotient, residuum, apply preconditioner, compute corrections norm
				size_t nrofconverged=0;


				size_t currentAdditionalCorrections = 0;
				size_t numCorrections=0;

				PROGRESS_START(prog, nEigenvalues, "PINVIT: compute corrections");
				for(size_t i=0; i<nEigenvalues; i++)
				{
					PROGRESS_UPDATE(prog, i);
	//				UG_LOG("compute rayleigh " << i << "...\n");

					double freq_old = 0;

					if(m_bPrintFrequencies_linear_elasticity){
						freq_old = 1.0/(2.0*PI)*sqrt(lambda[i]/m_dDensity_linear_elasticity);
					}

					compute_rayleigh_and_defect(px(i), lambda[i], *spDefect, vDefectNorm[i]);

					if(m_bPrintFrequencies_linear_elasticity){
						freq_change[i] = freq_old-1.0/(2.0*PI)*sqrt(lambda[i]/m_dDensity_linear_elasticity);
					}

					if((m_iteration != 0 && m_bRelativePrecision && vDefectNorm[i]/lambda[i] < m_dPrecision)
							|| (!m_bRelativePrecision && vDefectNorm[i] < m_dPrecision))
						bConverged[i] = true;
					else
						bConverged[i] = false;

					if(bConverged[i])
					{
						nrofconverged++;
						if(!m_bRelativePrecision &&
								m_bUseAdditionalCorrections && currentAdditionalCorrections < m_currentAdditionalEigenvectors)
						{
							double additionalLambda, additionalDefect;
							compute_rayleigh_and_defect(vAdditional(currentAdditionalCorrections),
									additionalLambda, *spDefect, additionalDefect);
							currentAdditionalCorrections++;
							if(additionalDefect > m_dPrecision && additionalDefect < 1e5)
							{
								calculate_correction(*vCorr[numCorrections++], *spDefect);
								vCorrectionName[numCorrections-1] = std::string("additional correction ") + ToString(currentAdditionalCorrections-1);
							}
						}
					}
					else
					{
	//					UG_LOG("Calculating correction from defect " << i << "\n");
						vCorrectionName[numCorrections] = std::string("correction ") + ToString(i);
						calculate_correction(*vCorr[numCorrections++], *spDefect);
					}

					if(m_bStoreDefects)
						m_defects[i].push_back(vDefectNorm[i]);
					if(m_bStoreLambdas)
						m_lambdas[i].push_back(lambda[i]);
				}
				PROGRESS_FINISH(prog);




				// output
				if(m_bPrintEigenvaluesAndDefect)
					print_eigenvalues_and_defect(m_iteration, vDefectNorm, oldXnorm, lambda, bConverged);

				if(nrofconverged==nEigenvalues)
				{
					UG_LOG("all eigenvectors converged\n");
					return true;
				}

				if(m_bAbortOnFreqConverged_linear_elasticity){
					bool s = true;
					for(unsigned k = 0; k < freq_change.size(); ++k){
						if(abs(freq_change[k]) >= m_dFreqPrecision){
							s = false;
							break;
						}
					}

					if(s){
						UG_LOG("all eigenvectors converged according to frequency change\n");
						return true;
					}
				}

				// 5. add Testvectors
				//UG_LOG("5. add Testvectors\nEigenvalues");

				get_testvectors(m_iteration, vCorr, vCorrectionName, numCorrections, vOldX, vAdditional, pTestVectors, vTestVectorDescription, bConverged);

				for(size_t i=0; i<nEigenvalues; i++)
				{
					write_debug(m_iteration, i, px(i), *spDefect, vCorr(i), bConverged[i]);
					if(vOldX.size() > i) write_debug_old(m_iteration, i, vOldX(i));
				}

				/*for(size_t i=0; i<vTestVectorDescription.size(); i++)
				{	UG_LOG(vTestVectorDescription[i] << "\nEigenvalues");	} */

				// 5. compute reduced Matrices rA, rB

				if(pTestVectors.size() < px.size()+1)
				{
					UG_LOG("ERROR: Projected Spaces has Dimension " << pTestVectors.size() << ", but we need at least " << px.size()+1 << " to make progress (for " << px.size() << " eigenvalues)\n");
					return false;
				}

				get_projected_eigenvalue_problem(rA, rB, pTestVectors, vTestVectorDescription);

				// 6. solve reduced eigenvalue problem
				size_t iNrOfTestVectors = pTestVectors.size();
				r_ev.resize(iNrOfTestVectors, iNrOfTestVectors);
				r_lambda.resize(iNrOfTestVectors);

				// solve rA x = lambda rB, --> r_ev, r_lambda
				GeneralizedEigenvalueProblemComplex(rA, r_ev, r_lambda, rB, true);

				if(m_bPrintProjectedEigenvalues)
					print_projected_eigenvectors(r_ev, r_lambda, vTestVectorDescription);


				if(m_bDebugCalcProjectedEigenvalues)
					debug_calc_projected_eigenvalues(r_ev, pTestVectors, m_iteration, false);
				set_new_approximations_and_save_old(r_ev, pTestVectors, vOldX, vAdditional);

				assert_real_positive(r_lambda);
			}UG_CATCH_THROW("PINVIT::apply iteration " << m_iteration << " failed.");
		}

		UG_LOG("not converged after" << m_maxIterations << " steps.\n");
		return false;
	}

	// only use this type of vector creation if you don't use the preconditioner
	void init_nonprecond_vector(vector_type &t)
	{
		CloneVector(t, px(0));
		t.set(0.0);
	}
	/**
	 * debug_calc_projected_eigenvalues
	 * @param r_ev
	 * @param pTestVectors
	 * @param iteration
	 * @param bWrite
	 */
	void debug_calc_projected_eigenvalues(DenseMatrix<VariableArray2<double> > &r_ev,
			SmartPtrVector<vector_type> &pTestVectors, int iteration, bool bWrite)
	{
		PROFILE_FUNC_GROUP("debug");
		try{
			if(debug_writer() == SPNULL) return;

			vector_type t, defect;
			init_nonprecond_vector(t);
			init_nonprecond_vector(defect);

	#ifdef UG_PARALLEL
			for(size_t c=0; c<pTestVectors.size(); c++)
				pTestVectors[c]->change_storage_type(PST_CONSISTENT);
	#endif

			for(size_t r=0; r<pTestVectors.size(); r++)
			{
	#ifdef UG_PARALLEL
				t.set_storage_type(PST_CONSISTENT);
				defect.set_storage_type(PST_ADDITIVE);
	#endif
				for(size_t i=0; i<t.size(); i++)
				{
					t[i] = 0.0;
					if(!vbDirichlet[i])
					{
						for(size_t c=0; c<pTestVectors.size(); c++)
							t[i] += r_ev(c, r) * (*pTestVectors[c])[i];
					}
				}


				t *= 1.0/B_norm(t);

	#ifdef UG_PARALLEL
				t.set_storage_type(PST_CONSISTENT);
				defect.set_storage_type(PST_ADDITIVE);
	#endif
				m_pA->apply(defect, t);

				double lambda = t.dotprod(defect);

				if(m_pB)
				{
	#ifdef UG_PARALLEL
					defect.change_storage_type(PST_ADDITIVE);
					t.change_storage_type(PST_CONSISTENT);
	#endif
					MatMultAdd(defect, 1.0, defect, -lambda, *m_pB, t);
				}
				else
					VecScaleAdd(defect, 1.0, defect, -lambda, t);



				UG_LOG(r << " lambda: " << std::setw(14) << lambda << " defect: " << std::setw(14) << defect.norm() << "\n");

				if(bWrite)
					write_debug(t, ("pinvit_it_" + ToString(iteration) + "_pev_" + ToString(r)).c_str());
			}
		}UG_CATCH_THROW("PINVIT::debug_calc_projected_eigenvalues failed.");
	}


	void set_new_approximations_and_save_old(DenseMatrix<VariableArray2<double> > &r_ev, SmartPtrVector<vector_type> &pTestVectors, SmartPtrVector<vector_type> &vOldX,
			SmartPtrVector<vector_type> &vAdditional)
	{
		PINVIT_PROFILE_FUNC();
		try{
	#ifdef UG_PARALLEL
			for(size_t i=0; i<pTestVectors.size(); i++)
				pTestVectors[i]->change_storage_type(PST_CONSISTENT);
			for(size_t i=0; i<px.size(); i++)
				px[i]->change_storage_type(PST_CONSISTENT);

			for(size_t i=0; i<vOldX.size(); i++)
				vOldX[i]->set_storage_type(PST_CONSISTENT);

			for(size_t i=0; i<vAdditional.size(); i++)
				vAdditional[i]->set_storage_type(PST_CONSISTENT);
	#endif
			// assume r_lambda is sorted

			size_t size = px[0]->size();
			THROW_IF_NOT_EQUAL(pTestVectors.size(), r_ev.num_rows());

			size_t numCalc = std::min(px.size()+m_additionalEigenvectorsToKeep, r_ev.num_cols()); //px.size()
			m_currentAdditionalEigenvectors = numCalc-px.size();
	//		UG_LOG("current aditional eigenvectors: " << m_currentAdditionalEigenvectors << "\n");

			std::vector<vector_value_type> x_tmp(numCalc);
			for(size_t i=0; i<size; i++)
			{
				if(vbDirichlet[i]) continue;

				// since vx can be part of the Testvectors, temporary safe result in x_tmp.
				for(size_t r=0; r < numCalc; r++)
				{
					x_tmp[r] = 0.0;
					for(size_t c=0; c<pTestVectors.size(); c++)
						x_tmp[r] += r_ev(c, r) * (*pTestVectors[c])[i];
				}

				// now overwrite
				for(size_t r=0; r<px.size(); r++)
				{
					// save old
					/*for(int k=0; k<m_nrOfOld; k++)
						vOldX[k][r][i] = vOldX[k][r][i];
					vOldX[0][r][i] = (px[r])[i];*/
					if(vOldX.size()>0)
						vOldX(r) [i] = px(r) [i];
					// store new
					px(r) [i] = x_tmp[r];
				}

				// store additional
				for(size_t r=px.size(); r<numCalc; r++)
					vAdditional(r-px.size()) [i] = x_tmp[r];

			}
		}UG_CATCH_THROW("PINVIT::set_new_approximations_and_save_old failed.");
	}

	void assert_real_positive(const std::vector<std::complex<double> > &r_lambda)
	{
		PINVIT_PROFILE_FUNC();
		try{
			//#define UG_ASSERT2(cond, txt) {if(!cond) { UG_LOG("Warning: " << txt << "\n"); }}

			for(size_t i=0; i<r_lambda.size(); i++) // px.size() or r_lambda.size()
			{
				if(r_lambda[i].imag() >= 1e-8)
					{UG_LOG("WARNING: eigenvalue " << i << " is imaginary (" << r_lambda[i] << ")\n"); }
				if(r_lambda[i].real() <= 1e-8)
					{UG_LOG("WARNING: eigenvalues " << i << "<= 0\n");}
			}

			if(m_bPrintProjectedEigenvalues && m_bPrintProjectedEigenvectors == false)
			{
				for(size_t i=0; i<r_lambda.size(); i++)
					UG_LOG(i << ".: " << r_lambda[i] << "\n");
				UG_LOG("\n");
			}

			for(size_t i=0; i<r_lambda.size(); i++) // px.size() or r_lambda.size()
			{
				UG_ASSERT(r_lambda[i].imag() < 1e-4 && r_lambda[i].imag() > -1e-4, "eigenvalue " << i << " = " << r_lambda[i] << " is imaginary (" << r_lambda[i] << ")\n");
				UG_ASSERT(r_lambda[i].real() > -1e-4, "eigenvalues " << i << " = " << r_lambda[i] << " < 0\n");
			}
		}UG_CATCH_THROW("PINVIT::assert_real_positive failed.");

	}

	void get_max_deflection_of_a_mode(vector_type& maxDeflect, vector_type& statSol,
			const vector_type& eigenVec, const matrix_type& massMat)
	{
		PINVIT_PROFILE_FUNC();
#ifdef UG_PARALLEL
		statSol.set_storage_type(PST_CONSISTENT);
		maxDeflect.set_storage_type(PST_CONSISTENT);
#endif

		SmartPtr<vector_type> Bv = statSol.clone_without_values();
		SmartPtr<vector_type> vwBv = statSol.clone_without_values();

		//	compute Bv = massMat * eigenVec;
		massMat.axpy(*Bv, 0.0, eigenVec, 1.0, eigenVec);

		//	vwBv = eigenVec * < statSol, Bv >
		number scalarProd = 0.0;
		for(size_t i = 0; i < statSol.size(); i++)
			scalarProd += VecProd((*Bv)[i], statSol[i]);

		for(size_t i = 0; i < statSol.size(); i++)
		{
			(*vwBv)[i] = eigenVec[i] * (-1.0) * scalarProd;
			//	maxDeflect += vwBv
			maxDeflect[i] = maxDeflect[i] + (*vwBv)[i];
		}

	}

private:
	void write_debug(int iteration, int i, vector_type &x, vector_type &defect, vector_type &corr, bool bConverged)
	{
		PROFILE_FUNC_GROUP("debug");
		write_debug(x, ("pinvit_it_" + ToString(iteration) + "_ev_" + ToString(i)).c_str());
		write_debug(defect, ("pinvit_it_" + ToString(iteration) + "_defect_" + ToString(i)).c_str());
		if(!bConverged)
			write_debug(corr, ("pinvit_it_" + ToString(iteration) + "_corr_" + ToString(i)).c_str());
	}

	void write_debug_old(int iteration, int i, vector_type &oldX)
	{
		write_debug(oldX, ("pinvit_it_" + ToString(iteration) + "_old_" + ToString(i)).c_str());
	}

	double B_norm(vector_type &x)
	{
		PINVIT_PROFILE_FUNC();
		if(m_pB != NULL)
			return EnergyNorm(x, *m_pB);
		else
			return x.norm();
	}

	void normalize_approximations()
	{
		PINVIT_PROFILE_FUNC();
		for(size_t i=0; i< px.size(); i++)
			px(i) *= 1/ B_norm(px(i));
	}

	void calculate_correction(vector_type &corr, vector_type &defect)
	{
		PINVIT_PROFILE_FUNC();
		try{
	#ifdef UG_PARALLEL
			defect.change_storage_type(PST_ADDITIVE);
			corr.set_storage_type(PST_CONSISTENT);
	#endif
			corr *= 0.0;
			// d. apply preconditioner
			if(1)
			{
				m_spPrecond->apply(corr, defect);
			}
			else
			{
				m_spPrecond->apply_update_defect(corr, defect);
			}
			corr *= 1/ B_norm(corr);
	#ifdef UG_PARALLEL
			corr.change_storage_type(PST_UNIQUE);
	#endif
		}UG_CATCH_THROW("PINVIT::calculate_correction failed.");
	}
	/**
	 * For a given eigenvalue approximation, computes the
	 * rayleigh quotient, the defect, the norm of the defect, and the correction calculated by the preconditioner
	 * @param[in]  x			current normalized eigenvalue approximation (<x,x> = 1)
	 * @param[out] lambda		lambda = <x, Ax> / <x, Bx>
	 * @param[out] defect		defect = lambda x - Ax
	 * @param[out] vDefectNorm 	vDefectNorm = | defect |_2
	 * @param[out] vCorr		P defect
	 */
	void compute_rayleigh_and_defect(vector_type &x, double &lambda, vector_type &defect, double &defectNorm)
	{
		PINVIT_PROFILE_FUNC();
		try{
	// a. compute rayleigh quotients

	#ifdef UG_PARALLEL
			x.change_storage_type(PST_CONSISTENT);
			defect.set_storage_type(PST_ADDITIVE);
	#endif
			m_pA->apply(defect, x);


	#ifdef UG_PARALLEL
			defect.change_storage_type(PST_UNIQUE);
			x.change_storage_type(PST_UNIQUE);
	#endif
			lambda = x.dotprod(defect); // / <px[i], Bpx[i]> = 1.0.
			//UG_LOG("lambda[" << i << "] = " << lambda << "\n");

	// b. calculate residuum
			// defect = A px[i] - lambda[i] B px[i]
			if(m_pB)
			{
	#ifdef UG_PARALLEL
				defect.change_storage_type(PST_ADDITIVE);
				x.change_storage_type(PST_CONSISTENT);
	#endif
				MatMultAdd(defect, 1.0, defect, -lambda, *m_pB, x);

			}
			else
				VecScaleAdd(defect, 1.0, defect, -lambda, x);

	// c. check if converged


	#ifdef UG_PARALLEL
			defect.change_storage_type(PST_UNIQUE);
	#endif
			defectNorm = defect.norm();
		}UG_CATCH_THROW("PINVIT::compute_rayleigh_and_defect failed.");
	}

	/**
	 * prints the current eigenvalues and convergence status
	 * @param[in] 		iteration		iteration number
	 * @param[in] 		vDefectNorm		vector of defect norms
	 * @param[in,out] 	vOldDefectNorm	vector of defect norms from previous iteration
	 * @param[in] 		vLambda			vector of eigenvalue approximations
	 */
	void print_eigenvalues_and_defect(int iteration, const std::vector<double> &vDefectNorm,
			std::vector<double> &vOldDefectNorm, const std::vector<double> &vLambda, std::vector<bool> &bConverged)
	{
		PINVIT_PROFILE_FUNC();
		UG_LOG("=====================================================================================\n");
		UG_LOG("iteration " << iteration << "\n");

		for(size_t i=0; i<vLambda.size(); i++)
		{
			UG_LOG(i << " lambda: " << std::setw(20) << vLambda[i]);
			if(m_bPrintFrequencies_linear_elasticity){
				UG_LOG(" freq: " << std::setw(20) << 1.0/(2.0*PI)*sqrt(vLambda[i]/m_dDensity_linear_elasticity) );
				UG_LOG(" change: " << std::setw(20) << freq_change[i]);
			}

			UG_LOG(" defect: " << std::setw(20) << vDefectNorm[i]);
			if(iteration != 0) { UG_LOG(" reduction: " << std::setw(14) << vDefectNorm[i]/vOldDefectNorm[i]); }
			if(m_bPrintFrequencies_linear_elasticity && abs(freq_change[i]) < m_dFreqPrecision) { UG_LOG(" (freq. converged) "); }
			if(bConverged[i]) { UG_LOG(" (converged)"); }
			UG_LOG("\n");
			vOldDefectNorm[i] = vDefectNorm[i];
		}
		UG_LOG("\n");
	}

	/**
	 * depending on the PINVIT-method, this function calculates the used testvectors
	 * - for PINVIT(1), projected space is  L^k = span_i < c^k_i - x^{k}_i>,
	 *   that is (current eigenvalue - correction)
	 * - PINVIT(s) for s>=2:
	 *   L^k = span_i < x^{k-s+2}_i , .. x^{k}_i, c^k_i>
	 *   that is the space spanned by the current eigenvalue, its correction, and the s-2 previous eigenvalue approximations
	 *
	 *  if an eigenvalue is converged, we don't calculate a correction for this and previous approximations will be
	 *  not too much different, so we also won't add previous approximations
	 * @param[in]  iteration				iteration number
	 * @param[in]  vCorr					correction for eigenvector approximation
	 * @param[in]  vOldX					previous eigenvector approximations
	 * @param[out] pTestVectors				vector in which we store the used test vectors for the projected eigenvalue problem
	 * @param[out] vTestVectorDescription	description of the vectors (ev, corr or oldEv)
	 * @param[in]  vDefectNorm				norm of the defects
	 */
	void get_testvectors(int iteration, SmartPtrVector<vector_type> &vCorr,
			std::vector<std::string> &vCorrectionName,
			size_t numCorrections, SmartPtrVector<vector_type> &vOldX, SmartPtrVector<vector_type> &vAdditional,
			SmartPtrVector<vector_type> &pTestVectors, std::vector<std::string> &vTestVectorDescription,
			std::vector<bool> &bConverged)
	{
		try{
			PINVIT_PROFILE_FUNC();
			pTestVectors.clear();
			vTestVectorDescription.clear();
			if(m_iPINVIT == 1)
			{
				UG_ASSERT(0, "todo");
				for(size_t i=0; i < px.size(); i++)
				{
					if(!bConverged[i])
					{
						VecScaleAdd(vCorr(i), -1.0, vCorr(i), 1.0, px(i));
						vTestVectorDescription.push_back(std::string("ev - vCorr [") + ToString(i) + std::string("]") );
					}
				}
			}
			else
			{
				for(size_t i=0; i < px.size(); i++)
				{
					if(IsFiniteAndNotTooBig(px(i)))
					{
						pTestVectors.push_back(px[i]);
						vTestVectorDescription.push_back(std::string("eigenvector [") + ToString(i) + std::string("]") );

						if(!bConverged[i] && m_iPINVIT >= 3)
						{
							if(iteration != 0)
							{
								UG_ASSERT(vOldX.size() > i, "?");
								pTestVectors.push_back(vOldX[i]);
								vTestVectorDescription.push_back(std::string("old eigenvalue[") + ToString(i) + std::string("]") );
							}

							/*if(iteration == 0)
							{
								for(size_t j=0; j<px[i]->size(); j++)
									vOldX[i][j] = (px[i])[j] * urand(-1.0, 1.0);
							}*/
						}
					}
					else
					{	UG_LOG("WARNING ev [" << i << "] not normal or too big\n");	}
				}

				for(size_t i=0; i<numCorrections; i++)
				{
					if(IsFiniteAndNotTooBig(vCorr(i)))
					{
						pTestVectors.push_back(vCorr[i]);
						vTestVectorDescription.push_back(vCorrectionName[i]);
					}
					else
					{
						UG_LOG("WARNING correction [" << i << "] not normal or too big\n");
					}
				}

			}
			for(size_t i=0; i<m_currentAdditionalEigenvectors; i++)
			{
				if(IsFiniteAndNotTooBig(vAdditional(i)))
				{
					pTestVectors.push_back(vAdditional[i]);
					vTestVectorDescription.push_back(std::string("additional [") + ToString(i) + std::string("]") );
				}
				else
				{	UG_LOG("WARNING additional [" << i << "] not normal or too big\n");	}
			}
		}UG_CATCH_THROW("PINVIT::get_testvectors failed.");
	}


	/**
	 * Calculates a maximal set of rows which are linear independent
	 * @param[in]  mat					the input matrix
	 * @param[out] bLinearIndependent	output vector (true if linear independent)
	 */
	void get_linear_independent_rows(DenseMatrix<VariableArray2<double> > mat, std::vector<bool> &bLinearIndependent)
	{
		PINVIT_PROFILE_FUNC();
		// Remove linear depended vectors
		bLinearIndependent.resize(mat.num_rows(), true);
		for(size_t i=0; i<mat.num_rows(); i++)
		{
			for(size_t j=0; j<i; j++)
			{
				if(!bLinearIndependent[i]) continue;
				double val = mat(i, j)/mat(j,j);
				mat(i,j) = 0;
				for(size_t k=j+1; k<mat.num_rows(); k++)
					mat(i,k) -= val*mat(j, k);
			}
			if(mat(i,i) < m_linearDependentEps) bLinearIndependent[i] = false;
			else bLinearIndependent[i] = true;
		}
	}

	/**
	 * remove all entries with vbUse[i]==false from vector i
	 * @param[in, out] 	v 		vector to contain result
	 * @param[in] 		vbUse	if vbUse[i] is true, add it to new vector
	 */
	template<typename T>
	void remove_unused(T &v, const std::vector<bool> vbUse)
	{
		PINVIT_PROFILE_FUNC();
		T tmp;
		tmp = v;
		v.clear();
		for(size_t i=0; i<tmp.size(); i++)
			if(vbUse[i])
				v.push_back(tmp[i]);

	}

	void print_used_testvectors(std::vector<std::string> &vTestVectorDescription, std::vector<bool> bUse)
	{
		UG_LOG("used testvectors:\n");
		for(size_t i=0; i<vTestVectorDescription.size(); i++)
			if(bUse[i]) { UG_LOG(vTestVectorDescription[i] << "\n"); }
		UG_LOG("unused testvectors:\n");
		for(size_t i=0; i<vTestVectorDescription.size(); i++)
			if(!bUse[i]) { UG_LOG(vTestVectorDescription[i] << "\n"); }
	}

	/**
	 * Calculate projected eigenvalue problem on space spanned by testvectors in pTestVectors
	 * 1. calculate W as a subset of the testvectors so that those are linear B-independent
	 * 2. rA = W^T A W
	 * 3. rB = W^T B W
	 * @param[out] rA							reduced eigenvalue problem matrix A
	 * @param[out] rB							reduced eigenvalue problem matrix B
	 * @param[in, out] pTestVectors				the testvectors, out: the used testvectors
	 * @param[in, out] vTestVectorDescription 	their description
	 */
	void get_projected_eigenvalue_problem(DenseMatrix<VariableArray2<double> > &rA,
			DenseMatrix<VariableArray2<double> > &rB, SmartPtrVector<vector_type> &pTestVectors,
			std::vector<std::string> &vTestVectorDescription)
	{
		try{
			PINVIT_PROFILE_FUNC();
			// 1. calculate W as a subset of the testvectors so that those are linear B-independent

			size_t iNrOfTestVectors = pTestVectors.size();
			rA.resize(iNrOfTestVectors, iNrOfTestVectors);
			rB.resize(iNrOfTestVectors, iNrOfTestVectors);

			if(m_pB)
				MultiEnergyProd(*m_pB, pTestVectors, rB, iNrOfTestVectors);
			else
				MultiScalProd(pTestVectors, rB, iNrOfTestVectors);


			// Remove linear depended vectors
			std::vector<bool> bUse;
			get_linear_independent_rows(rB, bUse);


			if(m_bPrintUsedTestvectors)
				print_used_testvectors(vTestVectorDescription, bUse);

			// save used testvectors
			remove_unused(pTestVectors, bUse);
			remove_unused(vTestVectorDescription, bUse);

			iNrOfTestVectors = pTestVectors.size();

			// 2. & 3. compute reduced Matrices rA, rB
			rA.resize(iNrOfTestVectors, iNrOfTestVectors);
			rB.resize(iNrOfTestVectors, iNrOfTestVectors);

	//		UG_LOG("iNrOfTestVectors = " << iNrOfTestVectors << "\n");

			if(m_pB)
				MultiEnergyProd(*m_pB, pTestVectors, rB, iNrOfTestVectors);
			else
				MultiScalProd(pTestVectors, rB, iNrOfTestVectors);

			MultiEnergyProd(*m_pA, pTestVectors, rA, iNrOfTestVectors);

			if(m_bPrintProjectedEigenproblem)
			{
				PrintMaple(rA, "rA");
				PrintMaple(rB, "rB");
			}
		}UG_CATCH_THROW("PINVIT::get_projected_eigenvalue_problem failed.");

	}
public:
	void set_laplacian(bool b)
	{
		m_bLaplacian = b;
	}

	void set_relative_precision(double precision)
	{
		m_dPrecision = precision;
		m_dMinimumDefectToCalcCorrection = precision;
		m_bRelativePrecision = true;
	}

};

} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__PINVIT_H__
