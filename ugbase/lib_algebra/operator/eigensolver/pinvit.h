/**
 * \file densematrix_impl.h
 *
 * \author Martin Rupp
 *
 * \date 01.11.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__LIB_ALGEBRA__PINVIT_H__
#define __H__UG__LIB_ALGEBRA__PINVIT_H__

#include "common/util/string_util.h"

// constructors
namespace ug{


/*template<typename mat_type, typename vec_type, typename densematrix_type>
void MultiEnergyProd(const SparseMatrix<mat_type> &A,
			Vector<vec_type> **x,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vec_type Ai_xc;
	rA = 0.0;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		for(size_t c=0; c<n; c++)
		{
			// Ai_xc = A[i] * x[c].
			Ai_xc = 0.0;
			A.mat_mult_add_row(i, Ai_xc, 1.0, (*x[c]));
			for(size_t r=0; r<n; r++)
				rA(c, r) += VecDot((*x[r])[i], Ai_xc);
		}
	}
}*/

template<typename vector_type, typename densematrix_type>
void MultiScalProd(vector_type **px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	UG_ASSERT(0, "");
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	for(size_t r=0; r<n; r++)
	{
		for(size_t c=r; c<n; c++)
		{
			rA(r, c) = px[c]->dotprod(*px[r]);
			UG_LOG("MultiScalProd : " << rA(r, c) << "\n");
		}
	}

	for(size_t r=0; r<n; r++)
		for(size_t c=0; c<r; c++)
			rA(r,c) = rA(c, r);
}


template<typename matrix_type, typename vector_type>
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2, vector_type &tmp)
{

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
	v2.change_storage_type(PST_CONSISTENT);
#endif
	A.apply(tmp, v2);
#ifdef UG_PARALLEL
	tmp.change_storage_type(PST_CONSISTENT);
#endif
	double a = v1.dotprod(tmp);
	//UG_LOG("EnergyProd " << a << "\n");

	return a;
}

template<typename matrix_type, typename vector_type>
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2)
{
	vector_type t;
	CloneVector(t, v1);
	return EnergyProd(v1, A, v2, t);
}

template<typename matrix_type, typename vector_type>
double EnergyNorm(vector_type &x, matrix_type &A, vector_type &tmp)
{
	return sqrt( EnergyProd(x, A, x, tmp) );
}

template<typename matrix_type, typename vector_type>
double EnergyNorm(vector_type &x, matrix_type &A)
{
	vector_type tmp;
	CloneVector(tmp, x);
	return sqrt( EnergyProd(x, A, x, tmp) );
}


template<typename matrix_type, typename vector_type, typename densematrix_type>
void MultiEnergyProd(matrix_type &A,
			vector_type **px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
#endif
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vector_type t;
	CloneVector(t, *px[0]);


	for(size_t r=0; r<n; r++)
	{
		// todo: why is SparseMatrix<T>::apply not const ?!?
#ifdef UG_PARALLEL
		px[r]->change_storage_type(PST_CONSISTENT);
#endif
		A.apply(t, (*px[r]));
#ifdef UG_PARALLEL
		t.change_storage_type(PST_CONSISTENT);
#endif
		for(size_t c=r; c<n; c++)
		{
			//rA(r, c) = VecProd((*px[c]), t);
			rA(r, c) = px[c]->dotprod(t);
			//UG_LOG("MultiEnergyProd : (" << r << ", " << c << ") = " << rA(r, c) << "\n");
		}
	}


	for(size_t r=0; r<n; r++)
		for(size_t c=0; c<r; c++)
			rA(r,c) = rA(c, r);
}


template<typename tvector>
void PrintStorageType(const tvector &v)
{
#ifdef UG_PARALLEL
	if(v.has_storage_type(PST_UNDEFINED))
		UG_LOG("PST_UNDEFINED ");
	if(v.has_storage_type(PST_CONSISTENT))
		UG_LOG("PST_CONSISTENT ");
	if(v.has_storage_type(PST_ADDITIVE))
		UG_LOG("PST_ADDITIVE ");
	if(v.has_storage_type(PST_UNIQUE))
		UG_LOG("PST_UNIQUE ");
#else
	UG_LOG("serial ");
#endif
}


template<typename matrix_type>
void PrintMatrix(const matrix_type &mat, const char *name)
{
	UG_LOG(name << ":\n" << name << " := matrix([\n");
	for(size_t r=0; r<mat.num_rows(); r++)
	{
		UG_LOG("[");
		for(size_t c=0; c<mat.num_cols(); c++)
		{
			UG_LOG(mat(r, c));
			if(c < mat.num_cols()-1) UG_LOG(",\t");
		}
		UG_LOG("]\n");
	}
	UG_LOG("]);\n");

}

template<typename matrix_type>
void PrintMaple(const matrix_type &mat, const char *name)
{
	UG_LOG(name << ":\n" << name << " := matrix([");
	for(size_t r=0; r<mat.num_rows(); r++)
	{
		UG_LOG("[");
		for(size_t c=0; c<mat.num_cols(); c++)
		{
			UG_LOG(mat(r, c));
			if(c < mat.num_cols()-1) UG_LOG(", ");
		}
		UG_LOG("]");
		if(r < mat.num_rows()-1) UG_LOG(", ");
	}
	UG_LOG("]);\n");

}

template<typename T>
void MemSwap(T &a, T &b)
{
	char c[sizeof(T)];
	memcpy(c, &a, sizeof(T));
	memcpy(&a, &b, sizeof(T));
	memcpy(&b, c, sizeof(T));
}

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



	std::vector<vector_type*> px;
	SmartPtr<ILinearIterator<vector_type> > m_spPrecond;

	size_t m_maxIterations;
	double m_dPrecision;
	size_t m_iPINVIT;

public:
	PINVIT()
	{
		m_pA = NULL;
		m_pB = NULL;
		m_iPINVIT = 3;
	}

	/**
	 * adds a vector which should be used for eigenvalue computation
	 * @param vec
	 */
	void add_vector(vector_type &vec)
	{
		px.push_back(&vec);
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

	void set_linear_operator_A(SmartPtr<ILinearOperator<vector_type> > A)
	{
		m_pA = A;
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

	/**
	 * perform the calculation
	 * @return
	 */
	int apply()
	{
		UG_LOG("Eigensolver\n");
		DenseMatrix<VariableArray2<double> > rA;
		DenseMatrix<VariableArray2<double> > rB;
		DenseMatrix<VariableArray2<double> > r_ev;
		DenseVector<VariableArray1<double> > r_lambda;
		std::vector<double> lambda;

		typedef typename vector_type::value_type value_type;
		vector_type defect;
		CloneVector(defect, *px[0]);
		size_t n = px.size();

		size_t size = px[0]->size();
		/*
		ParallelMatrix<SparseMatrix<double> > B;
		B.resize(size, size);
		for(size_t i=0; i<size; i++)
			B(i,i) = 0.00390625;
		B.set_storage_type(PST_ADDITIVE);*/

		vector_type tmp;
		CloneVector(tmp, *px[0]);
		std::vector<vector_type> vCorr;

		std::vector<vector_type> vOldX;


		lambda.resize(n);
		vCorr.resize(n);
		vOldX.resize(n);
		for(size_t i=0; i<n; i++)
		{
			UG_ASSERT(px[0]->size() == px[i]->size(), "all vectors must have same size");
			CloneVector(vCorr[i], *px[0]);
			CloneVector(vOldX[i], *px[0]);

			//PrintStorageType(*px[i]);
			//PrintStorageType(vCorr[i]);
			//PrintStorageType(vOldX[i]);
		}

		std::vector<vector_type *> pTestVectors;

		std::vector<double> vDefectNorm(n, m_dPrecision*10);
		std::vector<double> oldXnorm(n);

		std::vector<std::string> vTestVectorDescription;

		m_spPrecond->init(m_pA);
#ifdef UG_PARALLEL
		pcl::ProcessCommunicator pc;
#endif

		for(size_t iteration=0; iteration<m_maxIterations; iteration++)
		{

			UG_LOG("iteration " << iteration << "\n");

			// 0. normalize
			normalize_approximations();


			// 1. before calculating new correction, save old correction
			save_old_approximations(vOldX);

			//  2. compute rayleigh quotient, residuum, apply preconditioner, compute corrections norm
			size_t nrofconverged=0;


			write_debug(iteration);
			for(size_t i=0; i<n; i++)
			{
				compute_rayleigh_and_new_correction(*px[i], lambda[i], defect, vDefectNorm[i], vCorr[i]);
				if(vDefectNorm[i] < m_dPrecision)
					nrofconverged++;
			}


			// output
			print_eigenvalues_and_defect(iteration, vDefectNorm, oldXnorm, lambda);

			if(nrofconverged==n)
			{
				UG_LOG("all eigenvectors converged\n");
				return true;
			}

			// 5. add Testvectors
			//UG_LOG("5. add Testvectors\n");

			get_testvectors(iteration, vCorr, vOldX, pTestVectors, vTestVectorDescription, vDefectNorm);


			/*for(size_t i=0; i<vTestVectorDescription.size(); i++)
			{	UG_LOG(vTestVectorDescription[i] << "\n");	} */

			// 5. compute reduced Matrices rA, rB

			get_projected_eigenvalue_problem(rA, rB, pTestVectors, vTestVectorDescription);


			// 6. solve reduced eigenvalue problem
			size_t iNrOfTestVectors = pTestVectors.size();
			r_ev.resize(iNrOfTestVectors, iNrOfTestVectors);
			r_lambda.resize(iNrOfTestVectors);

			// solve rA x = lambda rB, --> r_ev, r_lambda
			GeneralizedEigenvalueProblem(rA, r_ev, r_lambda, rB, true);

			size_t nrzero;
			for(nrzero=0; nrzero<iNrOfTestVectors; nrzero++)
				if(r_lambda[nrzero] > 1e-15)
					break;

			if(nrzero)
			{
				UG_LOG("Lambda < 0: \n");
				for(size_t i=0; i<nrzero; i++)
				{
					UG_LOG(i << ".: " << r_lambda[i] << "\n");

				}
			}
			UG_LOG("Lambda > 0: \n");
			for(size_t i=nrzero; i<r_lambda.size(); i++)
				UG_LOG(i << ".: " << r_lambda[i] << "\n");
			UG_LOG("\n");

			/*UG_LOG("evs: \n");
			for(size_t c=nrzero; c < std::min(nrzero+n, r_ev.num_cols()); c++)
			{
				UG_LOG("ev [" << c << "]:\n");
				for(size_t r=0; r<r_ev.num_rows(); r++)
					if(dabs(r_ev(r, c)) > 1e-9 )
						UG_LOG(std::setw(14) << r_ev(r, c) << "   " << vTestVectorDescription[r] << "\n");
				UG_LOG("\n\n");
			}*/

#ifdef UG_PARALLEL
			for(size_t i=0; i<iNrOfTestVectors; i++)
				pTestVectors[i]->change_storage_type(PST_UNIQUE);
			for(size_t i=0; i<n; i++)
				px[i]->change_storage_type(PST_UNIQUE);
#endif
			// assume r_lambda is sorted
			std::vector<typename vector_type::value_type> x_tmp(n);
			for(size_t i=0; i<size; i++)
			{
				// since x is part of the Testvectors, temporary safe result in x_tmp.
				for(size_t r=0; r<n; r++)
				{
					x_tmp[r] = 0.0;
					for(size_t c=0; c<iNrOfTestVectors; c++)
						x_tmp[r] += r_ev(c, r+nrzero) * (*pTestVectors[c])[i];
				}

				// now overwrite
				for(size_t r=0; r<n; r++)
					(*px[r])[i] = x_tmp[r];

			}
		}

		UG_LOG("not converged after" << m_maxIterations << " steps.\n");
		return false;
	}

private:
	void write_debug(int iteration)
	{
		for(size_t i=0; i<px.size(); i++)
		{
			string name = "pinvit_it_" + ToString(iteration) + "_ev_" + ToString(i);
			write_debug(*px[i], name.c_str());
		}
	}

	double B_norm(vector_type &x)
	{
		if(m_pB != NULL)
			return EnergyNorm(x, *m_pB);
		else
			return x.norm();
	}

	void normalize_approximations()
	{
		for(size_t i=0; i< px.size(); i++)
			(*px[i]) *= 1/ (B_norm(*px[i]));
	}

	/**
	 * For a given eigenvalue approximation, computes the
	 * rayleigh quotient, the defect, the norm of the defect, and the correction calculated by the preconditioner
	 * @param[in]  x			current normalized eigenvalue approximation (<x,x> = 1)
	 * @param[out] lambda		lambda = <x, Ax> / <x, x>
	 * @param[out] defect		defect = lambda x - Ax
	 * @param[out] vDefectNorm 	vDefectNorm = | defect |_2
	 * @param[out] vCorr		P defect
	 */
	void compute_rayleigh_and_new_correction(vector_type &x, double &lambda, vector_type &defect, double &vDefectNorm, vector_type &vCorr)
	{
// a. compute rayleigh quotients
		// lambda = <x, Ax>/<x,x>
		// todo: replace with MatMult
//				UG_LOG("m_pA has storage type "); PrintStorageType(*m_pA); UG_LOG(", and vector px[" << i << "] has storage type"); PrintStorageType(*px[i]); UG_LOG("\n");
		// px can be set to unique because of norm

#ifdef UG_PARALLEL
		x.change_storage_type(PST_CONSISTENT);
		defect.set_storage_type(PST_ADDITIVE);
#endif
		m_pA->apply(defect, x);


#ifdef UG_PARALLEL
		defect.change_storage_type(PST_UNIQUE);
		x.change_storage_type(PST_UNIQUE);
#endif
		lambda = x.dotprod(defect); // / <px[i], px[i]> = 1.0.
		//UG_LOG("lambda[" << i << "] = " << lambda << "\n");

// b. calculate residuum
		// defect = A px[i] - lambda[i] B px[i]
		if(m_pB)
		{
			// todo: replace with MatMultAdd
			//MatMultAddDirect(defect, 1.0, defect, -lambda[i], *m_pB, *px[i]);
#ifdef UG_PARALLEL
			defect.change_storage_type(PST_ADDITIVE);
			x.change_storage_type(PST_CONSISTENT);
#endif
			MatMultAddDirect(defect, 1.0, defect, -lambda, *m_pB, x);
		}
		else
			VecScaleAdd(defect, 1.0, defect, -lambda, x);

// c. check if converged


#ifdef UG_PARALLEL
		defect.change_storage_type(PST_UNIQUE);
#endif
		vDefectNorm = defect.norm();
#ifdef UG_PARALLEL
		defect.change_storage_type(PST_UNIQUE);
#endif
		if(vDefectNorm < 1e-12)
			return;

// d. apply preconditioner
		m_spPrecond->apply(vCorr, defect);
		vCorr *= 1/ B_norm(vCorr);
#ifdef UG_PARALLEL
		vCorr.change_storage_type(PST_UNIQUE);
#endif

	}

	/**
	 * prints the current eigenvalues and convergence status
	 * @param[in] 		iteration		iteration number
	 * @param[in] 		vDefectNorm		vector of defect norms
	 * @param[in,out] 	vOldDefectNorm	vector of defect norms from previous iteration
	 * @param[in] 		vLambda			vector of eigenvalue approximations
	 */
	void print_eigenvalues_and_defect(int iteration, const std::vector<double> &vDefectNorm,
			std::vector<double> &vOldDefectNorm, const std::vector<double> &vLambda)
	{
		UG_LOG("================================================\n");
		UG_LOG("iteration " << iteration << "\n");

		for(size_t i=0; i<vLambda.size(); i++)
		{
			UG_LOG(i << " lambda: " << std::setw(14) << vLambda[i] << " defect: " <<
					std::setw(14) << vDefectNorm[i]);
			if(iteration != 0) { UG_LOG(" reduction: " << std::setw(14) << vDefectNorm[i]/vOldDefectNorm[i]); }
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
	void get_testvectors(int iteration, std::vector<vector_type> &vCorr, std::vector<vector_type> &vOldX,
			std::vector<vector_type *> &pTestVectors, std::vector<std::string> &vTestVectorDescription,
			const std::vector<double> &vDefectNorm)
	{
		pTestVectors.clear();
		vTestVectorDescription.clear();
		for(size_t i=0; i < px.size(); i++)
		{
			if(m_iPINVIT == 1)
			{
				VecScaleAdd(vCorr[i], -1.0, vCorr[i], 1.0, *px[i]);
				vTestVectorDescription.push_back(std::string("ev - vCorr [") + ToString(i) + std::string("]") );
			}
			else
			{
				pTestVectors.push_back(px[i]);
				vTestVectorDescription.push_back(std::string("eigenvector [") + ToString(i) + std::string("]") );

				if(vDefectNorm[i] < m_dPrecision)
						continue;
				pTestVectors.push_back(&vCorr[i]);
				vTestVectorDescription.push_back(std::string("correction [") + ToString(i) + std::string("]") );

				if(m_iPINVIT >= 3)
				{
					pTestVectors.push_back(&vOldX[i]);
					vTestVectorDescription.push_back(std::string("old correction [") + ToString(i) + std::string("]") );

					if(iteration == 0)
					{
						for(size_t j=0; j<px[i]->size(); j++)
							vOldX[i][j] = (*px[i])[j] * urand(-1.0, 1.0);
					}
				}
			}
		}
	}

	/**
	 * save current eigenvector approximation into vector old
	 * @param[out] old vector to save approximations to
	 */
	void save_old_approximations( std::vector<vector_type> &old)
	{
		for(size_t i=0; i<px.size(); i++)
			old[i] = *px[i];
	}

	/**
	 * Calculates a maximal set of rows which are linear independent
	 * @param[in]  mat					the input matrix
	 * @param[out] bLinearIndependent	output vector (true if linear independent)
	 */
	void get_linear_independent_rows(DenseMatrix<VariableArray2<double> > mat, std::vector<bool> &bLinearIndependent)
	{
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
			if(mat(i,i) < 1e-8) bLinearIndependent[i] = false;
			else bLinearIndependent[i] = true;
		}
	}

	/**
	 * remove all entries with vbUse[i]==false from vector i
	 * @param[in, out] 	v 		vector to contain result
	 * @param[in] 		vbUse	if vbUse[i] is true, add it to new vector
	 */
	template<typename T>
	void remove_unused(std::vector<T> &v, const std::vector<bool> vbUse)
	{
		std::vector<T> tmp = v;
		v.clear();
		for(size_t i=0; i<tmp.size(); i++)
			if(vbUse[i])
				v.push_back(tmp[i]);

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
			DenseMatrix<VariableArray2<double> > &rB, std::vector<vector_type *> &pTestVectors,
			std::vector<std::string> &vTestVectorDescription)
	{
		// 1. calculate W as a subset of the testvectors so that those are linear B-independent

		size_t iNrOfTestVectors = pTestVectors.size();
		rA.resize(iNrOfTestVectors, iNrOfTestVectors);
		rB.resize(iNrOfTestVectors, iNrOfTestVectors);

		if(m_pB)
			MultiEnergyProd(*m_pB, &pTestVectors[0], rB, iNrOfTestVectors);
		else
			MultiScalProd(&pTestVectors[0], rB, iNrOfTestVectors);


		// Remove linear depended vectors
		std::vector<bool> bUse;
		get_linear_independent_rows(rB, bUse);


		// save used testvectors
		remove_unused(pTestVectors, bUse);
		remove_unused(vTestVectorDescription, bUse);

		iNrOfTestVectors = pTestVectors.size();

		// 2. & 3. compute reduced Matrices rA, rB
		rA.resize(iNrOfTestVectors, iNrOfTestVectors);
		rB.resize(iNrOfTestVectors, iNrOfTestVectors);

		if(m_pB)
			MultiEnergyProd(*m_pB, &pTestVectors[0], rB, iNrOfTestVectors);
		else
			MultiScalProd(&pTestVectors[0], rB, iNrOfTestVectors);

		MultiEnergyProd(*m_pA, &pTestVectors[0], rA, iNrOfTestVectors);

		PrintMaple(rA, "rA");
		PrintMaple(rB, "rB");
	}


};

} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__PINVIT_H__
