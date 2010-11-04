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

// constructors
namespace ug{

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
}

template<typename matrix_type, typename vector_type, typename densematrix_type>
void MultiEnergyProd(const matrix_type &A,
			vector_type **x,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	UG_ASSERT(n == rA.num_rows() && n == rA.num_cols(), "");
	vec_type Ai_xc;
	vector_type t;
	t.resize(x[0]->size());

	for(size_t r=0; r<n; r++)
	{
		A.apply(t, (*x[r]));
		for(size_t c=r; c<n; c++)
			rA(r, c) = VecProd((*x[c]), t);
	}

	for(size_t r=0; r<n; r++)
		for(size_t c=0; c<r; c++)
			rA(r,c) = rA(c, r);
}



template<typename TAlgebra>
class PINVIT
{
public:
	PINVIT()
	{
		m_A = NULL;
		m_B = NULL;
	}


	template<typename function_type>
	int Eigensolver(ILinearOperator<vector_type, vector_type> &A,
			ILinearOperator<vector_type, vector_type> B,
			vector_type *x,
			vector_type *corr,
			vector_type *oldcorr,
			double *lambda,
			size_t n,
			ILinearizedIteratorOperator<function_type, function_type> &precond,
			int maxIterations,
			double precision)
	{
		UG_LOG("Eigensolver\n");
		DenseMatrix<VariableArray2<double, ColMajor> > rA;
		DenseMatrix<VariableArray2<double, ColMajor> > rB;
		DenseMatrix<VariableArray2<double, ColMajor> > r_ev;
		DenseVector<VariableArray1<double> > r_lambda;

		typedef double vec_type;
		function_type defect;
		defect.clone_pattern(x[0]);

		std::vector<Vector<vec_type> *> pTestVectors;

		std::vector<vec_type> x_tmp(n);

		std::vector<double> corrnorm(n, precision*10);

		std::vector<string> testVectorDescription;

		for(int iteration=0; iteration<maxIterations; iteration++)
		{
			// 0. normalize
			for(int i=0; i<n; i++)
				x[i] *= 1/ sqrt( EnergyProd(x[i], B, x[i]));

			// 1. compute rayleigh quotients

			for(size_t i=0; i<n; i++)
				lambda[i] = EnergyProd(x[i], A, x[i]);


			// 2. before calculating new correction, save old correction
			for(size_t i=0; i<n; i++)
				swap(oldcorr[i], corr[i]);

			// 3. compute defects, apply preconditioner, compute corrections norm
			size_t nrofconverged=0;
			for(size_t i=0; i<n; i++)
			{
				/*if(corrnorm[i] < precision)
				{
					nrofconverged++;
					continue;
				}*/

				//for(size_t j=0; j<A.num_rows(); j++)
					//defect.get_vector()[j] = A[j]*x[i].get_vector() - lambda[i] * /*B*/ x[i].get_vector()[j];


				// 2. calculate residuum
				// defect = A x[i] - lambda[i] B x[i]
				MatMult(defect, 1.0, A, x[i]);
				MatMultAdd(defect, 1.0, defect, -lambda[i], B, x[i]);

				// 3. check if converged
				//corrnorm[i]= defect.two_norm();

				defect.set_storage_type(PST_ADDITIVE);
				corrnorm[i] = defect.two_norm();
				if(corrnorm[i] < precision)
				{
					nrofconverged++;
					continue;
				}


				// ?????
				precond.prepare(defect, defect, defect); // why 3 parameters? why not 2?
				precond.apply(defect, corr[i], false);
				//corrnorm[i] = corr[i].two_norm();
			}

			// output
			cout << "================================================" << endl;
			cout << "iteration " << iteration << endl;
			for(size_t i=0; i<n; i++)
				cout << i << " lambda: " << setw(14) << lambda[i] << " defect: " << setw(14) << corrnorm[i] << endl;
			cout << endl;

			if(nrofconverged==n)
			{
				UG_LOG("all eigenvectors converged" << endl);
				return true;
			}

			// 5. add Testvectors

			pTestVectors.clear();
			testVectorDescription.clear();
			for(size_t i=0; i<n; i++)
			{
				pTestVectors.push_back(&x[i].get_vector());
				testVectorDescription.push_back(string("eigenvector [") + ToString(i) + string("]") );
				if(corrnorm[i] < precision)
					continue;
				pTestVectors.push_back(&corr[i].get_vector());
				testVectorDescription.push_back(string("correction [") + ToString(i) + string("]") );
				pTestVectors.push_back(&oldcorr[i].get_vector());
				testVectorDescription.push_back(string("old correction [") + ToString(i) + string("]") );
			}

			size_t iNrOfTestVectors = pTestVectors.size();

			// 5. compute reduced Matrices rA, rB
			rA.resize(iNrOfTestVectors, iNrOfTestVectors);
			rB.resize(iNrOfTestVectors, iNrOfTestVectors);

			MultiEnergyProd(B, &pTestVectors[0], rB, iNrOfTestVectors);
			//MultiScalProd(&pTestVectors[0], rB, iNrOfTestVectors);
			// todo: remove doubled corrections by gauss-eliminating rB

			MultiEnergyProd(A, &pTestVectors[0], rA, iNrOfTestVectors);

			/*cout << "A: " << endl;
			cout << "A := matrix([";
			for(size_t r=0; r<rA.num_rows(); r++)
			{
				cout << "[";
				for(size_t c=0; c<rA.num_cols(); c++)
				{
					cout << rA(r, c);
					if(c < rA.num_cols()-1) cout << ", ";
				}
				cout << "]";
				if(r < rA.num_rows()-1) cout << ", ";
			}
			cout << "]);" << endl;

			cout << "B: " << endl;
			cout << "B := matrix([";
			for(size_t r=0; r<rB.num_rows(); r++)
			{
				cout << "[";
				for(size_t c=0; c<rB.num_cols(); c++)
				{
					cout << rB(r, c);
					if(c < rB.num_cols()-1) cout << ", ";
				}
				cout << "]";
				if(r < rB.num_rows()-1) cout << ", ";
			}
			cout << "]);" << endl;

			cout.flush();*/
			// output

			// 6. solve reduced eigenvalue problem

			r_ev.resize(iNrOfTestVectors, iNrOfTestVectors);
			r_lambda.resize(iNrOfTestVectors);

			GeneralizedEigenvalueProblem(rA, r_ev, r_lambda, rB, true);

			size_t nrzero;
			for(nrzero=0; nrzero<iNrOfTestVectors; nrzero++)
				if(r_lambda[nrzero] > 1e-15)
					break;

			if(nrzero)
			{
				cout << "Lambda < 0: " << endl;
				for(size_t i=0; i<nrzero; i++)
					cout << i << ".: " << r_lambda[i] << endl;
			}

			cout << "Lambda > 0: " << endl;
			for(size_t i=nrzero; i<r_lambda.size(); i++)
				cout << i << ".: " << r_lambda[i] << endl;

			cout << endl;


			cout << "evs: " << endl;
			for(size_t c=nrzero; c < min(nrzero+n, r_ev.num_cols()); c++)
			{
				cout << "ev [" << c << "]:" << endl;
				for(size_t r=0; r<r_ev.num_rows(); r++)
					if(dabs(r_ev(r, c)) > 1e-5 )
						cout << setw(14) << r_ev(r, c) << "   " << testVectorDescription[r] << endl;
				cout << endl << endl;
			}


			// assume r_lambda is sorted
			for(size_t i=0; i<A.num_rows(); i++)
			{

				// since x is part of the Testvektors, temporary safe result in x_tmp.
				for(size_t r=0; r<n; r++)
				{
					x_tmp[r] = 0.0;
					for(size_t c=0; c<iNrOfTestVectors; c++)
						x_tmp[r] += r_ev(c, r+nrzero) * (*pTestVectors[c])[i];
				}

				// now overwrite
				for(size_t r=0; r<n; r++)
					x[r].get_vector()[i] = x_tmp[r];

			}
		}

		UG_LOG("not converged after" << maxIterations << " steps." << endl);
		return false;
	}



private:
	ILinearOperator<vector_type,vector_type>* m_A;
	ILinearOperator<vector_type,vector_type>* m_B;

}



#endif // __H__UG__LIB_ALGEBRA__PINVIT_H__
