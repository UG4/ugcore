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

template<typename mat_type, typename vec_type, typename densematrix_type>
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
	typename vector_type::value_type Ai_xc;
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
// 	Algebra type
	typedef TAlgebra algebra_type;

// 	Vector type
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

private:
	//ILinearOperator<vector_type,vector_type>* m_pA;
	//ILinearOperator<vector_type,vector_type>* m_pB;
	matrix_type *m_pA;
	matrix_type *m_pB;

	std::vector<vector_type*> px;
	ILinearIterator<vector_type, vector_type> &m_pPrecond;

	size_t m_maxIterations;
	double m_precision;

public:
	PINVIT()
	{
		m_A = NULL;
		m_B = NULL;
	}

	void add_vector(vector_type &vec)
	{
		px.push_back(&vec);
	}

	void set_preconditioner(ILinearIterator<vector_type, vector_type> &precond)
	{
		m_pPrecond = &precond;
	}

	bool set_linear_operator_A(ILinearOperator<vector_type, vector_type> &A)
	{
		MatrixOperator<vector_type, vector_type, matrix_type>* Op =
					dynamic_cast<MatrixOperator<vector_type, vector_type, matrix_type>*>(&J);

	//	Check that matrix if of correct type
		if(Op == NULL)
		{
			UG_LOG("ERROR in '" << name() << "::init': Passed Operator is not based on matrix.\n"
					"This Preconditioner can only handle matrix-based operators. Aborting.\n");
			return false;
		}

		m_pA = &Op->get_matrix();
	}

	bool set_linear_operator_B(ILinearOperator<vector_type, vector_type> &B)
	{
		MatrixOperator<vector_type, vector_type, matrix_type>* Op =
							dynamic_cast<MatrixOperator<vector_type, vector_type, matrix_type>*>(&J);

		//	Check that matrix if of correct type
		if(Op == NULL)
		{
			UG_LOG("ERROR in '" << name() << "::init': Passed Operator is not based on matrix.\n"
					"This Preconditioner can only handle matrix-based operators. Aborting.\n");
			return false;
		}

		m_pB = &Op->get_matrix();
	}

	void set_max_iterations(size_t maxIterations)
	{
		m_maxIterations = maxIterations;
	}

	void set_precision(double precision)
	{
		m_precision = precision;
	}


	int apply()
	{
		UG_LOG("Eigensolver\n");
		DenseMatrix<VariableArray2<double> > rA;
		DenseMatrix<VariableArray2<double> > rB;
		DenseMatrix<VariableArray2<double> > r_ev;
		DenseVector<VariableArray1<double> > r_lambda;

		typedef typename vector_type::value_type value_type;
		vector_type defect;
		size_t size = px[0]->size();
		size_t n = px.size();
		defect.create(size);

		std::vector<vector_type> corr;
		std::vector<vector_type> oldcorr;
		for(size_t i=0; i<n; i++)
		{
			UG_ASSERT(size == px[i]->size(), "all vectors must have same size");
			corr[i].create(size);
			oldcorr[i].create(size);
		}

		std::vector<vector_type *> pTestVectors;

		std::vector<double> corrnorm(n, m_precision*10);

		std::vector<string> testVectorDescription;

		m_pPrecond->init(*m_pA);

		for(int iteration=0; iteration<m_maxIterations; iteration++)
		{
			// 0. normalize
			if(m_B)
			{
				UG_ASSERT(0, "implement");
				//for(int i=0; i<n; i++)
					//(*x[i]) *= 1/ sqrt( EnergyProd(*x[i], m_B, *x[i]));
			}
			else
			{
				for(int i=0; i<n; i++)
					(*x[i]) *= 1/ x[i]->two_norm();
			}



			// 1. before calculating new correction, save old correction
			for(size_t i=0; i<n; i++)
				swap(oldcorr[i], corr[i]);

			//  2. compute rayleigh quotient, residuum, apply preconditioner, compute corrections norm
			size_t nrofconverged=0;
			for(size_t i=0; i<n; i++)
			{
				// 2a. compute rayleigh quotients
				MatMult(defect, 1.0, A, *x[i]);
				lambda[i] = VecProd(*x[i], defect);

				// 2b. calculate residuum
				// defect = A x[i] - lambda[i] B x[i]
				if(m_B)
				{
					UG_ASSERT(0, "implement");
					//MatMultAdd(defect, 1.0, defect, -lambda[i], B, *x[i]);
				}
				else
					VecScaleAdd(defect, 1.0, defect, -lambda[i], *x[i]);

				// 2c. check if converged
				//corrnorm[i]= defect.two_norm();

				defect.set_storage_type(PST_ADDITIVE);
				corrnorm[i] = defect.two_norm();
				if(corrnorm[i] < m_precision)
				{
					nrofconverged++;
					continue;
				}

				// 2d. apply preconditioner
				n_pPrecond->apply(corr[i], defect);
			}

			// output
			UG_LOG("================================================" << endl);
			UG_LOG("iteration " << iteration << endl);

			for(size_t i=0; i<n; i++)
				UG_LOG(i << " lambda: " << setw(14) << lambda[i] << " defect: " << setw(14) << corrnorm[i] << endl);
			UG_LOG(endl);

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
				pTestVectors.push_back(x[i]);
				testVectorDescription.push_back(string("eigenvector [") + ToString(i) + string("]") );
				if(corrnorm[i] < precision)
					continue;
				pTestVectors.push_back(&corr[i]);
				testVectorDescription.push_back(string("correction [") + ToString(i) + string("]") );

				//pTestVectors.push_back(&oldcorr[i]);
				//testVectorDescription.push_back(string("old correction [") + ToString(i) + string("]") );
			}

			size_t iNrOfTestVectors = pTestVectors.size();

			// 5. compute reduced Matrices rA, rB
			rA.resize(iNrOfTestVectors, iNrOfTestVectors);
			rB.resize(iNrOfTestVectors, iNrOfTestVectors);

			if(pB)
				MultiEnergyProd(*pB, &pTestVectors[0], rB, iNrOfTestVectors);
			else
				MultiScalProd(&pTestVectors[0], rB, iNrOfTestVectors);
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
					(*x[r])[i] = x_tmp[r];

			}
		}

		UG_LOG("not converged after" << maxIterations << " steps." << endl);
		return false;
	}



};

} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__PINVIT_H__
