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

template<typename T>
inline std::string ToString(const T &t)
{
	std::stringstream out;
	out << t;
	return out.str();
}

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
double EnergyProd(vector_type &v1, matrix_type &A, vector_type &v2)
{
	pcl::ProcessCommunicator pc;

	vector_type t;
	CloneVector(t, v1);
#ifdef UG_PARALLEL
	v2.change_storage_type(PST_CONSISTENT);
#endif
	A.apply(t, v2);
#ifdef UG_PARALLEL
	t.change_storage_type(PST_CONSISTENT);
#endif
	double a = v1.dotprod(t);
	//UG_LOG("EnergyProd " << a << "\n");

	return a;
}

template<typename matrix_type, typename vector_type, typename densematrix_type>
void MultiEnergyProd(matrix_type &A,
			vector_type **px,
			DenseMatrix<densematrix_type> &rA, size_t n)
{
	pcl::ProcessCommunicator pc;
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

	void add_vector(vector_type &vec)
	{
		px.push_back(&vec);
	}

	void set_preconditioner(SmartPtr<ILinearIterator<vector_type> > precond)
	{
		m_spPrecond = precond;
	}

	bool set_linear_operator_A(SmartPtr<ILinearOperator<vector_type> > A)
	{
		m_pA = A;
		return true;
	}

	bool set_linear_operator_B(matrix_operator_type &B)
	{
		m_pB = &B;
		return true;
	}

	void set_max_iterations(size_t maxIterations)
	{
		m_maxIterations = maxIterations;
	}

	void set_precision(double precision)
	{
		m_dPrecision = precision;
	}

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

		matrix_operator_type &B = *m_pB;

		vector_type tmp;
		CloneVector(tmp, *px[0]);
		std::vector<vector_type> corr;
		std::vector<vector_type> oldcorr;
		lambda.resize(n);
		corr.resize(n);
		oldcorr.resize(n);
		for(size_t i=0; i<n; i++)
		{
			UG_ASSERT(px[0]->size() == px[i]->size(), "all vectors must have same size");
			CloneVector(corr[i], *px[0]);
			CloneVector(oldcorr[i], *px[0]);

			//PrintStorageType(*px[i]);
			//PrintStorageType(corr[i]);
			//PrintStorageType(oldcorr[i]);
		}

		std::vector<vector_type *> pTestVectors;

		std::vector<double> corrnorm(n, m_dPrecision*10);

		std::vector<std::string> testVectorDescription;

		m_spPrecond->init(m_pA);
		pcl::ProcessCommunicator pc;

		for(size_t iteration=0; iteration<m_maxIterations; iteration++)
		{
			UG_LOG("iteration " << iteration << "\n");
			// 0. normalize

			if(m_pB)
			{
				for(size_t i=0; i<n; i++)
				{
					double d = sqrt( EnergyProd(*px[i], B, *px[i]));
					//UG_LOG("sqrt(<x[" << i << "], Bx[" << i << "]>) = " << d << "\n");
					(*px[i]) *= 1/d;
				}
			}
			else
			{
				for(size_t i=0; i<n; i++)
				{
					double d = px[i]->two_norm();
					//UG_LOG("sqrt(<x[" << i << "], x[" << i << "]>) = " << d << "\n");
					(*px[i]) *= 1 / d;
				}
			}



			// 1. before calculating new correction, save old correction
			for(size_t i=0; i<n; i++)
				std::swap(oldcorr[i], corr[i]);

			//  2. compute rayleigh quotient, residuum, apply preconditioner, compute corrections norm
			size_t nrofconverged=0;


			for(size_t i=0; i<n; i++)
			{
				// 2a. compute rayleigh quotients
				// lambda = <x, Ax>/<x,x>
				// todo: replace with MatMult
//				UG_LOG("m_pA has storage type "); PrintStorageType(*m_pA); UG_LOG(", and vector px[" << i << "] has storage type"); PrintStorageType(*px[i]); UG_LOG("\n");
				// px can be set to unique because of two_norm

#ifdef UG_PARALLEL
				px[i]->change_storage_type(PST_CONSISTENT);
				defect.set_storage_type(PST_ADDITIVE);
#endif
				m_pA->apply(defect, *px[i]);


#ifdef UG_PARALLEL
				defect.change_storage_type(PST_UNIQUE);
				px[i]->change_storage_type(PST_UNIQUE);
#endif
				lambda[i] = px[i]->dotprod(defect); // / <px[i], px[i]> = 1.0.
				//UG_LOG("lambda[" << i << "] = " << lambda[i] << "\n");

				// 2b. calculate residuum
				// defect = A px[i] - lambda[i] B px[i]
				if(m_pB)
				{
					// todo: replace with MatMultAdd
					//MatMultAddDirect(defect, 1.0, defect, -lambda[i], *m_pB, *px[i]);
#ifdef UG_PARALLEL
					defect.change_storage_type(PST_ADDITIVE);
					px[i]->change_storage_type(PST_CONSISTENT);
#endif
					MatMultAddDirect(defect, 1.0, defect, -lambda[i], B, *px[i]);
				}
				else
					VecScaleAdd(defect, 1.0, defect, -lambda[i], *px[i]);

				// 2c. check if converged


#ifdef UG_PARALLEL
				defect.change_storage_type(PST_UNIQUE);
#endif
				corrnorm[i] = defect.two_norm();
#ifdef UG_PARALLEL
				defect.change_storage_type(PST_UNIQUE);
#endif
				//UG_LOG("corrnorm[" << i << "] = " << corrnorm[i] << "\n");
				if(corrnorm[i] < m_dPrecision)
					nrofconverged++;
				if(corrnorm[i] < 1e-12)
					continue;
				// 2d. apply preconditioner
				m_spPrecond->apply(corr[i], defect);
				corr[i] *= 1/ sqrt( EnergyProd(corr[i], B, corr[i]));
#ifdef UG_PARALLEL
				corr[i].change_storage_type(PST_UNIQUE);
#endif
			}

			// output
			UG_LOG("================================================\n");
			UG_LOG("iteration " << iteration << "\n");

			for(size_t i=0; i<n; i++)
				UG_LOG(i << " lambda: " << std::setw(14) << lambda[i] << " defect: " <<
						std::setw(14) << corrnorm[i] << "\n");
			UG_LOG("\n");

			if(nrofconverged==n)
			{
				UG_LOG("all eigenvectors converged\n");
				return true;
			}

			// 5. add Testvectors
			//UG_LOG("5. add Testvectors\n");

			pTestVectors.clear();
			testVectorDescription.clear();
			for(size_t i=0; i<n; i++)
			{
				if(m_iPINVIT == 1)
				{
					VecScaleAdd(corr[i], -1.0, corr[i], 1.0, *px[i]);
					testVectorDescription.push_back(std::string("ev - corr [") + ToString(i) + std::string("]") );
				}
				else
				{
					pTestVectors.push_back(px[i]);
					testVectorDescription.push_back(std::string("eigenvector [") + ToString(i) + std::string("]") );
					if(corrnorm[i] < m_dPrecision)
							continue;
					pTestVectors.push_back(&corr[i]);
					testVectorDescription.push_back(std::string("correction [") + ToString(i) + std::string("]") );

					if(m_iPINVIT >= 3)
					{
						pTestVectors.push_back(&oldcorr[i]);
						testVectorDescription.push_back(std::string("old correction [") + ToString(i) + std::string("]") );

						if(iteration == 0)
						{
							for(size_t j=0; j<px[i]->size(); j++)
								oldcorr[i][j] = (*px[i])[j] * urand(-1.0, 1.0);
						}
					}
				}
			}
			/*for(size_t i=0; i<testVectorDescription.size(); i++)
			{	UG_LOG(testVectorDescription[i] << "\n");	} */

			size_t iNrOfTestVectors = pTestVectors.size();

			// 5. compute reduced Matrices rA, rB
			//UG_LOG("5. compute reduced Matrices rA, rB\n");
			rA.resize(iNrOfTestVectors, iNrOfTestVectors);
			rB.resize(iNrOfTestVectors, iNrOfTestVectors);

			if(m_pB)
				MultiEnergyProd(*m_pB, &pTestVectors[0], rB, iNrOfTestVectors);
			else
				MultiScalProd(&pTestVectors[0], rB, iNrOfTestVectors);


			/////// Remove linear depended vectors //////////////
			std::vector<bool> bUse;
			bUse.resize(iNrOfTestVectors, true);

			{
				DenseMatrix<VariableArray2<double> > mat = rB;

				for(size_t i=0; i<mat.num_rows(); i++)
				{
					for(size_t j=0; j<i; j++)
					{
						if(!bUse[i]) continue;
						double val = mat(i, j)/mat(j,j);
						mat(i,j) = 0;
						for(size_t k=j+1; k<mat.num_rows(); k++)
							mat(i,k) -= val*mat(j, k);
					}
					if(mat(i,i) < 1e-8) bUse[i] = false;
					else bUse[i] = true;
				}

				//for(size_t i=0; i<iNrOfTestVectors; i++)
					//UG_LOG("Using " << testVectorDescription[i] << (bUse[i] ?"\n":" NOT\n"));
				//PrintMatrix(rB, "rB");
				//PrintMatrix(mat, "mat");
			}


			std::vector<vector_type *> pTestVectors2 = pTestVectors;
			std::vector<std::string> testVectorDescription2 = testVectorDescription;
			pTestVectors.clear();
			testVectorDescription.clear();
			for(size_t i=0; i<iNrOfTestVectors; i++)
			{
				if(bUse[i])
				{
					pTestVectors.push_back(pTestVectors2[i]);
					testVectorDescription.push_back(testVectorDescription2[i]);
				}
			}

			iNrOfTestVectors = pTestVectors.size();

			// 5. compute reduced Matrices rA, rB
			rA.resize(iNrOfTestVectors, iNrOfTestVectors);
			rB.resize(iNrOfTestVectors, iNrOfTestVectors);

			if(m_pB)
				MultiEnergyProd(*m_pB, &pTestVectors[0], rB, iNrOfTestVectors);
			else
				MultiScalProd(&pTestVectors[0], rB, iNrOfTestVectors);

			MultiEnergyProd(*m_pA, &pTestVectors[0], rA, iNrOfTestVectors);

			PrintMaple(rA, "rA");
			PrintMaple(rB, "rB");

			// output

			// 6. solve reduced eigenvalue problem

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
					UG_LOG(i << ".: " << r_lambda[i] << "\n");
			}

			UG_LOG("Lambda > 0: \n");
			for(size_t i=nrzero; i<r_lambda.size(); i++)
				UG_LOG(i << ".: " << r_lambda[i] << "\n");

			UG_LOG("\n");


			UG_LOG("evs: \n");
			for(size_t c=nrzero; c < std::min(nrzero+n, r_ev.num_cols()); c++)
			{
				UG_LOG("ev [" << c << "]:\n");
				for(size_t r=0; r<r_ev.num_rows(); r++)
					if(dabs(r_ev(r, c)) > 1e-9 )
						UG_LOG(std::setw(14) << r_ev(r, c) << "   " << testVectorDescription[r] << "\n");
				UG_LOG("\n\n");
			}


			for(size_t i=0; i<iNrOfTestVectors; i++)
				pTestVectors[i]->change_storage_type(PST_UNIQUE);
			for(size_t i=0; i<n; i++)
				px[i]->change_storage_type(PST_UNIQUE);
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

	void set_pinvit(size_t iPINVIT)
	{
		m_iPINVIT = iPINVIT;
		UG_ASSERT(iPINVIT >=1 && iPINVIT <= 3, "i has to be >= 1 and <= 3, but is " << iPINVIT);
	}

};

} // namespace ug


#endif // __H__UG__LIB_ALGEBRA__PINVIT_H__
