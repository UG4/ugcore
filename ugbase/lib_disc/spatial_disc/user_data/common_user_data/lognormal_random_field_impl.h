/*
 * Lognormal_Random_Field_impl.h
 *
 *  Created on: 17.12.2014
 *      Author: Ivo Muha, Martin Rupp
 *
 *
 *  Method described in Comput Visual Sci (2006) 9: 1Ð10, DOI 10.1007/s00791-006-0012-2
 *  Simulation of lognormal random fields with varying resolution scale and local average for Darcy flow
 *  For the generation of the d-dimensional Gaussian random field we choose
 *  a simple spectral method which can sim- ulate anisotropicly correlated fields.
 *  The Gaussian random field f (x) is realized as a superposition of a large number
 *  of randomly chosen harmonic modes following the method introduced by Kraichnan
 *
 */
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_RANDOM_FIELD_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_RANDOM_FIELD_IMPL__

#include "common/util/typename.h"
#include "lognormal_random_field.h"
#include "common/math/misc/math_util.h" // urand

#ifndef M_PI
#define M_PI    3.14159265358979323846264338327950288   /* pi */
#endif

namespace ug{

template <typename TData, int dim, typename TRet>
TRet LognormalRandomField<TData,dim,TRet>::evaluate(TData& D, const MathVector<dim>& x, number time, int si) const
{
	double k = eval_K(x);
	for(size_t i = 0; i < dim; ++i)
	{
		for(size_t j = 0; j < dim; ++j)
		{
			D[i][j] = 0.0;
			if (i==j)
				D[i][j] = k;
		}
	}
	return;
}

template <typename TData, int dim, typename TRet>
double LognormalRandomField<TData,dim,TRet>::gasdev()
{
	// from Numerical Recipes
	int iset;
	static double gset;
	double fac, rsq, v1, v2, x1, x2;

	iset = 0;
	if (iset == 0)
	{
		do
		{
			x1 = urand(0.0, 1.0);
			x2 = urand(0.0, 1.0);
			v1 = 2.0 * x1 - 1.0;
			v2 = 2.0 * x2 - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}
	else
	{
		iset = 0;
		return gset;
	}

}

template<typename TData, int dim, typename TRet>
double LognormalRandomField<TData, dim, TRet>::undev()
{
	return urand(0.0, 1.0);
}

template <typename TData, int dim, typename TRet>
double LognormalRandomField<TData,dim,TRet>::eval_K(const MathVector<dim> &x) const
{

	double result = 0.0;
	for(int i = 0; i < m_N; i++)
		result += cos(VecDot(m_vRandomQvec[i], x) + m_vRandomAlpha[i]);


	double f = m_dMean_f + sqrt(2*m_dSigma_f*m_dSigma_f/m_N)*result;

	if(m_bNoExp)
		return f;
	else
		return exp(f);
}

template <typename TData, int dim, typename TRet>
void LognormalRandomField<TData,dim,TRet>::set_config(size_t N, double mean_f, double sigma_f, double sigma)
{
	m_N = N;
	for(int j=0; j<dim; j++)
		m_sigma[j] = sqrt(1.0/sigma);
	m_dMean_f = mean_f;
	m_dSigma_f = sigma_f;
	m_dSigma = sigma;

	m_vRandomQvec.clear();
	m_vRandomAlpha.clear();

	MathVector<dim> q;

	m_vRandomQvec.resize(N);
	for(int j=0; j<dim; j++)
		for(int i = 0; i < m_N; i++)
			m_vRandomQvec[i][j] = m_sigma[j]*gasdev();


	for(int i = 0; i < m_N; i++)
		m_vRandomAlpha.push_back(undev()*2*M_PI);

//	UG_LOG("corrx = " << m_sigma[0] << " corry = " << m_sigma[1] << "\n");
//	PRINT_VECTOR(m_vRandomQvec, "m_vRandomQvec");
//	PRINT_VECTOR(m_vRandomAlpha, "m_vRandomAlpha");

}

template <typename TData, int dim, typename TRet>
std::string LognormalRandomField<TData,dim,TRet>::config_string() const
{
	std::stringstream ss;

	//ss << "LognormalRandomField < TData =  " << TypeName<TData>() << ", dim = " << dim << ", TRet = " << TypeName<TRet>() << " >\n";
	ss << " LognormalRandomField<" << dim << "d> ( m_N = " << m_N << ", m_dMean_f = " << m_dMean_f <<
			", m_dSigma_f = " << m_dSigma_f << ", sigma = " << m_dSigma << ", m_bNoExp = " << TrueFalseString(m_bNoExp) << ")";
	return ss.str();

}


} // ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__LOGNORMAL_RANDOM_FIELD_IMPL__
