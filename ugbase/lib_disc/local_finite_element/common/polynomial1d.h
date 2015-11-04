
#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__COMMON__POLYNOMIAL1D__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__COMMON__POLYNOMIAL1D__

#include "common/math/ugmath.h"
#include <vector>

namespace ug{

/// \addtogroup lib_discretization
///	@{

/** base class for one dimensional polynomials
 * This class is used to represent polynomials in one variable. For the
 * evaluation the horner scheme is used. Note that using this representation
 * the computation of higher order derivatives turns out easier than by
 * hard coded implementations.
 */
class Polynomial1D
{
	public:
	///	Constructor producing zero polynomial of degree 'degree'
		Polynomial1D(size_t degree = 0)
			: m_vCoeff(degree+1, 0.0)
		{}

	///	Constructor passing coefficients for the polynomial
		Polynomial1D(const std::vector<number>& a)
			: m_vCoeff(a)
		{
		//	check that at least constant of polynomial set
			if(m_vCoeff.empty())
				m_vCoeff.resize(1, 0.0);
		};

	/**	returns the degree of the polynomial.
	 * This function returns the degree of the polynomial, i.e. the
	 * highest coefficient stored. Note that no checking is performed if
	 * the leading coefficient is zero.
	 */
		size_t degree() const {return m_vCoeff.size() - 1;}

	///	evaluate the value of the polynom at x
		number value(const number x) const
		{
		//	get degree of polynomial (is >= 0 by construction)
			const size_t deg = m_vCoeff.size() - 1;

		//	loop horner scheme
			number val = m_vCoeff[deg];
			for(size_t i = deg; i > 0; --i)
				val = m_vCoeff[i-1] + val * x;

		//	we're done
			return val;
		}

	///	returns the derivative of this polynomial as a polynomial
		Polynomial1D derivative() const
		{
		//	if only constant present, return empty Polynomial
			if(degree() == 0)
				return Polynomial1D();

		//	create empty polynomial of with correct size
			Polynomial1D tmpPol(degree() - 1);

		//	differentiate
			for(size_t i = 0; i <= tmpPol.degree(); ++i)
				tmpPol.m_vCoeff[i] = (i+1) * m_vCoeff[i+1];

		//	return derivative by copy
			return tmpPol;
		}

	///	multiply by a polynomial
		Polynomial1D& operator *=(const Polynomial1D& v)
		{
		//	new size of polynomial
			size_t newDeg = degree() + v.degree();

		//	create new coefficients
			std::vector<number> vNewCoeff(newDeg+1, 0.0);

		//	multiply
			for(size_t i = 0; i <= degree(); ++i)
				for(size_t j = 0; j <= v.degree(); ++j)
					vNewCoeff[i+j] += m_vCoeff[i] * v.m_vCoeff[j];

		//	Copy new coeffs
			m_vCoeff = vNewCoeff;

		//	we're done
			return *this;
		}

	///	multiply by a scalar
		Polynomial1D& operator *=(number scale)
		{
		//	multiply
			for(size_t i = 0; i <= degree(); ++i)
				m_vCoeff[i] *= scale;

		//	we're done
			return *this;
		}

	//	output
		friend std::ostream& operator<< (std::ostream& outStream, Polynomial1D& v);

	protected:
		void set_coefficients(const std::vector<number>& a)
		{
		//	assign coefficients
			m_vCoeff = a;

		//	check that at least constant of polynomial set
			if(m_vCoeff.empty())
				m_vCoeff.resize(1, 0.0);
		};

	private:
	//	vector holding the coefficients of the polynom
	//	An empty vector is the Polynomial p = 0;
	//	else we have p(x) = sum_i m_vCoeff[i] *x^i
		std::vector<number> m_vCoeff;
};

inline std::ostream& operator<< (std::ostream& outStream, Polynomial1D& v)
{
	for(size_t i = 0; i <= v.degree(); ++i)
	{
		outStream << v.m_vCoeff[i] << " *x^" << i;
		if(i != v.degree()) outStream << " + ";
	}
	return outStream;
}

/// @}
} // end namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__COMMON__POLYNOMIAL1D__ */
