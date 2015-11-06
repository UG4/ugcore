/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__COMMON__LAGRANGE1D__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__COMMON__LAGRANGE1D__

#include "./polynomial1d.h"
#include "common/math/ugmath.h"
#include <vector>

namespace ug{

/** Lagrange Polynomial for arbitrary points
 *
 */
class Lagrange1D
	: public Polynomial1D
{
	public:
	/**	constructor for lagrange polynomial i using interpolation points pos
	 * This constructor creates a lagrange polynomial with interpolation points
	 * given in pos for the i-th point, i.e. value(pos_i) == 1, value(pos_j) == 0
	 * for j != i. Therefore, it must hold that 0 <= i < pos.size()
	 */
		Lagrange1D(const size_t i, const std::vector<number>& vPos)
		{
		//	compute coefficients
			compute_coeffs(i, vPos);

		//	remember positions
			m_vPos = vPos;
		}

	///	returns the position of the i'th interpolation point
		number position(const size_t i) const
		{
			UG_ASSERT(i < m_vPos.size(), "Invalid index");
			return m_vPos[i];
		}

	protected:
	/// computes the coefficients for passed interpolation points
		void compute_coeffs(const size_t i, const std::vector<number>& vPos)
		{
		//	start coefficients
			std::vector<number> vStart(1, 1.0);

		//	start polynomial
			this->set_coefficients(vStart);

		//	help polynomial used for each factor
			std::vector<number> vFactor(2,1.0);

		//	scaling of polynom
			number scale = 1.0;

		//	fill coefficients
			for(size_t j = 0; j < vPos.size(); ++j)
			{
				if(j == i) continue;

			//	set first coefficent to minus position
				vFactor[0] = -vPos[j];
			//	create polynom for (-vPos, 1)
				Polynomial1D tmpPol(vFactor);
			//	multiply
				(*this) *= tmpPol;
			// 	multiply scale
				scale *= 1./(vPos[i]-vPos[j]);
			}

		//	multiply by scale
			(*this) *= scale;
		}

	private:
		std::vector<number> m_vPos;
};

/** EquiDistant Lagrange Function
 *
 */
class EquidistantLagrange1D
	: public Polynomial1D
{
	public:
	/** creates a lagrange polynomial with equidistant interpolation points
	 * \param[in] 	i		number of interpolation point, where polynom is 1
	 * \param[in]	degree	degree of polynom
	 */
		EquidistantLagrange1D(const size_t i, const size_t degree)
		{
			UG_ASSERT(i <= degree, "Only #degree shape functions.");
			compute_coeffs(i, degree);
		}

	///	returns the position of the i'th interpolation point
		static number position(const size_t i, const size_t degree)
		{
			UG_ASSERT(i <= degree, "Invalid index");
			return (number)i/(number)degree;
		}

	protected:
	/// computes the coefficients for passed interpolation points
		void compute_coeffs(const int i, const int p)
		{
		//	start coefficients
			std::vector<number> vStart(1, 1.0);

		//	start polynomial
			this->set_coefficients(vStart);

		//	help polynomial used for each factor
			std::vector<number> vFactor(2, p);

		//	scaling of polynom
			number scale = 1.0;

		//	fill coefficients
			for(int j = 0; j <= p; ++j)
			{
				if(j == i) continue;

			//	set first coefficent to minus position
				vFactor[0] = -j;
			//	create polynom for (-j, p)
				Polynomial1D tmpPol(vFactor);
			//	multiply
				(*this) *= tmpPol;
			// 	multiply scale
				scale *= 1./(i-j);
			}

		//	multiply by scale
			(*this) *= scale;
		}
};

/** Truncated EquiDistant Lagrange Function
 *
 * Creates for given order <tt>p</tt> and iterpolation point <tt>i</tt> the
 * polynomial \f[ \prod_{j=0}^{i-1} \frac{x - \frac{j}{p}}{\frac{i}{p} -
 * \frac{j}{p}} \f]
 */
class TruncatedEquidistantLagrange1D
	: public Polynomial1D
{
	public:
	/** creates a lagrange polynomial with equidistant interpolation points
	 * \param[in] 	i		number of interpolation point, where polynom is 1
	 * \param[in]	degree	degree of polynom
	 */
		TruncatedEquidistantLagrange1D(const size_t i, const size_t degree)
		{
			UG_ASSERT(i <= degree, "Only #degree shape functions.");

			compute_coeffs(i, degree);
		}

	///	returns the position of the i'th interpolation point
		static number position(const size_t i, const size_t degree)
		{
			UG_ASSERT(i <= degree, "Invalid index");
			return (number)i/(number)degree;
		}

	protected:
	/// computes the coefficients for passed interpolation points
		void compute_coeffs(const int i, const int p)
		{
		//	start coefficients
			std::vector<number> vStart(1, 1.0);

		//	start polynomial
			this->set_coefficients(vStart);

		//	help polynomial used for each factor
			std::vector<number> vFactor(2, p);

		//	scaling of polynom
			number scale = 1.0;

		//	fill coefficients
			for(int j = 0; j < i; ++j)
			{
			//	set first coefficent to minus position
				vFactor[0] = -j;
			//	create polynom for (-j, p)
				Polynomial1D tmpPol(vFactor);
			//	multiply
				(*this) *= tmpPol;
			// 	multiply scale
				scale *= 1./(i-j);
			}

		//	multiply by scale
			(*this) *= scale;
		}
};

/** Bounded EquiDistant Lagrange Function
 *
 * Creates for given order <tt>p</tt>, interpolation point <tt>i</tt> and upper
 * bound <tt>0 <= b <= p</tt> the polynomial
 * \f[ \prod_{\substack{j=0\\j\neq i}}^{b} \frac{x - \frac{j}{p}}{\frac{i}{p} -
 * \frac{j}{p}} \f]
 * Thus, it is a polynomial of order b.
 */
class BoundedEquidistantLagrange1D
	: public Polynomial1D
{
	public:
	/** creates a lagrange polynomial with equidistant interpolation points
	 * \param[in] 	i		number of interpolation point, where polynom is 1
	 * \param[in]	degree	degree of polynom
	 * \param[in]	bound	Point until lagrange points are included
	 */
		BoundedEquidistantLagrange1D(const size_t i, const size_t degree,
		                             const size_t bound)
		{
			UG_ASSERT(i <= bound, "Only #bound shape functions.");
			UG_ASSERT(bound <= degree, "Only #bound shape functions.");

		//	init coefficients
			compute_coeffs(i, degree, bound);
		}

	///	returns the position of the i'th interpolation point
		static number position(const size_t i, const size_t degree)
		{
			UG_ASSERT(i <= degree, "Invalid index");
			return (number)i/(number)degree;
		}

	protected:
	/// computes the coefficients for passed interpolation points
		void compute_coeffs(const int i, const int p, const int b)
		{
		//	start coefficients
			std::vector<number> vStart(1, 1.0);

		//	start polynomial
			this->set_coefficients(vStart);

		//	help polynomial used for each factor
			std::vector<number> vFactor(2, p);

		//	scaling of polynom
			number scale = 1.0;

		//	fill coefficients
			for(int j = 0; j <= b; ++j)
			{
				if(j == i) continue;

			//	set first coefficent to minus position
				vFactor[0] = -j;
			//	create polynom for (-j, p)
				Polynomial1D tmpPol(vFactor);
			//	multiply
				(*this) *= tmpPol;
			// 	multiply scale
				scale *= 1./(i-j);
			}

		//	multiply by scale
			(*this) *= scale;
		}
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__COMMON__LAGRANGE1D__ */
