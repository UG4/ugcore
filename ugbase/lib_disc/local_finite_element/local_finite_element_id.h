/*
 * local_finite_element_id.h
 *
 *  Created on: 16.11.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LOCAL_FINITE_ELEMENT_ID__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LOCAL_FINITE_ELEMENT_ID__

#include <sstream>

namespace ug{


// Doxygen group
///////////////////////////////////////////////////////////////////////////////
/**
 * \brief provides Local Finite Elements.
 *
 * The Local Finite Element section is used to describe finite element spaces
 * by their definition on reference elements.
 *
 * A Finite Element is defined as a triplet \f$ \{ K, P, \Sigma \} \f$
 * (See e.g.
 * Ciarlet, P., "Basis Error Estimates for Elliptic Problems", North-Holland,
 * Amsterdam, 1991, p. 93;
 * or Ern, A. and Guermond J.L., "Theory and Practice of Finite Elements",
 *  Springer, 2004, p. 19), where
 *
 * <ol>
 * <li> \f$ K \f$ is a compact, connected, Lipschitz subset of \f$\mathbb{R}^d\f$
 * 		with non-empty interior
 * <li> \f$ P \f$ is a vector space of functions \f$p: K \mapsto \mathbb{R}^m \f$
 * 		with an integer \f$m > 0 \f$ (usually \f$m=1\f$ or \f$m=d\f$)
 * <li> \f$\Sigma\f$ is a set of \f$ n_{sh} \f$ linear forms \f$ \sigma_1, \dots,
 * 		\sigma_{n_{sh}} \f$ acting on the elements of \f$ P \f$, such that the
 * 		linear mapping
 * 		\f[
 * 			p \mapsto ( \sigma_1(p), \dots, \sigma_{n_{sh}}(p)) \in \mathbb{R}^{n_{sh}}
 * 		\f]
 * 		is bijective. These linear forms are called local degrees of freedom.
 * </ol>
 *
 * Since the mapping is bijective, there exist a basis \f$\{\phi_1, \dots,
 * \phi_{n_{sh}} \subset P\f$ such that
 * \f[
 * 	\sigma_i (\phi_j) = \delta_{ij}, \qquad 1 \leq i,j \leq n_{sh}.
 * \f]
 *
 * This set is called the set of local shape functions. The implemented
 * counterpart is the class LocalShapeFunctionSet.
 *
 * The set of local degrees of freedom finds its counterpart in the class
 * ILocalDoFSet.
 *
 * \defgroup lib_disc_local_finite_elements Local Finite Elements
 * \ingroup lib_discretization
 */
///////////////////////////////////////////////////////////////////////////////

/// \ingroup lib_disc_local_finite_elements
/// @{

/// Identifier for Local Finite Elements
/**
 * This Class is used to distinguish between different local finite elements (lfe).
 * Each lfe has a unique Space Type (e.g. Lagrange, DG) and the order of trial
 * functions. It the function space is p-adaptive, the enum ADAPTIVE is set as
 * order.
 */
class LFEID
{
	public:
	/// Space Type
		enum SpaceType
		{
			NONE = -1,
			LAGRANGE = 0,
			CROUZEIX_RAVIART,
			PIECEWISE_CONSTANT,
			DG,
			MINI,
			USER_DEFINED,
			NUM_SPACE_TYPES
		};

	///	special possibilities for order
		enum {ADAPTIV = -1, INVALID = -10};

	public:
	///	default constructor
		LFEID() : m_type(NONE), m_order(INVALID) {}

	///	constructor with values
		LFEID(SpaceType type, int order) : m_type(type), m_order(order) {}

	///	returns the order of the local finite element
		int order() const {return m_order;}

	///	returns the type of the local finite element
		SpaceType type() const {return m_type;}

	///	equality check
		bool operator==(const LFEID& v) const
		{
			return (m_type == v.m_type && m_order==v.m_order);
		}

	///	inequality check
		bool operator!=(const LFEID& v) const {return !((*this)==v);}

	///	operator <
	/**	returns comparison which set id is regarded lesser
	 * The Local Finite Elements are ordered by type first and
	 * then by increasing order.
	 */
		bool operator<(const LFEID& v) const
		{
			if(m_type != v.m_type) return m_type < v.m_type;
			else return m_order < v.m_order;
		}

	///	operator >
		bool operator>(const LFEID& v) const
		{
			if(m_type != v.m_type) return m_type > v.m_type;
			else return m_order > v.m_order;
		}

	///	operator <=
		bool operator<=(const LFEID& v) const
		{
			return (*this < v || *this == v);
		}

	///	operator >=
		bool operator>=(const LFEID& v) const
		{
			return (*this > v || *this == v);
		}

		friend std::ostream& operator<<(std::ostream& out,	const LFEID& v);

	private:
	///	Space type
		SpaceType m_type;

	///	Order
		int m_order;
};

/// writes the Identifier to the output stream
std::ostream& operator<<(std::ostream& out,	const LFEID& v);

///	returns the LFEID for a combination of Space and order
LFEID ConvertStringToLFEID(const char* type, int order);

///	returns the LFEID
LFEID ConvertStringToLFEID(const char* type);

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LOCAL_FINITE_ELEMENT_ID__ */
