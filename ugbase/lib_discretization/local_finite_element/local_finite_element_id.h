/*
 * local_finite_element_id.h
 *
 *  Created on: 16.11.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_FINITE_ELEMENT_ID__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_FINITE_ELEMENT_ID__

#include <sstream>

namespace ug{

/// \ingroup lib_discretization_finite_elements
/// @{

/// Identifier for local finite elements
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
			DG,
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

		friend std::ostream& operator<<(std::ostream& out,	const LFEID& v);

	private:
	///	Space type
		SpaceType m_type;

	///	Order
		int m_order;
};

/// writes the Identifier to the output stream
inline std::ostream& operator<<(std::ostream& out,	const LFEID& v)
{
	std::stringstream ss;
	if(v.m_order >= 0) ss << v.m_order;
	else if(v.m_order == LFEID::ADAPTIV) ss << "adaptive";
	else ss << "invalid";

	switch(v.m_type)
	{
		case LFEID::LAGRANGE: out << "(Lagrange, " << ss.str() << ")"; break;
		case LFEID::DG: out << "(DG, " << ss.str() << ")"; break;
		case LFEID::USER_DEFINED: out << "(User defined, " << ss.str() << ")"; break;
		default: out << "(unknown, " << ss.str() << ")";
	}
	return out;
}

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_FINITE_ELEMENT_ID__ */
