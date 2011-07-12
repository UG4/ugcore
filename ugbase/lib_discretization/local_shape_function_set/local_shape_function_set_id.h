/*
 * local_shape_function_set_id.h
 *
 *  Created on: 16.11.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_SHAPE_FUNCTION_SET_ID__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_SHAPE_FUNCTION_SET_ID__

#include <sstream>

namespace ug{

/// \ingroup lib_discretization_local_shape_function_set
/// @{

/// Identifier for local shape function sets
/**
 * This Class is used to distingush between different local shape function sets.
 * Each set has a unique Space Type (e.g. Lagrange, DG) and the order of trial
 * functions. It the function space is p-adaptive, the enum ADAPTIVE is set as
 * order.
 */
class LSFSID
{
	public:
	/// Space Type
		enum SpaceType
		{
			NONE = -1,
			LAGRANGE = 0,
			DG,
			NUM_SPACE_TYPES
		};

	///	special possibilities for order
		enum {ADAPTIV = -1, INVALID = -10};

	public:
	///	default constructor
		LSFSID() : m_type(NONE), m_order(INVALID) {}

	///	constructor with values
		LSFSID(SpaceType type, int order) : m_type(type), m_order(order) {}

	///	returns the order of the shape function set
		int order() const {return m_order;}

	///	returns the type of the shape function set
		SpaceType type() const {return m_type;}

	///	equality check
		bool operator==(const LSFSID& v) const
		{
			return (m_type == v.m_type && m_order==v.m_order);
		}

	///	inequality check
		bool operator!=(const LSFSID& v) const {return !((*this)==v);}

	///	operator <
	/**	returns comparison which set id is reguarded lesser
	 * The Local Shape Function Sets are ordered by type first and
	 * then by increasing order.
	 */
		bool operator<(const LSFSID& v) const
		{
			if(m_type != v.m_type) return m_type < v.m_type;
			else return m_order < v.m_order;
		}

		friend std::ostream& operator<<(std::ostream& out,	const LSFSID& v);

	private:
	///	Space type
		SpaceType m_type;

	///	Order
		int m_order;
};

/// writes the Identifier to the output stream
inline std::ostream& operator<<(std::ostream& out,	const LSFSID& v)
{
	std::stringstream ss;
	if(v.m_order >= 0) ss << v.m_order;
	else if(v.m_order == LSFSID::ADAPTIV) ss << "adaptive";
	else ss << "invalid";

	switch(v.m_type)
	{
		case LSFSID::LAGRANGE: out << "(Lagrange, " << ss.str() << ")"; break;
		case LSFSID::DG: out << "(DG, " << ss.str() << ")"; break;
		default: out << "(unknown, " << ss.str() << ")";
	}
	return out;
}

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_SHAPE_FUNCTION_SET_ID__ */
