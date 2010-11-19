/*
 * local_shape_function_set_id.h
 *
 *  Created on: 16.11.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_SHAPE_FUNCTION_SET_ID__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_SHAPE_FUNCTION_SET_ID__

namespace ug{

/// \ingroup lib_discretization_local_shape_function_set
/// @{

/** Identifier for local shape function sets
 * This Class is used to distingush between different
 * local shape function sets. Each set has a unique
 * Space Type (e.g. Lagrange, DG) and the order of trial
 * functions.
 */
class LocalShapeFunctionSetID
{
	public:
	/// Space Type
		enum SpaceType
		{
			LAGRANGE,
			DG,
		};

	public:
	///	default constructor
		LocalShapeFunctionSetID()
		{}

	///	constructor with values
		LocalShapeFunctionSetID(SpaceType type, size_t order)
			: m_type(type), m_order(order)
		{}

	///	equality check
		bool operator==(const LocalShapeFunctionSetID& v) const
		{
			return (m_type == v.m_type && m_order==v.m_order);
		}

	///	inequality check
		bool operator!=(const LocalShapeFunctionSetID& v) const
		{
			return !((*this)==v);
		}

	///	operator <
	/**	returns comparison which set id is reguarded lesser
	 * The Local Shape Function Sets are ordered by type first and
	 * then by increasing order.
	 */
		bool operator<(const LocalShapeFunctionSetID& v) const
		{
			if(m_type != v.m_type)
				return m_type < v.m_type;
			else
				return m_order < v.m_order;
		}

	/// writes the Identifier to the output stream
		friend std::ostream& operator<<(std::ostream& out,
		                                const LocalShapeFunctionSetID& v)
		{
			switch(v.m_type)
			{
				case LAGRANGE: out << "(Lagrange, " << v.m_order << ")"; break;
				case DG: out << "(DG, " << v.m_order << ")"; break;
			}
			return out;
		}

	private:
	//	Space type
		SpaceType m_type;

	//	Order
		size_t m_order;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LOCAL_SHAPE_FUNCTION_SET_ID__ */
