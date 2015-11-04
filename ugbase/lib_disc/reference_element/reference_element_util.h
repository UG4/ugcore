/*
 * reference_element_util.h
 *
 *  Created on: 03.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_UTIL__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_UTIL__


namespace ug{

/// returns the reference element dimension at run-time
inline int ReferenceElementDimension(ReferenceObjectID roid)
{
	switch(roid)
	{
		case ROID_VERTEX: return 0;
		case ROID_EDGE: return 1;
		case ROID_TRIANGLE: return 2;
		case ROID_QUADRILATERAL: return 2;
		case ROID_TETRAHEDRON: return 3;
		case ROID_OCTAHEDRON: return 3;
		case ROID_PYRAMID: return 3;
		case ROID_PRISM: return 3;
		case ROID_HEXAHEDRON: return 3;
		default:UG_THROW("ReferenceElementDimension: ReferenceObjectId "
						 <<roid<<" not found.");
	}
}

/// returns the Center of a reference element at run-time
template <int dim>
inline MathVector<dim> ReferenceElementCenter(ReferenceObjectID roid);

template <>
inline MathVector<1> ReferenceElementCenter<1>(ReferenceObjectID roid)
{
	switch(roid)
	{
		case ROID_EDGE: return MathVector<1>(0.5);
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 1.");
	}
}

template <>
inline MathVector<2> ReferenceElementCenter<2>(ReferenceObjectID roid)
{
	switch(roid)
	{
		case ROID_TRIANGLE: return MathVector<2>(1.0/3.0, 1.0/3.0);
		case ROID_QUADRILATERAL: return MathVector<2>(0.5, 0.5);
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 2.");
	}
}

template <>
inline MathVector<3> ReferenceElementCenter<3>(ReferenceObjectID roid)
{
	switch(roid)
	{
		case ROID_TETRAHEDRON: return MathVector<3>(0.25, 0.25, 0.25);
		case ROID_PYRAMID: return MathVector<3>(2./5., 2./5., 1./5.);
		case ROID_PRISM: return MathVector<3>(2./6., 2./6., 0.5);
		case ROID_HEXAHEDRON: return MathVector<3>(0.5, 0.5, 0.5);
		case ROID_OCTAHEDRON: return MathVector<3>(2./6., 2./6., 0.0);
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 3.");
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_UTIL__ */
