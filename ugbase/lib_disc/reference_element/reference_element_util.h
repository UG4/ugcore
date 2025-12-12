/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_UTIL__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_UTIL__


namespace ug {

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

#endif