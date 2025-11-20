/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

/* Include this file if you want to perform serialization, since it defines
 * basic serialization routines for common types. You won't need this include
 * if you're just writing serialization routines for your own types. Include
 * 'boost_serialization.h' in this case.*/

#ifndef __H__UG_boost_serialization_routines
#define __H__UG_boost_serialization_routines

#include "boost_serialization.h"
#include "math/ugmath_types.h"
 
BOOST_CLASS_IMPLEMENTATION(ug::vector1, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(ug::vector2, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(ug::vector3, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(ug::vector4, boost::serialization::object_serializable);

namespace boost{
namespace serialization{

	template <typename Archive>
	void serialize(Archive& ar, ug::vector1& v, const unsigned int version)
	{
		ar & ug::make_nvp("x", v[0]);
	}

	template <typename Archive>
	void serialize(Archive& ar, ug::vector2& v, const unsigned int version)
	{
		ar & ug::make_nvp("x", v[0]);
		ar & ug::make_nvp("y", v[1]);
	}

	template <typename Archive>
	void serialize(Archive& ar, ug::vector3& v, const unsigned int version)
	{
		ar & ug::make_nvp("x", v[0]);
		ar & ug::make_nvp("y", v[1]);
		ar & ug::make_nvp("z", v[2]);
	}

	template <typename Archive>
	void serialize(Archive& ar, ug::vector4& v, const unsigned int version)
	{
		ar & ug::make_nvp("x", v[0]);
		ar & ug::make_nvp("y", v[1]);
		ar & ug::make_nvp("z", v[2]);
		ar & ug::make_nvp("w", v[3]);
	}
}

}//	end of namespace

#endif
