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

#ifndef __H__UG_boost_serialization
#define __H__UG_boost_serialization

//#include <boost/serialization/access.hpp>
//#include <boost/serialization/export.hpp>
//#include <boost/serialization/level.hpp>
#include <boost/serialization/nvp.hpp>
//#include <boost/serialization/version.hpp>

/* Include this file if you want to add serialization functionality to your types.
 * If you want to perform serialization, please make sure to also include
 * 'boost_serialization_routines.h', which defines serialization routines for
 * common types.*/


/* If a derived class is serialized and the base shall not be serialized,
 * this macro should be used during serialization of the derived class to declare
 * the base class nevertheless.*/
#define UG_EMPTY_BASE_CLASS_SERIALIZATION( clsDerived, clsBase)\
		boost::serialization::void_cast_register<clsDerived, clsBase>(\
				static_cast<clsDerived *>(NULL), static_cast<clsBase *>(NULL));

namespace ug{

/**	make_nvp is used to create a name-value pair during serialization. Given an
 * archive 'ar', one serializes data using 'ar & make_nvp("somename", someVar);'
 */
using boost::serialization::make_nvp;


///	different types of archives
enum ArchiveType {
	AT_DATA,
	AT_GUI
};

///	Provides custom information for different archives
template <class TArchive> // (ø)[[unused template argument]] // todo remove template parameter or fix if unintended
struct ArchiveInfo {
	static const ArchiveType TYPE = AT_DATA;
};


}//	end of namespace

#endif	//__H__UG_boost_serialization
