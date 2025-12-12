/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_BASE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_BASE__


namespace ug {

/// helper class to store dimension and id of a midpoint of a sub-element
struct MidID
{
		MidID() : dim(0), id(0) {};
		MidID(size_t dim_, size_t id_) : dim(dim_), id(id_) {};
		size_t dim;
		size_t id;
};

///	a singleton class that returns a new id for each type
class UniqueFVGeomIDProvider{
	public:
		static UniqueFVGeomIDProvider& inst(){
			static UniqueFVGeomIDProvider instance;
			return instance;
		}

		size_t new_id()	{return ++m_id;}

	private:
		UniqueFVGeomIDProvider() : m_id(0)	{}
		size_t m_id;
};

///	This method associates a unique unsigned integer value with each type.
template <typename TType>
size_t GetUniqueFVGeomID()
{
	static size_t typeID = UniqueFVGeomIDProvider::inst().new_id();
	//** Unkomment the next line to see the correspondence of the id's and the types
	//UG_LOG ("--> Assigning id " << typeID << " to type " << typeid(TType).name () << " <--\n");
	//**
	return typeID;
}


/// base class for all FVGeometries
class FVGeometryBase {};

} // end namespace ug

#endif