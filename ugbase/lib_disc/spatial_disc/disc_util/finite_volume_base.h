/*
 * finite_volume_base.h
 *
 *  Created on: 11.01.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_BASE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_BASE__


namespace ug{

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
template <class TType>
size_t GetUniqueFVGeomID()
{
	static size_t typeID = UniqueFVGeomIDProvider::inst().new_id();
	return typeID;
}


/// base class for all FVGeometries
class FVGeometryBase {};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_BASE__ */
