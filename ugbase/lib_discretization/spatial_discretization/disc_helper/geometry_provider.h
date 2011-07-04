/*
 * geometry_provider.h
 *
 *  Created on: 06.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__GEOMETRY_PROVIDER__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__GEOMETRY_PROVIDER__

namespace ug{

/// Singleton, holding a Geometry for Discretization
/**
 * This class is used to wrap discretization geometries into a singleton, such
 * that every building-block uses the same geometry if they need the same. By
 * this double computations is avoided.
 */
class GeomProvider {

	// 	private constructor
		GeomProvider();

	// 	disallow copy and assignment (intentionally left unimplemented)
		GeomProvider(const GeomProvider&);
		GeomProvider& operator=(const GeomProvider&);

	// 	private destructor
		~GeomProvider(){};

	// 	geometry provider, holding the instance
		template <typename TFVGeom>
		inline static TFVGeom& inst()
		{
			static TFVGeom myInst;
			return myInst;
		};

	public:
	///	returns access to the singleton
		template <typename TFVGeom>
		inline static TFVGeom& get() {	return inst<TFVGeom>();}
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__GEOMETRY_PROVIDER__ */
