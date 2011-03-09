/*
 * diffusion_length.h
 *
 *  Created on: 29.10.2010
 *      Author: josefdubsky, andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__DIFFUSION_LENGTH__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__DIFFUSION_LENGTH__

// other ug4 modules
#include "common/common.h"

namespace ug{

/**
 *
 * \param[out]		DiffLengthSqInv		Inverse of squared Diffusion Length for each integration point
 * \param[in]		geo					Finite Volume Geometry
 */
template <typename TFVGeometry>
bool NSDiffLengthFivePoint(number DiffLengthSqInv[], const TFVGeometry& geo)
{
    //	dimension of element
	static const size_t dim = TFVGeometry::dim;

    // Check the SCVF - number of IPs (order of elements)
    UG_ASSERT(geo.scvf(0).num_ip() == 1, "Only implemented for first order.");

	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename TFVGeometry::SCVF& scvf = geo.scvf(i);

	//	get associated SubControlVolumes
		const typename TFVGeometry::SCV& scvFrom = geo.scv(scvf.from());
		const typename TFVGeometry::SCV& scvTo = geo.scv(scvf.to());

	//	Norm of Normal to SCVF squared
		const number normNormalSq = VecTwoNormSq(scvf.normal());

	//	average squared size of associated SCV
		number areaSCVSq = 0.5 * (scvFrom.volume() + scvTo.volume());
		areaSCVSq *= areaSCVSq;

		if(dim == 2)
		{
			DiffLengthSqInv[i] = 2 * normNormalSq/ areaSCVSq + 8 / normNormalSq;
		}
		else if (dim == 3)
		{
		//	Distance between edge midpoint and center of element, that are part
		//	of the SCVF
			const number distSq = VecDistanceSq(scvf.global_corner(0),
												scvf.global_corner(2));

			DiffLengthSqInv[i] = 2 * normNormalSq/ areaSCVSq + 8 * distSq / normNormalSq;
		}
		else
		{
			UG_LOG("ERROR in 'NSDiffLenghtAnsatz1':"
					"Not implemented for dimension " << dim << ".\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TFVGeometry>
bool NSDiffLengthRaw(number DiffLengthSqInv[], const TFVGeometry& geo)
{
    //	dimension of element
	static const size_t dim = TFVGeometry::dim;

    // Check the SCVF - number of IPs (order of elements)
    UG_ASSERT(geo.scvf(0).num_ip() == 1, "Only implemented for first order.");


	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename TFVGeometry::SCVF& scvf = geo.scvf(i);

	//	get associated SubControlVolumes
		const typename TFVGeometry::SCV& scvFrom = geo.scv(scvf.from());
		const typename TFVGeometry::SCV& scvTo = geo.scv(scvf.to());

	//	Norm of Normal to SCVF squared
		const number normNormalSq = VecTwoNormSq(scvf.normal());

	//	average squared size of associated SCV
		number areaSCVSq = 0.5 * (scvFrom.volume() + scvTo.volume());
		areaSCVSq *= areaSCVSq;

		if(dim == 2)
		{
			DiffLengthSqInv[i] = 1 / (0.5*areaSCVSq/normNormalSq + 3*normNormalSq/8);
		}
		else if (dim == 3)
		{
		//	Distance between edge midpoint and center of element, that are part
		//	of the SCVF
			const number distSq = VecDistanceSq(scvf.global_corner(0),
												scvf.global_corner(2));

			DiffLengthSqInv[i] = 1 / (0.5*areaSCVSq/normNormalSq + 8 * 3*distSq/8);
		}
		else
		{
			UG_LOG("ERROR in 'NSDiffLenghtAnsatz1':"
					"Not implemented for dimension " << dim << ".\n");
			return false;
		}
    }

//	we're done
	return true;
}

template <typename TFVGeometry>
bool NSDiffLengthCor(number DiffLengthSqInv[], const TFVGeometry& geo)
{
    //	dimension of element
	//static const size_t dim = TFVGeometry::dim;

    // Check the SCVF - number of IPs (order of elements)
    UG_ASSERT(geo.scvf(0).num_ip() == 1, "Only implemented for first order.");


	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
        // TODO: Must be implemented
        DiffLengthSqInv[i]=1.0;
    }

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__DIFFUSION_LENGTH__ */
