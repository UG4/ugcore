/*
 * consistent_gravity.h
 *
 *  Created on: 26.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONSISTENT_GRAVITY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONSISTENT_GRAVITY__

// other ug4 modules
#include "common/common.h"
#include <vector>

namespace ug{

template <int dim>
inline bool PrepareConsistentGravity(	MathVector<dim>* vConsGravity,
										int* standard_gravity,
										const int coe,
										const MathVector<dim>* vCorners,
										const number* vDensity,
										const MathVector<dim>& PhysicalGravity)
{
	UG_ASSERT(0, "Not implemented.");
	return false;
}

template <>
inline bool PrepareConsistentGravity<2>(	MathVector<2>* vConsGravity,
											int* standard_gravity,
											const int coe,
											const MathVector<2>* vCorners,
											const number* vDensity,
											const MathVector<2>& PhysicalGravity)
{


//	static constants
	static const size_t _X_ = 0;
	static const size_t _Y_ = 1;

//	Reference Mappings
	static ReferenceMapping<ReferenceQuadrilateral, 2> QuadMapping;
	static ReferenceMapping<ReferenceTriangle, 2> TriangleMapping;

//	Number of corners
//	const int coe = vCorners.size();

	MathVector<2> LocalPoint,ShiftedGlobalPoint;
	MathVector<2> LocalGravityAt000, LocalGravityAt110;
	MathMatrix<2,2> JT;

	double DensityIP, Diff;
	double gradient;

	if(PhysicalGravity.x()==0.0)
		*standard_gravity = 1;
	else
		*standard_gravity = 0;

	if(coe==3 && *standard_gravity==1)
	{
		//new consistent gravity for TRIANGLE (cf. Frolkovic98)

		TriangleMapping.update(&vCorners[0]);
		TriangleMapping.jacobian_transposed_inverse(JT,vCorners[0]);

		for(int i = 1; i < coe; ++i)
		{
			VecSubtract(ShiftedGlobalPoint,vCorners[i],vCorners[0]);
			Diff = ShiftedGlobalPoint.y();

			ShiftedGlobalPoint.y() = 0.0;
			TransposedMatVecMult(LocalPoint,JT,ShiftedGlobalPoint);

			DensityIP = vDensity[0] * (1 - LocalPoint.x() - LocalPoint.y()); //GN(0)=1-local(0)-local(1)

			for(int j = 1; j < coe; ++j)
			{
				DensityIP += vDensity[j] * LocalPoint[j-1];//GN(j)=local(j-1) for j=1,2
			}

			gradient = (JT[0][0]*(vDensity[1]-vDensity[0]) + JT[0][1]*(vDensity[2]-vDensity[0]));
			vConsGravity[i][0] = Diff * PhysicalGravity.y() * gradient;
			vConsGravity[i][1] = Diff * PhysicalGravity.y() * (vDensity[i] + DensityIP)*0.5;
		}
	}
	else
	{
		//generalizied consistent gravity of Voss & Souza for local elements
		switch(coe)
		{
		case 4:
			/* QUADRILATERAL */
			QuadMapping.update(&vCorners[0]);

			/* compute the local gravity at local corner (0,0) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 0.0;
			QuadMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt000, JT, PhysicalGravity);

			/* compute the local gravity at local corner (1,1) */
			LocalPoint.x() = 1.0; LocalPoint.y() = 1.0;
			QuadMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt110, JT, PhysicalGravity);

			vConsGravity[0][_X_] = 0.0; vConsGravity[3][_X_] = 0.0;
			vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
			vConsGravity[2][_X_] = LocalGravityAt110[_X_]*(vDensity[2] + vDensity[3])*0.5;

			vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
			vConsGravity[2][_Y_] = LocalGravityAt110[_Y_]*(vDensity[1] + vDensity[2])*0.5;
			vConsGravity[3][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[3])*0.5;

			break;

		case 3:
			/* TRIANGLE */
			TriangleMapping.update(&vCorners[0]);

			/* compute the local gravity at local corner (0,0) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 0.0;
			TriangleMapping.jacobian_transposed(JT, LocalPoint);  //Transposed hier falsch?
			MatVecMult(LocalGravityAt000, JT, PhysicalGravity);

			vConsGravity[0][_X_] = 0.0; vConsGravity[2][_X_] = 0.0;
			vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;

			vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
			vConsGravity[2][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[2])*0.5;

			break;

		default:
			UG_LOG("In 'PrepareConsistentGravity': Element type not found.\n");
			return false;
		}
	}

 	return true;
}

template <>
inline bool PrepareConsistentGravity<3>(	MathVector<3>* vConsGravity,
											int* standard_gravity,
											const int coe,
											const MathVector<3>* vCorners,
											const number* vDensity,
											const MathVector<3>& PhysicalGravity)
{
//	static constants
	static const size_t _X_ = 0;
	static const size_t _Y_ = 1;
	static const size_t _Z_ = 2;

//	Reference Mappings
	static ReferenceMapping<ReferenceHexahedron, 3> HexMapping;
	static ReferenceMapping<ReferenceTetrahedron, 3> TetMapping;
	static ReferenceMapping<ReferencePrism, 3> PrismMapping;
	static ReferenceMapping<ReferencePyramid, 3> PyramidMapping;

//	Number of corners
//	const int coe = vCorners.size();

	MathVector<3> LocalPoint,ShiftedGlobalPoint;
	MathVector<3> 	LocalGravityAt000, LocalGravityAt110,
					LocalGravityAt101, LocalGravityAt011;
	MathMatrix<3,3> JT;

	double DensityIP, Diff;
	double gradient;

	if(coe==3 && *standard_gravity==1)
	{
		//new consistent gravity for TETRAHEDRON (cf. Frolkovic98)
		//TODO new consistent gravity for TETRAHEDRON needs to be tested (compared to UG3 results)

		TetMapping.update(&vCorners[0]);
		TetMapping.jacobian_transposed_inverse(JT,vCorners[0]);

		for(int i = 1; i < coe; ++i)
		{
			VecSubtract(ShiftedGlobalPoint,vCorners[i],vCorners[0]);
			Diff = ShiftedGlobalPoint.z();

			ShiftedGlobalPoint.z() = 0.0;
			TransposedMatVecMult(LocalPoint,JT,ShiftedGlobalPoint);

			DensityIP = vDensity[0] * (1 - LocalPoint.x() - LocalPoint.y() - LocalPoint.z()); //GN(0)=1-local(0)-local(1)-local(2)

			for(int j = 1; j < coe; ++j)
			{
				DensityIP += vDensity[j] * LocalPoint[j-1];//GN(j)=local(j-1) for j=1,2,3
			}
			for(int j = 0; j < 2; ++j)
			{
				gradient = JT[j][0]*(vDensity[1]-vDensity[0]) + JT[j][1]*(vDensity[2]-vDensity[0]) + JT[j][2]*(vDensity[3]-vDensity[0]);
				vConsGravity[i][j] = Diff * PhysicalGravity.z() * gradient;
			}
			vConsGravity[i][2] = Diff * PhysicalGravity.z() * (vDensity[i] + DensityIP)*0.5;
		}
	}
	else
	{
		switch(coe)
		{
		case 8:
			/*  HEXAHEDRON  */
			/****************/
			HexMapping.update(&vCorners[0]);

			/* compute the local gravity at local corner (0,0,0) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
			HexMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt000, JT, PhysicalGravity);

			/* compute the local gravity at local corner (1,1,0) */
			LocalPoint.x() = 1.0; LocalPoint.y() = 1.0; LocalPoint.z() = 0.0;
			HexMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt110, JT, PhysicalGravity);

			/* compute the local gravity at local corner (1,0,1) */
			LocalPoint.x() = 1.0; LocalPoint.y() = 0.0; LocalPoint.z() = 1.0;
			HexMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt101, JT, PhysicalGravity);

			/* compute the local gravity at local corner (0,1,1) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 1.0; LocalPoint.z() = 1.0;
			HexMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt011, JT, PhysicalGravity);

			vConsGravity[0][_X_] = 0.0; vConsGravity[3][_X_] = 0.0;
			vConsGravity[4][_X_] = 0.0; vConsGravity[7][_X_] = 0.0;
			vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
			vConsGravity[2][_X_] = LocalGravityAt110[_X_]*(vDensity[2] + vDensity[3])*0.5;
			vConsGravity[5][_X_] = LocalGravityAt101[_X_]*(vDensity[4] + vDensity[5])*0.5;
			vConsGravity[6][_X_] = LocalGravityAt011[_X_]*(vDensity[6] + vDensity[7])*0.5;

			vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
			vConsGravity[4][_Y_] = 0.0; vConsGravity[5][_Y_] = 0.0;
			vConsGravity[2][_Y_] = LocalGravityAt110[_Y_]*(vDensity[1] + vDensity[2])*0.5;
			vConsGravity[3][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[3])*0.5;
			vConsGravity[6][_Y_] = LocalGravityAt101[_Y_]*(vDensity[5] + vDensity[6])*0.5;
			vConsGravity[7][_Y_] = LocalGravityAt011[_Y_]*(vDensity[4] + vDensity[7])*0.5;

			vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0;
			vConsGravity[2][_Z_] = 0.0; vConsGravity[3][_Z_] = 0.0;
			vConsGravity[4][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[4])*0.5;
			vConsGravity[5][_Z_] = LocalGravityAt101[_Z_]*(vDensity[1] + vDensity[5])*0.5;
			vConsGravity[6][_Z_] = LocalGravityAt110[_Z_]*(vDensity[2] + vDensity[6])*0.5;
			vConsGravity[7][_Z_] = LocalGravityAt011[_Z_]*(vDensity[3] + vDensity[7])*0.5;

			break;

        case 6:
			/*  PRISM  */
			/***********/
    		PrismMapping.update(&vCorners[0]);

			/* compute the local gravity at local corner (0,0,0) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
			PrismMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt000, JT, PhysicalGravity);

			/* compute the local gravity at local corner (1,0,1) */
			LocalPoint.x() = 1.0; LocalPoint.y() = 0.0; LocalPoint.z() = 1.0;
			PrismMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt101, JT, PhysicalGravity);

			/* compute the local gravity at local corner (0,1,1) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 1.0; LocalPoint.z() = 1.0;
			PrismMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt011, JT, PhysicalGravity);

			vConsGravity[0][_X_] = 0.0; vConsGravity[2][_X_] = 0.0;
			vConsGravity[3][_X_] = 0.0; vConsGravity[5][_X_] = 0.0;
			vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
			vConsGravity[4][_X_] = LocalGravityAt011[_X_]*(vDensity[3] + vDensity[4])*0.5;

			vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0;
			vConsGravity[3][_Y_] = 0.0; vConsGravity[4][_Y_] = 0.0;
			vConsGravity[2][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[2])*0.5;
			vConsGravity[5][_Y_] = LocalGravityAt011[_Y_]*(vDensity[3] + vDensity[5])*0.5;

			vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0; vConsGravity[2][_Z_] = 0.0;
			vConsGravity[3][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[3])*0.5;
			vConsGravity[4][_Z_] = LocalGravityAt101[_Z_]*(vDensity[1] + vDensity[4])*0.5;
			vConsGravity[5][_Z_] = LocalGravityAt011[_Z_]*(vDensity[2] + vDensity[5])*0.5;
			break;

        case 5:
			/*  PYRAMID  */
			/*************/
        	PyramidMapping.update(&vCorners[0]);

 			/* compute the local gravity at local corner (0,0,0) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
			TetMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt000, JT, PhysicalGravity);

			/* compute the local gravity at local corner (1,1,0) */
			LocalPoint.x() = 1.0; LocalPoint.y() = 1.0; LocalPoint.z() = 0.0;
			TetMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt110, JT, PhysicalGravity);

			vConsGravity[0][_X_] = 0.0; vConsGravity[3][_X_] = 0.0; vConsGravity[4][_X_] = 0.0;
			vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;
			vConsGravity[2][_X_] = LocalGravityAt110[_X_]*(vDensity[2] + vDensity[3])*0.5;

			vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0; vConsGravity[4][_Y_] = 0.0;
			vConsGravity[2][_Y_] = LocalGravityAt110[_Y_]*(vDensity[1] + vDensity[2])*0.5;
			vConsGravity[3][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[3])*0.5;

			vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0;
			vConsGravity[2][_Z_] = 0.0; vConsGravity[3][_Z_] = 0.0;
			vConsGravity[4][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[4])*0.5;
			break;

        case 4:
			/* TETRAHEDRON  */
			/****************/
    		TetMapping.update(&vCorners[0]);

			/* compute the local gravity at local corner (0,0,0) */
			LocalPoint.x() = 0.0; LocalPoint.y() = 0.0; LocalPoint.z() = 0.0;
			TetMapping.jacobian_transposed(JT, LocalPoint);
			MatVecMult(LocalGravityAt000, JT, PhysicalGravity);

			vConsGravity[0][_X_] = 0.0; vConsGravity[2][_X_] = 0.0; vConsGravity[3][_X_] = 0.0;
			vConsGravity[1][_X_] = LocalGravityAt000[_X_]*(vDensity[0] + vDensity[1])*0.5;

			vConsGravity[0][_Y_] = 0.0; vConsGravity[1][_Y_] = 0.0; vConsGravity[3][_Y_] = 0.0;
			vConsGravity[2][_Y_] = LocalGravityAt000[_Y_]*(vDensity[0] + vDensity[2])*0.5;

			vConsGravity[0][_Z_] = 0.0; vConsGravity[1][_Z_] = 0.0; vConsGravity[2][_Z_] = 0.0;
			vConsGravity[3][_Z_] = LocalGravityAt000[_Z_]*(vDensity[0] + vDensity[3])*0.5;
			break;

        default:
        	UG_LOG("In 'PrepareConsistentGravity': Element type not found.\n");
			return false;
		}
	}

	return true;
}

template <int dim>
bool ComputeConsistentGravity(	MathVector<dim>& ConsistentGravity,
								const size_t coe,
								const MathVector<dim>& LocalCoord,
								const MathMatrix<dim, dim>& JTInv,
								const MathVector<dim>* vLocalGrad,
								const MathVector<dim>* vConsGravity,
								const int standard_gravity)
{
	//TODO: Überprüfe für Pyramiden!!!!!! siehe UG3


//	Clear ConsistentGravity
	MathVector<dim> LocalGravity;
	VecSet(LocalGravity, 0.0);

	int dm1;
	dm1 = dim - 1;

	if(coe==dim+1 && standard_gravity==1)
	{
		//new consistent gravity for TRIANGLE and TETRAHEDRON (cf. Frolkovic98)

		for(size_t d = 0; d < dim; ++d)
		{
			LocalGravity[d] = vConsGravity[d+1][dm1];
		}

	//	Multiply by JacobianTransposedInverse
		MatVecMult(ConsistentGravity, JTInv, LocalGravity);

		ConsistentGravity[0] -= vConsGravity[1][0]*LocalCoord[0];  //GN(j)=local(j-1) for j=1,2
		ConsistentGravity[0] -= vConsGravity[2][0]*LocalCoord[1];
		if(dim==3)
		{
			ConsistentGravity[0] -= vConsGravity[3][0]*LocalCoord[2];	//GN(j)=local(j-1) for j=1,2,3!
			ConsistentGravity[1] -= vConsGravity[1][1]*LocalCoord[0];
			ConsistentGravity[1] -= vConsGravity[2][1]*LocalCoord[1];
			ConsistentGravity[1] -= vConsGravity[3][1]*LocalCoord[2];
		}
	}
	else
	{
		//generalized consistent gravity of Voss & Souza for local elements

		//	Loop shape functions
		for(size_t sh = 0; sh < coe; ++sh)
		{
			for(size_t d = 0; d < dim; ++d)
			{
				LocalGravity[d] += vConsGravity[sh][d] * vLocalGrad[sh][d];
			}
		}

		//	Multiply by JacobianTransposedInverse
		MatVecMult(ConsistentGravity, JTInv, LocalGravity);
	}

		return true;
	}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONSISTENT_GRAVITY__ */
