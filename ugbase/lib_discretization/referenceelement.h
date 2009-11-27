/*
 * referenceelement.hpp
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCEELEMENT__
#define __H__LIBDISCRETIZATION__REFERENCEELEMENT__

#include <cassert>
#include <iostream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "../lib_grid/lib_grid.h"

namespace ug{


class ReferenceElement {

	public:

		/*** GENERAL INFORMATIONS ***/
		/* dimension in which reference element lives */
		virtual int dimension() const = 0;

		/* size of reference element */
		virtual double size() const = 0;

		/* coordinates of corner of reference element (i=0..numberOfCorners) */
		virtual const vector2& coordsOfReferenceCorner(const int i) const = 0;

		/*** NUMBER OF GEOMETRIC OBJECTS OF REFERENCE ELEMENT ***/
		/* number of GeomObjects of dimension 'dim' of reference element */
		virtual unsigned int numberOfGeomObjectsOfRefElem(const int dim) const = 0;

		/* number of corners of reference element */
		virtual unsigned int numberOfCornersOfRefElem() const = 0;

		/* number of edges of reference element */
		virtual unsigned int numberOfEdgesOfRefElem() const = 0;

		/* number of faces of reference element */
		virtual unsigned int numberOfFacesOfRefElem() const = 0;

		/* number of volumes of reference element */
		virtual unsigned int numberOfVolumesOfRefElem() const = 0;

		/*** NUMBER OF GEOMETRIC OBJECTS OF SUB GEOMETRIC ELEMENTS ***/
		/* number of Geometric Objects of dimension dim_i of GeomObject with dim_j, nr j (dim_i < dim_j) */
		virtual unsigned int numberOfGeomObjectsOfGeomObject(const int dim_i, const int dim_j, const int j) const = 0;

		virtual unsigned int numberOfCornersOfEdge(const int j) const = 0;
		virtual unsigned int numberOfCornersOfFace(const int j) const = 0;
		virtual unsigned int numberOfCornersOfVolume(const int j) const = 0;

		virtual unsigned int numberOfEdgesOfFace(const int j) const = 0;
		virtual unsigned int numberOfEdgesOfVolume(const int j) const = 0;

		virtual unsigned int numberOfFacesOfVoume(const int j) const = 0;

		/*** ID's OF GEOMETRIC OBJECTS OF REFERENCE ELEMENT ***/
		/* Id of GeometricObjects of dimension i, nr. i, of Geometric object of dimension j, nr j */
		virtual int IdOfGeomObjectOfGeomObject(const int dim_i, const int i, const int dim_j, const int j) const = 0;

		/* Id of Corner i of Edge j */
		virtual int CornerIdOfEdge(const int i, const int j) const = 0;

		/* Id of Corner i of Face j */
		virtual int CornerIdOfFace(const int i, const int j) const = 0;

		/* Id of Corner i of Volume j */
		virtual int CornerIdOfVolume(const int i, const int j) const = 0;

		/* Id of Edge i of Face j */
		virtual int EdgeIdOfFace(const int i, const int j) const = 0;

		/* Id of Edge i of Volume j */
		virtual int EdgeIdOfVolume(const int i, const int j) const = 0;

		/* Id of Face i of Volume j */
		virtual int FaceIdOfVolume(const int i, const int j) const = 0;

		/*** MAPPINGS ***/
		/* Mapping Local to Global
		 *
		 * GlobalCorners = Coordinates of Global Corners
		 * Local = Local Coordinate
		 * Global = Global Coordinate
		 *
		 * return: true if ok, false if error occured
		 * */
		//virtual bool mapLocalToGlobal(const vector3 GlobalCorners[], const vector2 Local, vector3& Global) const = 0;
		//virtual bool mapLocalToGlobal(const vector3 GlobalCorners[], const vector2 Local[], vector3 Global[], const int n) const = 0;
		//virtual bool Trafo(const vector3 GlobalCorners[], const vector2 Local, MathMatrix<2,2>& InvTrafo, number& det) const = 0;

		/* virtual destructor */
		virtual ~ReferenceElement()
		{}

		bool printInfo()
		{
			std::string GeomObjects[4] ={"Corner", "Edge", "Face", "Volume"};

			std::cout << "Reference Element Info: " << std::endl;
			std::cout << "----------------------- " << std::endl;

			std::cout << "Size: " << this->size() << std::endl;
			std::cout << "Dimension where Reference Element lives: " << this->dimension() << std::endl;
			std::cout << "Reference Corners: ";
			for(unsigned int i = 0; i< this->numberOfCornersOfRefElem(); i++)
			{
				std::cout << "(" << this->coordsOfReferenceCorner(i)[0] << ","<< this->coordsOfReferenceCorner(i)[1] <<") ";
			}
			std::cout << std::endl;

			for(int i = this->dimension(); i>=0 ;i--)
			{
				std::cout << "Number of " << GeomObjects[i] << "s: " << this->numberOfGeomObjectsOfRefElem(i) << std::endl;
			}

			for(int dim_i = this->dimension(); dim_i>=0 ;dim_i--)
			{
				for(unsigned int i=0; i < this->numberOfGeomObjectsOfRefElem(dim_i); i++)
				{
					std::cout << GeomObjects[dim_i] << " with id '" << i << "' contains the following GeomObjects:" << std::endl;
					for(int dim_j=dim_i; dim_j>= 0; dim_j--)
					{
						std::cout << this->numberOfGeomObjectsOfGeomObject(dim_i,i,dim_j) << " " << GeomObjects[dim_j] << "s with id: ";
						for(unsigned int j=0; j< this->numberOfGeomObjectsOfGeomObject(dim_i,i,dim_j); j++)
						{
							std::cout << this->IdOfGeomObjectOfGeomObject(dim_i,i,dim_j,j) << " ";
						}
						std::cout << std::endl;
					}
				}
			}

			return true;
		}

	private:


};

template<class TElem>
class ReferenceElementFor : public ReferenceElement{

private:
	ReferenceElementFor()
	{};
};

template <>
class ReferenceElementFor<Triangle> : public ReferenceElement{
	private:
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME = 3, MAXDIM};
		enum{MYDIM = FACE};
		enum{MAXOBJECTS = 3};


	public:
		static const size_t Dim = FACE;


		ReferenceElementFor()
		{
			initializeArrays();
		}

		/* Dimension where reference element lives */
		int dimension() const
		{
			return m_dimension;
		}

		/* size of reference triangle */
		double size() const
		{
			return m_size;
		}

		/* coordinates of reference corner Nr i (i=0..numberOfCorners) */
		const vector2& coordsOfReferenceCorner(const int i) const
		{
			return m_coordsOfReferenceCorner[i];
		}

		/* number of corners of reference triangle */
		unsigned int numberOfCornersOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[POINT];
		}

		unsigned int numberOfEdgesOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[EDGE];
		}

		unsigned int numberOfFacesOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[FACE];
		}

		unsigned int numberOfVolumesOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[VOLUME];
		}

		unsigned int numberOfGeomObjectsOfRefElem(const int dim) const
		{
			return m_numberOfGeomObjectsOfRefElem[dim];
		}

		unsigned int numberOfGeomObjectsOfGeomObject(const int dim_i, const int i, const int dim_j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[dim_i][i][dim_j];
		}

		unsigned int numberOfCornersOfEdge(const int j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[EDGE][j][POINT];
		}
		unsigned int numberOfCornersOfFace(const int j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[FACE][j][POINT];
		}
		unsigned int numberOfCornersOfVolume(const int j) const
		{
			return 0;
		}
		unsigned int numberOfEdgesOfFace(const int j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[FACE][j][EDGE];
		}
		unsigned int numberOfEdgesOfVolume(const int j) const
		{
			return 0;
		}
		unsigned int numberOfFacesOfVoume(const int j) const
		{
			return 0;
		}

		int IdOfGeomObjectOfGeomObject(const int dim_i, const int i, const int dim_j, const int j) const
		{
			assert((dim_i <= MYDIM) && (dim_i>= dim_j) && "ERROR in ReferenceElement::IdOfGeomObjectOfGeomObject(const int dim_i, const int i, const int dim_j, const int j)"
					": You must have: 'dim_i < Dimension Reference Element' and 'dim_j <= dim_i'");
			return m_id[dim_i][i][dim_j][j];
		}

		int CornerIdOfEdge(const int i, const int j) const
		{
			return m_id[EDGE][i][POINT][j];
		}

		int CornerIdOfFace(const int i, const int j) const
		{
			return m_id[FACE][i][POINT][j];
		}

		int CornerIdOfVolume(const int i, const int j) const
		{
			return -1;
		}

		int EdgeIdOfFace(const int i, const int j) const
		{
			return m_id[FACE][i][EDGE][j];
		}

		int EdgeIdOfVolume(const int i, const int j) const
		{
			return -1;
		}

		int FaceIdOfVolume(const int i, const int j) const
		{
			return -1;
		}


		bool mapLocalToGlobal(const MathVector<3> GlobalCorners[], const vector2 Local, MathVector<3>& Global) const
		{
			MathVector<3> a10, a20;

			VecSubtract(a10, GlobalCorners[1], GlobalCorners[0]);
			VecSubtract(a20, GlobalCorners[2], GlobalCorners[0]);

			VecScale(a10, a10, Local[0]);
			VecScale(a20, a20, Local[1]);

			Global = GlobalCorners[0];
			VecAdd(Global, Global, a10);
			VecAdd(Global, Global, a20);

			return true;
		}

		bool mapLocalToGlobal(const MathVector<3> GlobalCorners[], const vector2 Local[], MathVector<3> Global[], const int n) const
		{
			MathVector<3> a10, a20, sa10, sa20;

			VecSubtract(a10, GlobalCorners[1], GlobalCorners[0]);
			VecSubtract(a20, GlobalCorners[2], GlobalCorners[0]);

			for(int i=0; i<n; ++i)
			{
				VecScale(sa10, a10, Local[i][0]);
				VecScale(sa20, a20, Local[i][1]);

				Global[i] = GlobalCorners[0];
				VecAdd(Global[i], Global[i], sa10);
				VecAdd(Global[i], Global[i], sa20);
			}
			return true;
		}

		bool Trafo(const MathVector<2> GlobalCorners[], const MathVector<2> Local, MathMatrix<2,2>& InvTrafo, number& det) const
		{

			MathMatrix<2,2> Trafo;

			Trafo[0][0] = GlobalCorners[1][0] - GlobalCorners[0][0];
			Trafo[0][1] = GlobalCorners[1][1] - GlobalCorners[0][1];

			Trafo[1][0] = GlobalCorners[2][0] - GlobalCorners[0][0];
			Trafo[1][1] = GlobalCorners[2][1] - GlobalCorners[0][1];

			det = Trafo[0][0]*Trafo[1][1] - Trafo[0][1]*Trafo[1][0];

			assert(det != 0.0);

			InvTrafo[0][0] = Trafo[1][1] / det;
			InvTrafo[1][0] = -Trafo[1][0] / det;
			InvTrafo[0][1] = -Trafo[0][1] / det;
			InvTrafo[1][1] = Trafo[0][0] / det;

			return true;
		}

		// Achtung: Hier sollte GlobalCorners eigentlich ein MathVector<2> sein.
		bool Trafo(const MathVector<3> GlobalCorners[], const MathVector<2> Local, MathMatrix<2,2>& InvTrafo, number& det) const
		{

			MathMatrix<3,3> Trafo;

			Trafo[0][0] = GlobalCorners[1][0] - GlobalCorners[0][0];
			Trafo[0][1] = GlobalCorners[1][1] - GlobalCorners[0][1];

			Trafo[1][0] = GlobalCorners[2][0] - GlobalCorners[0][0];
			Trafo[1][1] = GlobalCorners[2][1] - GlobalCorners[0][1];

			det = Trafo[0][0]*Trafo[1][1] - Trafo[0][1]*Trafo[1][0];

			assert(det != 0.0);

			InvTrafo[0][0] = Trafo[1][1] / det;
			InvTrafo[1][0] = -Trafo[1][0] / det;
			InvTrafo[0][1] = -Trafo[0][1] / det;
			InvTrafo[1][1] = Trafo[0][0] / det;

			return true;
		}

		~ReferenceElementFor()
		{}

	private:
		/* Dimension, in which the reference triangle lives */
		unsigned int m_dimension;
		/* size of reference triangle */
		double m_size;
		/* number of Geometric Objects of Reference Element
		 * (m_numberOfGeomObjectsOfRefElem[dim] = number of GeomObjects of dimension dim) */
		unsigned int m_numberOfGeomObjectsOfRefElem[MAXDIM];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int m_numberOfGeomObjectsOfGeomObject[MYDIM+1][MAXOBJECTS][MYDIM+1];
		/* coordinates of Reference Corner */
		vector2 m_coordsOfReferenceCorner[3];
		// indices of GeomObjects
		int m_id[MYDIM+1][MAXOBJECTS][MYDIM+1][MAXOBJECTS];

		void initializeArrays()
		{
			// dimension, where reference triangle lives
			m_dimension = MYDIM;

			// size of reference element
			m_size = 1./2.;

			//number of Geometric Objects
			m_numberOfGeomObjectsOfRefElem[POINT] = 3;
			m_numberOfGeomObjectsOfRefElem[EDGE] = 3;
			m_numberOfGeomObjectsOfRefElem[FACE] = 1;
			m_numberOfGeomObjectsOfRefElem[VOLUME] = 0;

			// number of Geometric Objects
			m_numberOfGeomObjectsOfGeomObject[FACE][0][POINT] = 3;
			m_numberOfGeomObjectsOfGeomObject[FACE][0][EDGE] = 3;
			m_numberOfGeomObjectsOfGeomObject[FACE][0][FACE] = 1;

			m_numberOfGeomObjectsOfGeomObject[EDGE][0][EDGE] = 1;
			m_numberOfGeomObjectsOfGeomObject[EDGE][1][EDGE] = 1;
			m_numberOfGeomObjectsOfGeomObject[EDGE][2][EDGE] = 1;

			m_numberOfGeomObjectsOfGeomObject[EDGE][0][POINT] = 2;
			m_numberOfGeomObjectsOfGeomObject[EDGE][1][POINT] = 2;
			m_numberOfGeomObjectsOfGeomObject[EDGE][2][POINT] = 2;

			m_numberOfGeomObjectsOfGeomObject[POINT][0][POINT] = 1;
			m_numberOfGeomObjectsOfGeomObject[POINT][1][POINT] = 1;
			m_numberOfGeomObjectsOfGeomObject[POINT][2][POINT] = 1;

			//reset m_id to -1
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(unsigned int k=0; k<=MYDIM; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
							m_id[i][j][k][l] = -1;
						}

			//self references:
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<m_numberOfGeomObjectsOfRefElem[i]; ++j)
				{
					m_id[i][j][i][0] = j;
				}

			//Edges of Face
			for(unsigned int i=0; i<m_numberOfGeomObjectsOfRefElem[EDGE]; ++i)
			{
				m_id[FACE][0][EDGE][i] = i;
			}

			// Points of Face
			for(unsigned int i=0; i<m_numberOfGeomObjectsOfRefElem[POINT]; ++i)
			{
				m_id[FACE][0][POINT][i] = i;
			}


			// Points of Edges
			// edge 0 = (0,1)
			m_id[EDGE][0][POINT][0] = 0;
			m_id[EDGE][0][POINT][1] = 1;
			// edge 1 = (1,2)
			m_id[EDGE][1][POINT][0] = 1;
			m_id[EDGE][1][POINT][1] = 2;
			// edge 2 = (2,0)
			m_id[EDGE][2][POINT][0] = 2;
			m_id[EDGE][2][POINT][1] = 0;

			// Reference Corners
			m_coordsOfReferenceCorner[0].x = 0.0;
			m_coordsOfReferenceCorner[0].y = 0.0;
			m_coordsOfReferenceCorner[1].x = 1.0;
			m_coordsOfReferenceCorner[1].y = 0.0;
			m_coordsOfReferenceCorner[2].x = 0.0;
			m_coordsOfReferenceCorner[2].y = 1.0;
		}


};


template <>
class ReferenceElementFor<Quadrilateral> : public ReferenceElement{
	private:
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME = 3, MAXDIM};
		enum{MYDIM = FACE};
		enum{MAXOBJECTS = 4};


	public:
		static const size_t Dim = FACE;


		ReferenceElementFor()
		{
			initializeArrays();
		}

		/* Dimension where reference element lives */
		int dimension() const
		{
			return m_dimension;
		}

		/* size of reference triangle */
		double size() const
		{
			return m_size;
		}

		/* coordinates of reference corner Nr i (i=0..numberOfCorners) */
		const vector2& coordsOfReferenceCorner(const int i) const
		{
			return m_coordsOfReferenceCorner[i];
		}

		/* number of corners of reference triangle */
		unsigned int numberOfCornersOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[POINT];
		}

		unsigned int numberOfEdgesOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[EDGE];
		}

		unsigned int numberOfFacesOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[FACE];
		}

		unsigned int numberOfVolumesOfRefElem() const
		{
			return m_numberOfGeomObjectsOfRefElem[VOLUME];
		}

		unsigned int numberOfGeomObjectsOfRefElem(const int dim) const
		{
			return m_numberOfGeomObjectsOfRefElem[dim];
		}

		unsigned int numberOfGeomObjectsOfGeomObject(const int dim_i, const int i, const int dim_j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[dim_i][i][dim_j];
		}

		unsigned int numberOfCornersOfEdge(const int j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[EDGE][j][POINT];
		}
		unsigned int numberOfCornersOfFace(const int j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[FACE][j][POINT];
		}
		unsigned int numberOfCornersOfVolume(const int j) const
		{
			return 0;
		}
		unsigned int numberOfEdgesOfFace(const int j) const
		{
			return m_numberOfGeomObjectsOfGeomObject[FACE][j][EDGE];
		}
		unsigned int numberOfEdgesOfVolume(const int j) const
		{
			return 0;
		}
		unsigned int numberOfFacesOfVoume(const int j) const
		{
			return 0;
		}

		int IdOfGeomObjectOfGeomObject(const int dim_i, const int i, const int dim_j, const int j) const
		{
			assert((dim_i <= MYDIM) && (dim_i>= dim_j) && "ERROR in ReferenceElement::IdOfGeomObjectOfGeomObject(const int dim_i, const int i, const int dim_j, const int j)"
					": You must have: 'dim_i < Dimension Reference Element' and 'dim_j <= dim_i'");
			return m_id[dim_i][i][dim_j][j];
		}

		int CornerIdOfEdge(const int i, const int j) const
		{
			return m_id[EDGE][i][POINT][j];
		}

		int CornerIdOfFace(const int i, const int j) const
		{
			return m_id[FACE][i][POINT][j];
		}

		int CornerIdOfVolume(const int i, const int j) const
		{
			return -1;
		}

		int EdgeIdOfFace(const int i, const int j) const
		{
			return m_id[FACE][i][EDGE][j];
		}

		int EdgeIdOfVolume(const int i, const int j) const
		{
			return -1;
		}

		int FaceIdOfVolume(const int i, const int j) const
		{
			return -1;
		}


		bool mapLocalToGlobal(const MathVector<3> GlobalCorners[], const vector2 Local, MathVector<3>& Global) const
		{
			VecScaleAdd(Global, (1.-Local[0])*(1.-Local[1]), GlobalCorners[0],
								Local[0]*(1.-Local[1])     , GlobalCorners[1],
								Local[0]*Local[1]          , GlobalCorners[2],
								(1.-Local[0])*Local[1]     , GlobalCorners[3]);

			return true;
		}

		bool mapLocalToGlobal(const MathVector<3> GlobalCorners[], const vector2 Local[], MathVector<3> Global[], const int n) const
		{
			for(int i=0; i<n; ++i)
			{
				VecScaleAdd(Global[i], (1.-Local[i][0])*(1.-Local[i][1]), GlobalCorners[0],
									   Local[i][0]*(1.-Local[i][1])     , GlobalCorners[1],
									   Local[i][0]*Local[i][1]          , GlobalCorners[2],
									   (1.-Local[i][0])*Local[i][1]     , GlobalCorners[3]);
			}
			return true;
		}

		bool Trafo(const MathVector<2> GlobalCorners[], const MathVector<2> Local, MathMatrix<2,2>& InvTrafo, number& det) const
		{

			return false;
		}

		// Achtung: Hier sollte GlobalCorners eigentlich ein MathVector<2> sein.
		bool Trafo(const MathVector<3> GlobalCorners[], const MathVector<2> Local, MathMatrix<2,2>& InvTrafo, number& det) const
		{
			number a;
			MathMatrix<3,3> Trafo;


			a = 1. - Local[1];
			Trafo[0][0] = a*(GlobalCorners[1][0] - GlobalCorners[0][0]) + Local[1]*(GlobalCorners[2][0] - GlobalCorners[3][0]);
			Trafo[0][1] = a*(GlobalCorners[1][1] - GlobalCorners[0][1]) + Local[1]*(GlobalCorners[2][1] - GlobalCorners[3][1]);

			a = 1. - Local[0];
			Trafo[1][0] = a*(GlobalCorners[3][0] - GlobalCorners[0][0]) + Local[0]*(GlobalCorners[2][0] - GlobalCorners[1][0]);
			Trafo[1][1] = a*(GlobalCorners[3][1] - GlobalCorners[0][1]) + Local[0]*(GlobalCorners[2][1] - GlobalCorners[1][1]);

			det = Trafo[0][0]*Trafo[1][1] - Trafo[0][1]*Trafo[1][0];

			assert(det != 0.0);

			InvTrafo[0][0] = Trafo[1][1] / det;
			InvTrafo[1][0] = -Trafo[1][0] / det;
			InvTrafo[0][1] = -Trafo[0][1] / det;
			InvTrafo[1][1] = Trafo[0][0] / det;

			return true;
		}

		~ReferenceElementFor()
		{}

	private:
		/* Dimension, in which the reference triangle lives */
		unsigned int m_dimension;
		/* size of reference triangle */
		double m_size;
		/* number of Geometric Objects of Reference Element
		 * (m_numberOfGeomObjectsOfRefElem[dim] = number of GeomObjects of dimension dim) */
		unsigned int m_numberOfGeomObjectsOfRefElem[MAXDIM];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int m_numberOfGeomObjectsOfGeomObject[MYDIM+1][MAXOBJECTS][MYDIM+1];
		/* coordinates of Reference Corner */
		vector2 m_coordsOfReferenceCorner[4];
		// indices of GeomObjects
		int m_id[MYDIM+1][MAXOBJECTS][MYDIM+1][MAXOBJECTS];

		void initializeArrays()
		{
			// dimension, where reference triangle lives
			m_dimension = MYDIM;

			// size of reference element
			m_size = 1.;

			//number of Geometric Objects
			m_numberOfGeomObjectsOfRefElem[POINT] = 4;
			m_numberOfGeomObjectsOfRefElem[EDGE] = 4;
			m_numberOfGeomObjectsOfRefElem[FACE] = 1;
			m_numberOfGeomObjectsOfRefElem[VOLUME] = 0;

			// number of Geometric Objects
			m_numberOfGeomObjectsOfGeomObject[FACE][0][POINT] = 4;
			m_numberOfGeomObjectsOfGeomObject[FACE][0][EDGE] = 4;
			m_numberOfGeomObjectsOfGeomObject[FACE][0][FACE] = 1;

			m_numberOfGeomObjectsOfGeomObject[EDGE][0][EDGE] = 1;
			m_numberOfGeomObjectsOfGeomObject[EDGE][1][EDGE] = 1;
			m_numberOfGeomObjectsOfGeomObject[EDGE][2][EDGE] = 1;
			m_numberOfGeomObjectsOfGeomObject[EDGE][3][EDGE] = 1;

			m_numberOfGeomObjectsOfGeomObject[EDGE][0][POINT] = 2;
			m_numberOfGeomObjectsOfGeomObject[EDGE][1][POINT] = 2;
			m_numberOfGeomObjectsOfGeomObject[EDGE][2][POINT] = 2;
			m_numberOfGeomObjectsOfGeomObject[EDGE][3][POINT] = 2;

			m_numberOfGeomObjectsOfGeomObject[POINT][0][POINT] = 1;
			m_numberOfGeomObjectsOfGeomObject[POINT][1][POINT] = 1;
			m_numberOfGeomObjectsOfGeomObject[POINT][2][POINT] = 1;
			m_numberOfGeomObjectsOfGeomObject[POINT][3][POINT] = 1;

			//reset m_id to -1
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(unsigned int k=0; k<=MYDIM; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
							m_id[i][j][k][l] = -1;
						}

			//self references:
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<m_numberOfGeomObjectsOfRefElem[i]; ++j)
				{
					m_id[i][j][i][0] = j;
				}

			//Edges of Face
			for(unsigned int i=0; i<m_numberOfGeomObjectsOfRefElem[EDGE]; ++i)
			{
				m_id[FACE][0][EDGE][i] = i;
			}

			// Points of Face
			for(unsigned int i=0; i<m_numberOfGeomObjectsOfRefElem[POINT]; ++i)
			{
				m_id[FACE][0][POINT][i] = i;
			}


			// Points of Edges
			// edge 0 = (0,1)
			m_id[EDGE][0][POINT][0] = 0;
			m_id[EDGE][0][POINT][1] = 1;
			// edge 1 = (1,2)
			m_id[EDGE][1][POINT][0] = 1;
			m_id[EDGE][1][POINT][1] = 2;
			// edge 2 = (2,3)
			m_id[EDGE][2][POINT][0] = 2;
			m_id[EDGE][2][POINT][1] = 3;
			// edge 3 = (3,0)
			m_id[EDGE][3][POINT][0] = 3;
			m_id[EDGE][3][POINT][1] = 0;

			// Reference Corners
			m_coordsOfReferenceCorner[0].x = 0.0;
			m_coordsOfReferenceCorner[0].y = 0.0;
			m_coordsOfReferenceCorner[1].x = 1.0;
			m_coordsOfReferenceCorner[1].y = 0.0;
			m_coordsOfReferenceCorner[2].x = 1.0;
			m_coordsOfReferenceCorner[2].y = 1.0;
			m_coordsOfReferenceCorner[3].x = 0.0;
			m_coordsOfReferenceCorner[3].y = 1.0;
		}

};


template <class TElem>
class reference_element_traits
{};

template <>
class reference_element_traits<Triangle>
{
	public:
		typedef MathVector<2> LocalMathVector;

		enum
		{
			RefDim = 2,
			NumberCorners = 3
		};
};

template <>
class reference_element_traits<Quadrilateral>
{
	public:
		typedef MathVector<2> LocalMathVector;

		enum
		{
			RefDim = 2,
			NumberCorners = 4
		};
};


}

#endif /* __H__LIBDISCRETIZATION__REFERENCEELEMENT__ */
