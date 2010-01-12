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


class BaseReferenceElement {

	public:

		/*** GENERAL INFORMATIONS ***/
		/* dimension in which reference element lives */
		virtual int dimension() const = 0;

		/* size of reference element */
		virtual double size() const = 0;

		/* coordinates of corner of reference element (i=0..numberOfCorners) */
		virtual const vector2& corner(int i) const = 0;

		/*** NUMBER OF GEOMETRIC OBJECTS OF REFERENCE ELEMENT ***/
		/* number of GeomObjects of dimension 'dim' of reference element */
		virtual unsigned int num_obj(int dim) const = 0;

		/* number of corners of reference element */
		virtual unsigned int num_corners() const = 0;

		/* number of edges of reference element */
		virtual unsigned int num_edges() const = 0;

		/* number of faces of reference element */
		virtual unsigned int num_faces() const = 0;

		/* number of volumes of reference element */
		virtual unsigned int num_volumes() const = 0;

		/*** NUMBER OF GEOMETRIC OBJECTS OF SUB GEOMETRIC ELEMENTS ***/
		/* number of Geometric Objects of dimension dim_i of GeomObject with dim_j, nr j (dim_i < dim_j) */
		virtual unsigned int num_obj_of_obj(int dim_i, int dim_j, int j) const = 0;

		virtual unsigned int num_corners_of_edge(int j) const = 0;
		virtual unsigned int num_corners_of_face(int j) const = 0;
		virtual unsigned int num_corners_of_volume(int j) const = 0;

		virtual unsigned int num_edges_of_face(int j) const = 0;
		virtual unsigned int num_edges_of_volume(int j) const = 0;

		virtual unsigned int num_faces_of_volume(int j) const = 0;

		/*** ID's OF GEOMETRIC OBJECTS OF REFERENCE ELEMENT ***/
		/* Id of GeometricObjects of dimension i, nr. i, of Geometric object of dimension j, nr j */
		virtual int id(int dim_i, int i, int dim_j, int j) const = 0;

		/* Id of Corner i of Edge j */
		virtual int id_corner_of_edge(int i, int j) const = 0;

		/* Id of Corner i of Face j */
		virtual int id_corner_of_face(int i, int j) const = 0;

		/* Id of Corner i of Volume j */
		virtual int id_corner_of_volume(int i, int j) const = 0;

		/* Id of Edge i of Face j */
		virtual int id_edge_of_face(int i, int j) const = 0;

		/* Id of Edge i of Volume j */
		virtual int id_edge_of_volume(int i, int j) const = 0;

		/* Id of Face i of Volume j */
		virtual int id_face_of_volume(int i, int j) const = 0;

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
		virtual ~BaseReferenceElement()
		{}

		bool printInfo()
		{
			std::string GeomObjects[4] ={"Corner", "Edge", "Face", "Volume"};

			std::cout << "Reference Element Info: " << std::endl;
			std::cout << "----------------------- " << std::endl;

			std::cout << "Size: " << this->size() << std::endl;
			std::cout << "Dimension where Reference Element lives: " << this->dimension() << std::endl;
			std::cout << "Reference Corners: ";
			for(unsigned int i = 0; i< this->num_corners(); i++)
			{
				std::cout << "(" << this->corner(i)[0] << ","<< this->corner(i)[1] <<") ";
			}
			std::cout << std::endl;

			for(int i = this->dimension(); i>=0 ;i--)
			{
				std::cout << "Number of " << GeomObjects[i] << "s: " << this->num_obj(i) << std::endl;
			}

			for(int dim_i = this->dimension(); dim_i>=0 ;dim_i--)
			{
				for(unsigned int i=0; i < this->num_obj(dim_i); i++)
				{
					std::cout << GeomObjects[dim_i] << " with id '" << i << "' contains the following GeomObjects:" << std::endl;
					for(int dim_j=dim_i; dim_j>= 0; dim_j--)
					{
						std::cout << this->num_obj_of_obj(dim_i,i,dim_j) << " " << GeomObjects[dim_j] << "s with id: ";
						for(unsigned int j=0; j< this->num_obj_of_obj(dim_i,i,dim_j); j++)
						{
							std::cout << this->id(dim_i,i,dim_j,j) << " ";
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
class ReferenceElement{

private:
	ReferenceElement()
	{};
};

template <>
class ReferenceElement<Triangle>{
	private:
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME = 3, MAXDIM};
		enum{MYDIM = FACE};
		enum{MAXOBJECTS = 3};


	public:
		static const size_t Dim = FACE;


		ReferenceElement()
		{
			initializeArrays();
		}

		/* Dimension where reference element lives */
		int dimension() const
		{
			return _dimension;
		}

		/* size of reference triangle */
		double size() const
		{
			return _size;
		}

		/* coordinates of reference corner Nr i (i=0..numberOfCorners) */
		const vector2& corner(int i) const
		{
			return _corner[i];
		}

		unsigned int num_obj(int dim) const
		{
			return _num_obj[dim];
		}

		/* number of corners of reference triangle */
		unsigned int num_corners() const
		{
			return _num_obj[POINT];
		}

		unsigned int num_edges() const
		{
			return _num_obj[EDGE];
		}

		unsigned int num_faces() const
		{
			return _num_obj[FACE];
		}

		unsigned int num_volumes() const
		{
			return _num_obj[VOLUME];
		}

		unsigned int num_obj_of_obj(int dim_i, int i, int dim_j) const
		{
			return _num_obj_of_obj[dim_i][i][dim_j];
		}

		unsigned int num_corners_of_edge(int j) const
		{
			return _num_obj_of_obj[EDGE][j][POINT];
		}
		unsigned int num_corners_of_face(int j) const
		{
			return _num_obj_of_obj[FACE][j][POINT];
		}
		unsigned int num_corners_of_volume(int j) const
		{
			return 0;
		}
		unsigned int num_edges_of_face(int j) const
		{
			return _num_obj_of_obj[FACE][j][EDGE];
		}
		unsigned int num_edges_of_volume(int j) const
		{
			return 0;
		}
		unsigned int num_faces_of_volume(int j) const
		{
			return 0;
		}

		int id(int dim_i, int i, int dim_j, int j) const
		{
			return _id[dim_i][i][dim_j][j];
		}

		int id_corner_of_edge(int i, int j) const
		{
			return _id[EDGE][i][POINT][j];
		}

		int id_corner_of_face(int i, int j) const
		{
			return _id[FACE][i][POINT][j];
		}

		int id_corner_of_volume(int i, int j) const
		{
			return -1;
		}

		int id_edge_of_face(int i, int j) const
		{
			return _id[FACE][i][EDGE][j];
		}

		int id_edge_of_volume(int i, int j) const
		{
			return -1;
		}

		int id_face_of_volume(int i, int j) const
		{
			return -1;
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<Dim>& Local, MathVector<d>& Global) const
		{
			assert(d>=Dim && "ERROR in mapLocalToGlobal: d should be greater than Dim");

			MathVector<d> a10, a20;

			VecSubtract(a10, GlobalCorners[1], GlobalCorners[0]);
			VecSubtract(a20, GlobalCorners[2], GlobalCorners[0]);

			VecScale(a10, a10, Local[0]);
			VecScale(a20, a20, Local[1]);

			Global = GlobalCorners[0];
			VecAdd(Global, Global, a10);
			VecAdd(Global, Global, a20);

			return true;
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<Dim> Local[], MathVector<d> Global[], int n) const
		{
			assert(d>=Dim && "ERROR in mapLocalToGlobal: d should be greater than Dim");

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

		bool Trafo(const MathVector<1> GlobalCorners[], const MathVector<Dim>& Local, MathMatrix<1,Dim>& InvTrafo, number& det) const
		{
			assert(0 && "ERROR in Trafo: Mapping form dim 'd' to dim '1' invalid");
		}

		bool Trafo(const MathVector<2> GlobalCorners[], const MathVector<Dim>& Local, MathMatrix<2,Dim>& InvTrafo, number& det) const
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

		bool Trafo(const MathVector<3> GlobalCorners[], const MathVector<Dim>& Local, MathMatrix<3,Dim>& InvTrafo, number& det) const
		{

			MathMatrix<Dim,3> Trafo;

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

		~ReferenceElement()
		{}

	private:
		/* Dimension, in which the reference triangle lives */
		unsigned int _dimension;
		/* size of reference triangle */
		double _size;
		/* number of Geometric Objects of Reference Element
		 * (_num_obj[dim] = number of GeomObjects of dimension dim) */
		unsigned int _num_obj[MAXDIM];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int _num_obj_of_obj[MYDIM+1][MAXOBJECTS][MYDIM+1];
		/* coordinates of Reference Corner */
		vector2 _corner[3];
		// indices of GeomObjects
		int _id[MYDIM+1][MAXOBJECTS][MYDIM+1][MAXOBJECTS];

		void initializeArrays()
		{
			// dimension, where reference triangle lives
		 	_dimension = MYDIM;

			// size of reference element
		 	_size = 1./2.;

			//number of Geometric Objects
		 	_num_obj[POINT] = 3;
		 	_num_obj[EDGE] = 3;
		 	_num_obj[FACE] = 1;
		 	_num_obj[VOLUME] = 0;

			// number of Geometric Objects
		 	_num_obj_of_obj[FACE][0][POINT] = 3;
		 	_num_obj_of_obj[FACE][0][EDGE] = 3;
		 	_num_obj_of_obj[FACE][0][FACE] = 1;

		 	for(unsigned int i = 0; i < _num_obj[EDGE]; ++i)
		 	{
			 	_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	_num_obj_of_obj[EDGE][i][FACE] = 1;
		 	}

		 	for(unsigned int i = 0; i < _num_obj[POINT]; ++i)
		 	{
			 	_num_obj_of_obj[POINT][i][POINT] = 1;
			 	_num_obj_of_obj[POINT][i][EDGE] = 2;
			 	_num_obj_of_obj[POINT][i][FACE] = 1;
		 	}

			//reset _id to -1
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(unsigned int k=0; k<=MYDIM; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
						 	_id[i][j][k][l] = -1;
						}

			//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<_num_obj[i]; ++j)
				{
				 	_id[i][j][i][0] = j;
				}

			//Edges <-> Face
			for(unsigned int i=0; i<_num_obj[EDGE]; ++i)
			{
			 	_id[FACE][0][EDGE][i] = i;
			 	_id[EDGE][i][FACE][0] = 0;
			}

			// Points <-> Face
			for(unsigned int i=0; i<_num_obj[POINT]; ++i)
			{
			 	_id[FACE][0][POINT][i] = i;
			 	_id[POINT][i][FACE][0] = 0;
			}

			// Points of Edges
			// edge 0 = (0,1)
		 	_id[EDGE][0][POINT][0] = 0;
		 	_id[EDGE][0][POINT][1] = 1;
		 	// edge 1 = (1,2)
		 	_id[EDGE][1][POINT][0] = 1;
		 	_id[EDGE][1][POINT][1] = 2;
			// edge 2 = (2,0)
		 	_id[EDGE][2][POINT][0] = 2;
		 	_id[EDGE][2][POINT][1] = 0;

		 	// Edges of Point
		 	_id[POINT][0][EDGE][0] = 2;
		 	_id[POINT][0][EDGE][1] = 0;

		 	_id[POINT][1][EDGE][0] = 0;
		 	_id[POINT][1][EDGE][1] = 1;

		 	_id[POINT][2][EDGE][0] = 1;
		 	_id[POINT][2][EDGE][1] = 2;


			// Reference Corners
		 	_corner[0].x = 0.0;
		 	_corner[0].y = 0.0;
		 	_corner[1].x = 1.0;
		 	_corner[1].y = 0.0;
		 	_corner[2].x = 0.0;
		 	_corner[2].y = 1.0;
		}


};


template <>
class ReferenceElement<Quadrilateral>{
	private:
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME = 3, MAXDIM};
		enum{MYDIM = FACE};
		enum{MAXOBJECTS = 4};


	public:
		static const size_t Dim = FACE;


		ReferenceElement()
		{
			initializeArrays();
		}

		/* Dimension where reference element lives */
		int dimension() const
		{
			return _dimension;
		}

		/* size of reference triangle */
		double size() const
		{
			return _size;
		}

		/* coordinates of reference corner Nr i (i=0..numberOfCorners) */
		const vector2& corner(int i) const
		{
			return _corner[i];
		}

		unsigned int num_obj(int dim) const
		{
			return _num_obj[dim];
		}

		/* number of corners of reference triangle */
		unsigned int num_corners() const
		{
			return _num_obj[POINT];
		}

		unsigned int num_edges() const
		{
			return _num_obj[EDGE];
		}

		unsigned int num_faces() const
		{
			return _num_obj[FACE];
		}

		unsigned int num_volumes() const
		{
			return _num_obj[VOLUME];
		}

		unsigned int num_obj_of_obj(int dim_i, int i, const int dim_j) const
		{
			return _num_obj_of_obj[dim_i][i][dim_j];
		}

		unsigned int num_corners_of_edge(int j) const
		{
			return _num_obj_of_obj[EDGE][j][POINT];
		}
		unsigned int num_corners_of_face(int j) const
		{
			return _num_obj_of_obj[FACE][j][POINT];
		}
		unsigned int num_corners_of_volume(int j) const
		{
			return 0;
		}
		unsigned int num_edges_of_face(int j) const
		{
			return _num_obj_of_obj[FACE][j][EDGE];
		}
		unsigned int num_edges_of_volume(int j) const
		{
			return 0;
		}
		unsigned int num_faces_of_volume(int j) const
		{
			return 0;
		}

		int id(int dim_i, int i, int dim_j, int j) const
		{
			return _id[dim_i][i][dim_j][j];
		}

		int id_corner_of_edge(int i, int j) const
		{
			return _id[EDGE][i][POINT][j];
		}

		int id_corner_of_face(int i, int j) const
		{
			return _id[FACE][i][POINT][j];
		}

		int id_corner_of_volume(int i, int j) const
		{
			return -1;
		}

		int id_edge_of_face(int i, int j) const
		{
			return _id[FACE][i][EDGE][j];
		}

		int id_edge_of_volume(int i, int j) const
		{
			return -1;
		}

		int id_face_of_volume(int i, int j) const
		{
			return -1;
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<Dim>& Local, MathVector<d>& Global) const
		{
			assert(d >= Dim && "ERROR in mapLocalToGlobal.");

			VecScaleAdd(Global, (1.-Local[0])*(1.-Local[1]), GlobalCorners[0],
								Local[0]*(1.-Local[1])     , GlobalCorners[1],
								Local[0]*Local[1]          , GlobalCorners[2],
								(1.-Local[0])*Local[1]     , GlobalCorners[3]);

			return true;
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<Dim> Local[], MathVector<d> Global[], int n) const
		{
			assert(d >= Dim && "ERROR in mapLocalToGlobal.");
			for(int i=0; i<n; ++i)
			{
				VecScaleAdd(Global[i], (1.-Local[i][0])*(1.-Local[i][1]), GlobalCorners[0],
									   Local[i][0]*(1.-Local[i][1])     , GlobalCorners[1],
									   Local[i][0]*Local[i][1]          , GlobalCorners[2],
									   (1.-Local[i][0])*Local[i][1]     , GlobalCorners[3]);
			}
			return true;
		}

		bool Trafo(const MathVector<1> GlobalCorners[], const MathVector<Dim>& Local, MathMatrix<1,Dim>& InvTrafo, number& det) const
		{
			assert(0 && "ERROR in Trafo: Mapping form dim 'd' to dim '1' invalid");
			return false;
		}
		bool Trafo(const MathVector<2> GlobalCorners[], const MathVector<Dim>& Local, MathMatrix<2,Dim>& InvTrafo, number& det) const
		{
			number a;
			MathMatrix<2,2> Trafo;

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
			return false;
		}

		bool Trafo(const MathVector<3> GlobalCorners[], const MathVector<Dim> Local, MathMatrix<3,Dim>& InvTrafo, number& det) const
		{
			assert(0 && "Not implemented yet.");
			return false;
		}

		~ReferenceElement()
		{}

	private:
		/* Dimension, in which the reference triangle lives */
		unsigned int _dimension;
		/* size of reference triangle */
		double _size;
		/* number of Geometric Objects of Reference Element
		 * (_num_obj[dim] = number of GeomObjects of dimension dim) */
		unsigned int _num_obj[MAXDIM];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int _num_obj_of_obj[MYDIM+1][MAXOBJECTS][MYDIM+1];
		/* coordinates of Reference Corner */
		MathVector<Dim> _corner[4];
		// indices of GeomObjects
		int _id[MYDIM+1][MAXOBJECTS][MYDIM+1][MAXOBJECTS];

		void initializeArrays()
		{
			// dimension, where reference triangle lives
		 	_dimension = MYDIM;

			// size of reference element
		 	_size = 1.;

			//number of Geometric Objects
		 	_num_obj[POINT] = 4;
		 	_num_obj[EDGE] = 4;
		 	_num_obj[FACE] = 1;
		 	_num_obj[VOLUME] = 0;

			// number of Geometric Objects
		 	_num_obj_of_obj[FACE][0][POINT] = 4;
		 	_num_obj_of_obj[FACE][0][EDGE] = 4;
		 	_num_obj_of_obj[FACE][0][FACE] = 1;

		 	for(unsigned int i = 0; i < _num_obj[EDGE]; ++i)
		 	{
		 		_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	_num_obj_of_obj[EDGE][i][FACE] = 1;
		 	}

		 	for(unsigned int i = 0; i < _num_obj[EDGE]; ++i)
		 	{
		 		_num_obj_of_obj[POINT][i][POINT] = 1;
		 		_num_obj_of_obj[POINT][i][EDGE] = 2;
		 		_num_obj_of_obj[POINT][i][FACE] = 1;
		 	}

			//reset _id to -1
			for(unsigned int i=0; i<=MYDIM; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(unsigned int k=0; k<=MYDIM; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
						 	_id[i][j][k][l] = -1;
						}

			//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
			for(unsigned int d=0; d<=MYDIM; ++d)
				for(unsigned int j=0; j<_num_obj[d]; ++j)
				{
				 	_id[d][j][d][0] = j;
				}

			//Edges <-> Face
			for(unsigned int i=0; i<_num_obj[EDGE]; ++i)
			{
			 	_id[FACE][0][EDGE][i] = i;
			 	_id[EDGE][i][FACE][0] = 0;
			}

			// Points <-> Face
			for(unsigned int i=0; i<_num_obj[POINT]; ++i)
			{
			 	_id[FACE][0][POINT][i] = i;
			 	_id[POINT][i][FACE][0] = 0;
			}

			// Points of Edges
			// edge 0 = (0,1)
		 	_id[EDGE][0][POINT][0] = 0;
		 	_id[EDGE][0][POINT][1] = 1;
			// edge 1 = (1,2)
		 	_id[EDGE][1][POINT][0] = 1;
		 	_id[EDGE][1][POINT][1] = 2;
			// edge 2 = (2,3)
		 	_id[EDGE][2][POINT][0] = 2;
		 	_id[EDGE][2][POINT][1] = 3;
			// edge 3 = (3,0)
		 	_id[EDGE][3][POINT][0] = 3;
		 	_id[EDGE][3][POINT][1] = 0;

		 	// Edges of Point
		 	_id[POINT][0][EDGE][0] = 3;
		 	_id[POINT][0][EDGE][1] = 0;

		 	_id[POINT][1][EDGE][0] = 0;
		 	_id[POINT][1][EDGE][1] = 1;

		 	_id[POINT][2][EDGE][0] = 1;
		 	_id[POINT][2][EDGE][1] = 2;

		 	_id[POINT][3][EDGE][0] = 2;
		 	_id[POINT][3][EDGE][1] = 3;

			// Reference Corners
		 	_corner[0].x = 0.0;
		 	_corner[0].y = 0.0;
		 	_corner[1].x = 1.0;
		 	_corner[1].y = 0.0;
		 	_corner[2].x = 1.0;
		 	_corner[2].y = 1.0;
		 	_corner[3].x = 0.0;
		 	_corner[3].y = 1.0;
		}

};



// Singleton, holding all Reference Elements available
template <typename TElem>
class ReferenceElements {

	private:

		static ReferenceElements& inst()
		{
			static ReferenceElements myInst;
			return myInst;
		};

		// private constructor
		ReferenceElements()
		{};

		inline static const ReferenceElement<TElem>& get_ReferenceElement()
		{
			static ug::ReferenceElement<TElem> refElem;
			return refElem;
		}

	public:
		static const ReferenceElement<TElem>& ReferenceElement()
		{
			return inst().get_ReferenceElement();
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

		const static int local_dim = 2;
		const static int num_corners = 3;
		const static int num_edges = 3;
		const static int num_faces = 1;
		const static int num_volumes = 0;
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

		const static int local_dim = 2;
		const static int num_corners = 4;
		const static int num_edges = 4;
		const static int num_faces = 1;
		const static int num_volumes = 0;

		enum
		{
			RefDim = 2,
			NumberCorners = 4
		};
};


}

#endif /* __H__LIBDISCRETIZATION__REFERENCEELEMENT__ */
