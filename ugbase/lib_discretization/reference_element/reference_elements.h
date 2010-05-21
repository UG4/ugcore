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
#include "lib_grid/lib_grid.h"

namespace ug{

enum ReferenceElementType {
	RET_INVALID = -1,
	RET_POINT = 0,
	RET_EDGE,
	RET_TRIANGLE,
	RET_QUADRILATERAL,
	NUM_REFERENCE_ELEMENTS
};

template <int d>
class BaseReferenceElement {

	public:
		static const int dim = d;

		/*** GENERAL INFORMATIONS ***/
		/* dimension in which reference element lives */
		int dimension() const {return dim;};

		/* size of reference element */
		virtual double size() const = 0;

		/* coordinates of corner of reference element (i=0..numberOfCorners) */
		virtual const MathVector<dim>& corner(int i) const = 0;

		/* number of reference elements with id of dimension 'dim' of this reference element */
		virtual unsigned int num_ref_elem(ReferenceElementType type) const = 0;

		/* reference element id of subObject nr i of dimension dim_i */
		virtual ReferenceElementType ref_elem_type(int dim_i, int i) const = 0;

		/* number of geom objects of dimension 'dim' of this reference element */
		virtual unsigned int num_obj(int dim) const = 0;

		/*** NUMBER OF GEOMETRIC OBJECTS OF SUB GEOMETRIC ELEMENTS ***/
		/* number of Geometric Objects of dimension dim_i of GeomObject with dim_j, nr j (dim_i < dim_j) */
		virtual unsigned int num_obj_of_obj(int dim_i, int dim_j, int j) const = 0;

		/*** ID's OF GEOMETRIC OBJECTS OF REFERENCE ELEMENT ***/
		/* Id of GeometricObjects of dimension i, nr. i, of Geometric object of dimension j, nr j */
		virtual int id(int dim_i, int i, int dim_j, int j) const = 0;

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

class ReferenceTriangle{
	public:
		static const int dim = 2;


		ReferenceTriangle()
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
		const MathVector<dim>& corner(int i) const
		{
			return _corner[i];
		}

		/* number of reference elements with id of dimension 'dim' of this reference element */
		unsigned int num_ref_elem(ReferenceElementType type) const
		{
			return m_ref_elem[type];
		}

		/* reference element type of subObject nr i of dimension dim_i */
		ReferenceElementType ref_elem_type(int dim_i, int i) const
		{
			return m_ref_elem_type[dim_i][i];
		}

		unsigned int num_obj(int dim) const
		{
			return _num_obj[dim];
		}

		unsigned int num_obj_of_obj(int dim_i, int i, int dim_j) const
		{
			return _num_obj_of_obj[dim_i][i][dim_j];
		}

		int id(int dim_i, int i, int dim_j, int j) const
		{
			return _id[dim_i][i][dim_j][j];
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<dim>& Local, MathVector<d>& Global) const
		{
			assert(d>=dim && "ERROR in mapLocalToGlobal: d should be greater than dim");

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
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<dim> Local[], MathVector<d> Global[], int n) const
		{
			assert(d>=dim && "ERROR in mapLocalToGlobal: d should be greater than dim");

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

		bool Trafo(const MathVector<1> GlobalCorners[], const MathVector<dim>& Local, MathMatrix<1,dim>& InvTrafo, number& det) const
		{
			assert(0 && "ERROR in Trafo: Mapping form dim 'd' to dim '1' invalid");
			return false;
		}

		bool Trafo(const MathVector<2> GlobalCorners[], const MathVector<dim>& Local, MathMatrix<2,dim>& InvTrafo, number& det) const
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

		bool Trafo(const MathVector<3> GlobalCorners[], const MathVector<dim>& Local, MathMatrix<3,dim>& InvTrafo, number& det) const
		{

			MathMatrix<dim,3> Trafo;

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

		~ReferenceTriangle()
		{}

	private:
		// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2};
		enum{MAXOBJECTS = 3};


		/* Dimension, in which the reference triangle lives */
		unsigned int _dimension;
		/* size of reference triangle */
		double _size;
		/* number of Geometric Objects of Reference Element
		 * (_num_obj[dim] = number of GeomObjects of dimension dim) */
		unsigned int _num_obj[dim+1];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int _num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
		/* coordinates of Reference Corner */
		vector2 _corner[3];
		// indices of GeomObjects
		int _id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

		unsigned int m_ref_elem[NUM_REFERENCE_ELEMENTS];
		ReferenceElementType m_ref_elem_type[dim+1][MAXOBJECTS];

		void initializeArrays()
		{
			// dimension, where reference triangle lives
		 	_dimension = dim;

			// size of reference element
		 	_size = 1./2.;

			//number of Geometric Objects
		 	_num_obj[POINT] = 3;
		 	_num_obj[EDGE] = 3;
		 	_num_obj[FACE] = 1;

			// number of Geometric Objects
		 	_num_obj_of_obj[FACE][0][POINT] = 3;
		 	_num_obj_of_obj[FACE][0][EDGE] = 3;
		 	_num_obj_of_obj[FACE][0][FACE] = 1;

		 	m_ref_elem_type[FACE][0] = RET_TRIANGLE;

		 	for(unsigned int i = 0; i < _num_obj[EDGE]; ++i)
		 	{
			 	_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	_num_obj_of_obj[EDGE][i][FACE] = 1;

			 	m_ref_elem_type[EDGE][i] = RET_EDGE;
		 	}

		 	for(unsigned int i = 0; i < _num_obj[POINT]; ++i)
		 	{
			 	_num_obj_of_obj[POINT][i][POINT] = 1;
			 	_num_obj_of_obj[POINT][i][EDGE] = 2;
			 	_num_obj_of_obj[POINT][i][FACE] = 1;

			 	m_ref_elem_type[POINT][i] = RET_POINT;
		 	}

			//reset _id to -1
			for(unsigned int i=0; i<=_dimension; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(unsigned int k=0; k<=_dimension; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
						 	_id[i][j][k][l] = -1;
						}

			//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
			for(unsigned int i=0; i<=_dimension; ++i)
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

		 	// Reference Element Types
		 	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
		 	{
				m_ref_elem[i] = 0;
		 	}
		 	m_ref_elem[RET_POINT] = 3;
		 	m_ref_elem[RET_EDGE] = 3;
		 	m_ref_elem[RET_TRIANGLE] = 1;
		}


};


class ReferenceQuadrilateral{
	public:
		// own dimension
		static const int dim = 2;


		ReferenceQuadrilateral()
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
		const MathVector<dim>& corner(int i) const
		{
			return _corner[i];
		}

		/* number of reference elements with id of dimension 'dim' of this reference element */
		unsigned int num_ref_elem(ReferenceElementType id) const
		{
			return m_ref_elem[id];
		}

		/* reference element type of subObject nr i of dimension dim_i */
		ReferenceElementType ref_elem_type(int dim_i, int i) const
		{
			return m_ref_elem_type[dim_i][i];
		}

		unsigned int num_obj(int dim) const
		{
			return _num_obj[dim];
		}

		unsigned int num_obj_of_obj(int dim_i, int i, const int dim_j) const
		{
			return _num_obj_of_obj[dim_i][i][dim_j];
		}

		int id(int dim_i, int i, int dim_j, int j) const
		{
			return _id[dim_i][i][dim_j][j];
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<dim>& Local, MathVector<d>& Global) const
		{
			assert(d >= dim && "ERROR in mapLocalToGlobal.");

			VecScaleAdd(Global, (1.-Local[0])*(1.-Local[1]), GlobalCorners[0],
								Local[0]*(1.-Local[1])     , GlobalCorners[1],
								Local[0]*Local[1]          , GlobalCorners[2],
								(1.-Local[0])*Local[1]     , GlobalCorners[3]);

			return true;
		}

		template <int d>
		bool mapLocalToGlobal(const MathVector<d> GlobalCorners[], const MathVector<dim> Local[], MathVector<d> Global[], int n) const
		{
			assert(d >= dim && "ERROR in mapLocalToGlobal.");
			for(int i=0; i<n; ++i)
			{
				VecScaleAdd(Global[i], (1.-Local[i][0])*(1.-Local[i][1]), GlobalCorners[0],
									   Local[i][0]*(1.-Local[i][1])     , GlobalCorners[1],
									   Local[i][0]*Local[i][1]          , GlobalCorners[2],
									   (1.-Local[i][0])*Local[i][1]     , GlobalCorners[3]);
			}
			return true;
		}

		bool Trafo(const MathVector<1> GlobalCorners[], const MathVector<dim>& Local, MathMatrix<1,dim>& InvTrafo, number& det) const
		{
			assert(0 && "ERROR in Trafo: Mapping form dim 'd' to dim '1' invalid");
			return false;
		}
		bool Trafo(const MathVector<2> GlobalCorners[], const MathVector<dim>& Local, MathMatrix<2,dim>& InvTrafo, number& det) const
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

		bool Trafo(const MathVector<3> GlobalCorners[], const MathVector<dim> Local, MathMatrix<3,dim>& InvTrafo, number& det) const
		{
			assert(0 && "Not implemented yet.");
			return false;
		}

		~ReferenceQuadrilateral()
		{}

	private:
		// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2};
		enum{MAXOBJECTS = 4};

		/* Dimension, in which the reference triangle lives */
		unsigned int _dimension;
		/* size of reference triangle */
		double _size;
		/* number of Geometric Objects of Reference Element
		 * (_num_obj[dim] = number of GeomObjects of dimension dim) */
		unsigned int _num_obj[dim+1];
		/* number of Geometric Objects contained in a (Sub-)Geometric Object of the Element */
		unsigned int _num_obj_of_obj[dim+1][MAXOBJECTS][dim+1];
		/* coordinates of Reference Corner */
		MathVector<dim> _corner[4];
		// indices of GeomObjects
		int _id[dim+1][MAXOBJECTS][dim+1][MAXOBJECTS];

		unsigned int m_ref_elem[NUM_REFERENCE_ELEMENTS];
		ReferenceElementType m_ref_elem_type[dim+1][MAXOBJECTS];

		void initializeArrays()
		{
			// dimension, where reference triangle lives
		 	_dimension = dim;

			// size of reference element
		 	_size = 1.;

			//number of Geometric Objects
		 	_num_obj[POINT] = 4;
		 	_num_obj[EDGE] = 4;
		 	_num_obj[FACE] = 1;

			// number of Geometric Objects
		 	_num_obj_of_obj[FACE][0][POINT] = 4;
		 	_num_obj_of_obj[FACE][0][EDGE] = 4;
		 	_num_obj_of_obj[FACE][0][FACE] = 1;
		 	m_ref_elem_type[FACE][0] = RET_QUADRILATERAL;

		 	for(unsigned int i = 0; i < _num_obj[EDGE]; ++i)
		 	{
		 		_num_obj_of_obj[EDGE][i][EDGE] = 1;
			 	_num_obj_of_obj[EDGE][i][POINT] = 2;
			 	_num_obj_of_obj[EDGE][i][FACE] = 1;

			 	m_ref_elem_type[EDGE][i] = RET_EDGE;
		 	}

		 	for(unsigned int i = 0; i < _num_obj[EDGE]; ++i)
		 	{
		 		_num_obj_of_obj[POINT][i][POINT] = 1;
		 		_num_obj_of_obj[POINT][i][EDGE] = 2;
		 		_num_obj_of_obj[POINT][i][FACE] = 1;

			 	m_ref_elem_type[POINT][i] = RET_POINT;
		 	}

			//reset _id to -1
			for(unsigned int i=0; i<=_dimension; ++i)
				for(unsigned int j=0; j<MAXOBJECTS; ++j)
					for(unsigned int k=0; k<=_dimension; ++k)
						for(unsigned int l=0; l<MAXOBJECTS; l++)
						{
						 	_id[i][j][k][l] = -1;
						}

			//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
			for(unsigned int d=0; d<=_dimension; ++d)
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

		 	// Reference Element Types
		 	for(int i = 0; i < NUM_REFERENCE_ELEMENTS; ++i)
		 	{
				m_ref_elem[i] = 0;
		 	}
		 	m_ref_elem[RET_POINT] = 4;
		 	m_ref_elem[RET_EDGE] = 4;
		 	m_ref_elem[RET_QUADRILATERAL] = 1;
		}

};

class ReferenceVertex
{
	// TODO: Implement
};

class ReferenceEdge
{
	// TODO: Implement
};


template <class TElem>
class reference_element_traits
{};

template <>
class reference_element_traits<ReferenceVertex>
{
	public:
		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_POINT;

		static const int dim = 0;
		static const int num_corners = 1;
		static const int num_edges = 0;
		static const int num_faces = 0;
		static const int num_volumes = 0;

		typedef MathVector<dim> position_type;
};

template <>
class reference_element_traits<ReferenceEdge>
{
	public:
		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_EDGE;

		static const int dim = 1;
		static const int num_corners = 2;
		static const int num_edges = 1;
		static const int num_faces = 0;
		static const int num_volumes = 0;

		typedef MathVector<dim> position_type;
};


template <>
class reference_element_traits<ReferenceTriangle>
{
	public:
		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_TRIANGLE;

		static const int dim = 2;
		static const int num_corners = 3;
		static const int num_edges = 3;
		static const int num_faces = 1;
		static const int num_volumes = 0;

		typedef MathVector<dim> position_type;
};

template <>
class reference_element_traits<ReferenceQuadrilateral>
{
	public:
		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_QUADRILATERAL;

		static const int dim = 2;
		static const int num_corners = 4;
		static const int num_edges = 4;
		static const int num_faces = 1;
		static const int num_volumes = 0;

		typedef MathVector<dim> position_type;
};

template <>
class reference_element_traits<VertexBase>
{
	public:
		typedef ReferenceVertex reference_element_type;
		typedef reference_element_traits<ReferenceVertex> reference_element_type_traits;

		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = reference_element_type_traits::REFERENCE_ELEMENT_TYPE;
		static const int dim = reference_element_type_traits::dim;
		static const int num_corners = reference_element_type_traits::num_corners;
		static const int num_edges = reference_element_type_traits::num_edges;
		static const int num_faces = reference_element_type_traits::num_faces;
		static const int num_volumes = reference_element_type_traits::num_volumes;
		typedef reference_element_type_traits::position_type position_type;
};

template <>
class reference_element_traits<Edge>
{
	public:
		typedef ReferenceEdge reference_element_type;
		typedef reference_element_traits<ReferenceEdge> reference_element_type_traits;

		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = reference_element_type_traits::REFERENCE_ELEMENT_TYPE;
		static const int dim = reference_element_type_traits::dim;
		static const int num_corners = reference_element_type_traits::num_corners;
		static const int num_edges = reference_element_type_traits::num_edges;
		static const int num_faces = reference_element_type_traits::num_faces;
		static const int num_volumes = reference_element_type_traits::num_volumes;
		typedef reference_element_type_traits::position_type position_type;
};

template <>
class reference_element_traits<Triangle>
{
	public:
		typedef ReferenceTriangle reference_element_type;
		typedef reference_element_traits<ReferenceTriangle> reference_element_type_traits;

		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = reference_element_type_traits::REFERENCE_ELEMENT_TYPE;
		static const int dim = reference_element_type_traits::dim;
		static const int num_corners = reference_element_type_traits::num_corners;
		static const int num_edges = reference_element_type_traits::num_edges;
		static const int num_faces = reference_element_type_traits::num_faces;
		static const int num_volumes = reference_element_type_traits::num_volumes;
		typedef reference_element_type_traits::position_type position_type;
};

template <>
class reference_element_traits<Quadrilateral>
{
	public:
		typedef ReferenceQuadrilateral reference_element_type;
		typedef reference_element_traits<ReferenceQuadrilateral> reference_element_type_traits;

		static const ReferenceElementType REFERENCE_ELEMENT_TYPE = reference_element_type_traits::REFERENCE_ELEMENT_TYPE;
		static const int dim = reference_element_type_traits::dim;
		static const int num_corners = reference_element_type_traits::num_corners;
		static const int num_edges = reference_element_type_traits::num_edges;
		static const int num_faces = reference_element_type_traits::num_faces;
		static const int num_volumes = reference_element_type_traits::num_volumes;
		typedef reference_element_type_traits::position_type position_type;
};


}

#endif /* __H__LIBDISCRETIZATION__REFERENCEELEMENT__ */
