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
#include "lib_grid/lg_base.h"

namespace ug{

enum ReferenceElementType {
	RET_INVALID = -1,
	RET_POINT = 0,
	RET_EDGE,
	RET_TRIANGLE,
	RET_QUADRILATERAL,
	RET_TETRAHEDRON,
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



template <typename TRefElem, int TWorldDim>
class ReferenceMapping
{
	public:
		static const int world_dim = TWorldDim;
		static const int dim = TRefElem::dim;

	public:
		ReferenceMapping()
		{}

		void update(const MathVector<world_dim>* corners)
		{}

		bool local_to_global(	const MathVector<dim> loc_pos,
								MathVector<world_dim>& glob_pos) const
		{return false;}

		bool jacobian_transposed(	const MathVector<dim> loc_pos,
									MathMatrix<dim, world_dim>& JT) const
		{return false;}

		bool jacobian_transposed_inverse(	const MathVector<dim> loc_pos,
											MathMatrix<world_dim, dim>& JTInv) const
		{return false;}

		bool jacobian_det(const MathVector<dim> loc_pos, number& det) const
		{return false;}
};



class ReferenceVertex
{
public:
	static const ReferenceElementType REFERENCE_ELEMENT_TYPE = RET_POINT;

	static const int dim = 0;
	static const int num_corners = 1;
	static const int num_edges = 0;
	static const int num_faces = 0;
	static const int num_volumes = 0;

	// TODO: Implement
};


template <class TElem>
class reference_element_traits
{};

template <>
class reference_element_traits<VertexBase>
{
	public:
		typedef ReferenceVertex reference_element_type;
};

}

#include "reference_edge.h"
#include "reference_triangle.h"
#include "reference_quadrilateral.h"
#include "reference_tetrahedron.h"

#endif /* __H__LIBDISCRETIZATION__REFERENCEELEMENT__ */
