/*
 * neumann_boundary.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spatial_discretization/disc_helper/geometry_provider.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"

namespace ug{

template<typename TDomain, typename TAlgebra>
class FVNeumannBoundaryElemDisc
	: public IDomainElemDisc<TDomain, TAlgebra>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain, TAlgebra> base_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	///	Algebra type
		typedef typename base_type::algebra_type algebra_type;

	///	Local matrix type
		typedef typename base_type::local_matrix_type local_matrix_type;

	///	Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Local index type
		typedef typename base_type::local_index_type local_index_type;

	protected:
	///	type of bnd number
		typedef typename IBoundaryData<number, dim>::functor_type BNDNumberFunctor;

	public:
	///	default constructor
		FVNeumannBoundaryElemDisc()
		 : m_numFct(0)
			{
				m_mBoundarySegment.clear();
				register_assemble_functions(Int2Type<dim>());
			}

	///	add a boundary value
		bool add_boundary_value(IBoundaryData<number, dim>& user, const char* function, const char* subsets);

	///	add a boundary value
		bool add_boundary_value(IBoundaryData<number, dim>& user, size_t fct, SubsetGroup bndSubsetGroup);

	private:
	//	Functor, function grouping
		struct UserDataFunction
		{
			UserDataFunction(size_t fct_, BNDNumberFunctor functor_)
				: loc_fct(fct_), functor(functor_) {}

			size_t loc_fct;
			BNDNumberFunctor functor;
		};

		std::map<int, std::vector<UserDataFunction> > m_mBoundarySegment;
		size_t m_numFct;

	public:
	///	number of functons required
		virtual size_t num_fct(){return m_numFct;}

	///	type of function required
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}


	private:
		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool prepare_element_loop();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool finish_element_loop();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		inline bool assemble_f(local_vector_type& d, number time=0.0);

	private:
	// 	position access
		const position_type* m_vCornerCoords;

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////
		// register for 1D
		void register_assemble_functions(Int2Type<1>)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
		}

		// register for 2D
		void register_assemble_functions(Int2Type<2>)
		{
			register_assemble_functions(Int2Type<1>());
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

		// register for 3D
		void register_assemble_functions(Int2Type<3>)
		{
			register_assemble_functions(Int2Type<2>());
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Pyramid>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

		// help function
		template<typename TElem>
		void register_all_assemble_functions(int id)
		{
			typedef FVNeumannBoundaryElemDisc T;

			register_prepare_element_loop_function(	id, &T::template prepare_element_loop<TElem, FV1Geometry>);
			register_prepare_element_function(		id, &T::template prepare_element<TElem, FV1Geometry>);
			register_finish_element_loop_function(	id, &T::template finish_element_loop<TElem, FV1Geometry>);
			register_assemble_JA_function(			id, &T::template assemble_JA<TElem, FV1Geometry>);
			register_assemble_JM_function(			id, &T::template assemble_JM<TElem, FV1Geometry>);
			register_assemble_A_function(			id, &T::template assemble_A<TElem, FV1Geometry>);
			register_assemble_M_function(			id, &T::template assemble_M<TElem, FV1Geometry>);
			register_assemble_f_function(			id, &T::template assemble_f<TElem, FV1Geometry>);
		}

};

} // end namespac ug

#include "neumann_boundary_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__*/
