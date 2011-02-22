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
#include "lib_discretization/spatial_discretization/disc_helper/disc_helper.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"

namespace ug{

template<template <	class TElem, int TWorldDim> class TFVGeom,
					typename TDomain,
					typename TAlgebra>
class FVNeumannBoundaryElemDisc
	: public IElemDisc<TAlgebra>
{
	public:
		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = TDomain::dim;

		// position type
		typedef typename TDomain::position_type position_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef LocalMatrix<typename TAlgebra::matrix_type::value_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::value_type> local_vector_type;

		// local index type
		//typedef typename algebra_type::vector_type::local_index_type local_index_type;
		typedef LocalIndices local_index_type;

	protected:
		typedef typename IBoundaryNumberProvider<dim>::functor_type BNDNumberFunctor;

	public:
		FVNeumannBoundaryElemDisc()
		 : m_numFct(0), m_pDomain(NULL)
			{
				m_mBoundarySegment.clear();
				register_assemble_functions(Int2Type<dim>());
			}

		void set_domain(domain_type& domain) {m_pDomain = &domain;}

		bool add_boundary_value(IBoundaryNumberProvider<dim>& user, const char* function, const char* subsets)
		{
		//	check that function pattern exists
			if(this->m_pPattern == NULL)
			{
				UG_LOG("FVNeumannBoundaryElemDisc:add_boundary_value: Function Pattern not set.\n");
				return false;
			}

		//	create Function Group and Subset Group
			FunctionGroup functionGroup;
			SubsetGroup subsetGroup;

		//	convert strings
			if(!ConvertStringToSubsetGroup(subsetGroup, *this->m_pPattern, subsets))
			{
				UG_LOG("ERROR while parsing Subsets.\n");
				return false;
			}
			if(!ConvertStringToFunctionGroup(functionGroup, *this->m_pPattern, function))
			{
				UG_LOG("ERROR while parsing Functions.\n");
				return false;
			}

		//	only one function allowed
			if(functionGroup.num_fct() != 1)
			{
				UG_LOG("FVNeumannBoundaryElemDisc:add_boundary_value: Exactly one function needed, but given '"<<function<<"' as functions.\n");
				return false;
			}

		//	forward request
			return add_boundary_value(user, functionGroup.unique_id(0), subsetGroup);
		}

		bool add_boundary_value(IBoundaryNumberProvider<dim>& user, size_t fct, SubsetGroup bndSubsetGroup)
		{
		//	check that function pattern exists
			if(this->m_pPattern == NULL)
			{
				UG_LOG("FVNeumannBoundaryElemDisc:add_boundary_value: Function Pattern not set.\n");
				return false;
			}

		//	get subsethandler
			const ISubsetHandler* pSH = this->m_pPattern->get_subset_handler();

		// 	check if function exist
			if(fct >= this->m_pPattern->num_fct())
			{
				UG_LOG("FVNeumannBoundaryElemDisc:add_boundary_value: Function "
						<< fct << " does not exist in pattern.\n");
				return false;
			}

		//	add function to function group if needed
		//	this will force, that the local vector contains the function
		//	note: 	The local indices for allready contained functions are not changed.
		//			Therefore an update in the UserDataFunctions is not necessary
			if(!this->m_FunctionGroup.contains(fct))
			{
				this->m_FunctionGroup.add(fct);
				m_numFct = this->m_FunctionGroup.num_fct();
			}

		//	get position of function in function group
			const size_t index = this->m_FunctionGroup.local_index(fct);

		// 	check that function is defined on inner subset
			if(!this->m_SubsetGroup.empty())
			{
				for(size_t si = 0; si < this->m_SubsetGroup.num_subsets(); ++si)
				{
					const int subsetIndex = this->m_SubsetGroup[si];
					if(!this->m_pPattern->is_def_in_subset(fct, subsetIndex))
					{
						UG_LOG("FVNeumannBoundaryElemDisc:add_boundary_value: Function "
								<< fct << " not defined in subset " << subsetIndex << ".\n");
						return false;
					}
				}
			}

		// 	loop subsets
			for(size_t si = 0; si < bndSubsetGroup.num_subsets(); ++si)
			{
			//	get subset index
				const int subsetIndex = bndSubsetGroup[si];

			//	check that subsetIndex is valid
				if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
				{
					UG_LOG("FVNeumannBoundaryElemDisc:add_boundary_value: Invalid subset Index "
							<< subsetIndex << ". (Valid is 0, .. , " << pSH->num_subsets() <<").\n");
					return false;
				}

			//	get Boundary segment from map
				std::vector<UserDataFunction>& vSegmentFunction = m_mBoundarySegment[subsetIndex];

			//	remember functor and function
				vSegmentFunction.push_back(UserDataFunction(index, user.get_functor()));
			}

		//	we're done
			return true;
		}

	private:
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
		virtual size_t num_fct(){return m_numFct;}

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
		template <typename TElem>
		inline bool prepare_element_loop();

		template <typename TElem>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		template <typename TElem>
		inline bool finish_element_loop();

		template <typename TElem>
		inline bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_f(local_vector_type& d, number time=0.0);

	private:
		// domain
		TDomain* m_pDomain;

		// position access
		std::vector<position_type> m_vCornerCoords;
		typename TDomain::position_accessor_type m_aaPos;

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
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_prepare_element_loop_function(	id, &FVNeumannBoundaryElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &FVNeumannBoundaryElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &FVNeumannBoundaryElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &FVNeumannBoundaryElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &FVNeumannBoundaryElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &FVNeumannBoundaryElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &FVNeumannBoundaryElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &FVNeumannBoundaryElemDisc::template assemble_f<TElem>);
		}

};

}

#include "neumann_boundary_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__*/
