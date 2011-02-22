/*
 * neumann_boundary.h
 *
 *  Created on: 26.02.2010
 *      Author: markusbreit
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__

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

/// Finite Volume Element Discretization for an inner BndCond that depends on the unknowns (on the bnd)
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the unknown-dependent Neumann-flux over an inner boundary.
 * The equation of this flux should be given on the script level.
 * 
 * \tparam	TFVGeom		Finite Volume Geometry used
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */


template<template <	class TElem, int TWorldDim> class TFVGeom,
					typename TDomain,
					typename TAlgebra>
class FVInnerBoundaryElemDisc
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
		typedef typename IElemDisc<TAlgebra>::local_matrix_type local_matrix_type;

		// local vector type
		typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;

		// local index type
		typedef typename IElemDisc<TAlgebra>::local_index_type local_index_type;
		//typedef LocalIndices local_index_type;

	protected:
		typedef typename IBoundaryNumberProvider<dim>::functor_type BNDNumberFunctor;

	public:
		FVInnerBoundaryElemDisc()
		 : m_numFct(0), m_pDomain(NULL)
			{
				m_mBoundarySegment.clear();
				register_assemble_functions(Int2Type<dim>());
			}

		void set_domain(domain_type& domain) {m_pDomain = &domain;}

		
		bool add_flux(const char* function, const char* subsets)
		{
			//	check that function pattern exists
			if(this->m_pPattern == NULL)
			{
				UG_LOG("FVInnerBoundaryElemDisc:add_flux: Function Pattern not set.\n");
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
			
			/*
			//	only one function allowed
			if(functionGroup.num_fct() != 1)
			{
				UG_LOG("FVInnerBoundaryElemDisc:add_flux: Exactly one function needed, but given '"<<function<<"' as functions.\n");
				return false;
			}
			*/
			
		//	forward request
			return add_flux(functionGroup, subsetGroup);
		}

		bool add_flux(FunctionGroup functionGroup, SubsetGroup bndSubsetGroup)
		{
			//	check that function pattern exists
			if(this->m_pPattern == NULL)
			{
				UG_LOG("FVInnerBoundaryElemDisc:add_flux: Function Pattern not set.\n");
				return false;
			}

			//	get subsethandler
			const ISubsetHandler* pSH = this->m_pPattern->get_subset_handler();
			
			
			for (size_t fct = 0; fct < functionGroup.num_fct(); fct++)
			{
				size_t fct_id = functionGroup.unique_id(fct);
				
				// 	check if functions exist
				if(fct_id >= this->m_pPattern->num_fct())
				{
					UG_LOG("FVInnerBoundaryElemDisc:add_flux: Function "
							<< fct_id << " does not exist in pattern.\n");
					return false;
				}
				
				//	add function to function group if needed
				//	this will force, that the local vector contains the function
				//	note: 	The local indices for allready contained functions are not changed.
				//			Therefore an update in the UserDataFunctions is not necessary
				if(!this->m_FunctionGroup.contains(fct_id))
				{
					this->m_FunctionGroup.add(fct_id);
					m_numFct = this->m_FunctionGroup.num_fct();
				}

				//	get position of function in function group
				size_t index = this->m_FunctionGroup.local_index(fct_id);
				
				
				// 	check that function is defined on inner subset
				if(!this->m_SubsetGroup.empty())
				{
					for (size_t si = 0; si < this->m_SubsetGroup.num_subsets(); ++si)
					{
						const int subsetIndex = this->m_SubsetGroup[si];
						if (!this->m_pPattern->is_def_in_subset(fct_id, subsetIndex))
						{
							UG_LOG("FVInnerBoundaryElemDisc:add_flux: Function "
									<< fct_id << " not defined in subset " << subsetIndex << ".\n");
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
						UG_LOG("FVInnerBoundaryElemDisc:add_flux: Invalid subset Index "
								<< subsetIndex << ". (Valid is 0, .. , " << pSH->num_subsets() <<").\n");
						return false;
					}

					//	get Boundary segment from map
					std::vector<size_t>& vSegmentFunction = m_mBoundarySegment[subsetIndex];
					
					//	remember functor and function
					vSegmentFunction.push_back(index);
				}
			}
		//	we're done
			return true;
		}
		
		
	private:
		std::map<int, std::vector<size_t> > m_mBoundarySegment;
		size_t m_numFct;

	
	public:	// inherited from IElemDisc
	///	number of functions used
		virtual size_t num_fct(){return m_numFct;}
		
	///	type of trial space for each function used
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
				UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}

	///	returns if hanging nodes are used
		virtual bool use_hanging() const {return TFVGeom<Edge, dim>::usesHangingNodes;}

	private:
	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem>
		inline bool prepare_element_loop();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

	///	finishes the loop over all elements
		template <typename TElem>
		inline bool finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem>
		inline bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem>
		inline bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

	///	assembles the stiffness part of the local defect
		template <typename TElem>
		inline bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

	///	assembles the mass part of the local defect
		template <typename TElem>
		inline bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

	///	assembles the local right hand side
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

		///	register all functions for on element type
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			typedef FVInnerBoundaryElemDisc T;
			register_prepare_element_loop_function(	id, &T::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &T::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &T::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &T::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &T::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &T::template assemble_A<TElem>);
			register_assemble_M_function(			id, &T::template assemble_M<TElem>);
			register_assemble_f_function(			id, &T::template assemble_f<TElem>);
		}

};

}

#include "inner_boundary_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__*/
