/*
 * constant_equation_impl.h
 *
 *  Created on: 13.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__TIME_NEUMANN_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__TIME_NEUMANN_BOUNDARY__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/ip_data/data_import_export.h"

namespace ug{

/// \ingroup lib_disc_elem_disc
/// @{

template<	typename TDomain>
class FV1TimeNeumannBoundary : public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	Own type
		typedef FV1TimeNeumannBoundary<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	///	Local matrix type
		typedef typename base_type::local_matrix_type local_matrix_type;

	///	Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Local index type
		typedef typename base_type::local_index_type local_index_type;

	public:
	///	Constructor
		FV1TimeNeumannBoundary();

	///	sets the capacity
		void set_capacity(number cap)	{m_capacity = cap;}

	public:
	///	number of functions used
		virtual size_t num_fct(){return 5;}

	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
		{
		//	check number
			if(vLfeID.size() != num_fct()) return false;

		//	check that Lagrange 1st order
			return vLfeID[0] == LFEID(LFEID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			register_all_fv1_funcs(bNonRegular);

		//	this disc supports both grids
			return true;
		}

		virtual bool use_hanging() const {return true;}

	private:
	///	prepares the discretization for time dependent discretization
	/**
	 * This function prepares the discretization for time-dependent problems.
	 * It sets the time in the imports.
	 *
	 * \param[in]	time	new time point
	 * \returns 	true	indicates, that old values are needed
	 */
		virtual bool time_point_changed(number time);

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
		bool prepare_element(TElem* elem, const local_vector_type& u);

	///	finishes the loop over all elements
		template <typename TElem>
		inline bool finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem>
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem>
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u);

	///	assembles the stiffness part of the local defect
		template <typename TElem>
		bool assemble_A(local_vector_type& d, const local_vector_type& u);

	///	assembles the mass part of the local defect
		template <typename TElem>
		bool assemble_M(local_vector_type& d, const local_vector_type& u);

	///	assembles the local right hand side
		template <typename TElem>
		bool assemble_f(local_vector_type& d);

	private:
	///	Corner Coordinates
		const position_type* m_vCornerCoords;

	///	abbreviation for the local solution
		static const size_t _PHI_ = 0;
		static const size_t _VM_ = 1;
		static const size_t _N_ = 2;
		static const size_t _M_ = 3;
		static const size_t _H_ = 4;

	///	Data import for the right-hand side
		number m_capacity;

	private:
		void register_all_fv1_funcs(bool bHang);

		template <typename TElem>
		void register_fv1_func();

};

/// @}

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__TIME_NEUMANN_BOUNDARY__*/
