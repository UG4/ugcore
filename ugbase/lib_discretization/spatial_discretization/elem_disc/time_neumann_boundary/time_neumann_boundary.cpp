/*
 * constant_equation_fv1.cpp
 *
 *  Created on: 13.05.2011
 *      Author: andreasvogel
 */

#include "time_neumann_boundary.h"

#include "common/util/provider.h"
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"
#include "lib_discretization/spatial_discretization/disc_util/hanging_finite_volume_geometry.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//  Provide a generic implementation for all elements
//  (since this discretization can be implemented in a generic way)
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
prepare_element_loop()
{
//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
finish_element_loop()
{
//	nothing to do
	return true;
}

template<typename TDomain>
bool FV1TimeNeumannBoundary<TDomain>::
time_point_changed(number time)
{
//	this disc does needs the old time solutions, thus, return true
	return true;
}


template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
prepare_element(TElem* elem, const local_vector_type& u)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

// 	Update Geometry for this element
	FV1ManifoldBoundary<TElem, dim>& geo = Provider<FV1ManifoldBoundary<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->get_subset_handler())))
	{
		UG_LOG("FV1TimeNeumannBoundary::prepare_element:"
				" Cannot update Finite Volume Geometry.\n"); return false;
	}

//	we're done
	return true;
}

template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
// 	get finite volume geometry
	const static FV1ManifoldBoundary<TElem, dim>& geo = Provider<FV1ManifoldBoundary<TElem,dim> >::get();

//	check that time dependent assembling, else meaningless disc
	if(!this->is_time_dependent())
	{
		UG_LOG("ERROR in 'FV1TimeNeumannBoundary::assemble_A': "
				" Only time-dependent discs supported.\n");
		return false;
	}

//	get and check current and old solution
	const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();

//	check
	if(vLocSol.size() != 2)
	{
		UG_LOG("ERROR in 'FV1TimeNeumannBoundary::assemble_A': "
				" Disc needs exactly two time points.\n");
		return false;
	}

	const local_vector_type& sol = vLocSol.solution(0);


//	remember local solutions
	number dt = vLocSol.time(0) - vLocSol.time(1);

// 	loop BF (i.e. the bf of the higher dim object)
	for(size_t ip = 0; ip < geo.num_bf(); ++ip)
	{
	// 	get current SCVF
		const typename FV1ManifoldBoundary<TElem, dim>::BF& bf = geo.bf(ip);

	//	get node id
		const size_t id = bf.node_id();


	//  add to local jacobian
		J(_PHI_, id, _VM_, id) += bf.volume() * (m_capacity / dt +36*pow(sol(_N_,id),4)
								+ 120*pow(sol(_M_,id),3)*sol(_H_,id) + 0.3*sol(_VM_,id) );
	}

//	done
	return true;
}


template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
//	no own contribution to jacobian (only due to linearized defect)
	return true;
}


template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
// 	get finite volume geometry
	const static FV1ManifoldBoundary<TElem, dim>& geo = Provider<FV1ManifoldBoundary<TElem,dim> >::get();

//	check that time dependent assembling, else meaningless disc
	if(!this->is_time_dependent())
	{
		UG_LOG("ERROR in 'FV1TimeNeumannBoundary::assemble_A': "
				" Only time-dependent discs supported.\n");
		return false;
	}

//	get and check current and old solution
	const LocalVectorTimeSeries& vLocSol = *this->local_time_solutions();

//	check
	if(vLocSol.size() != 2)
	{
		UG_LOG("ERROR in 'FV1TimeNeumannBoundary::assemble_A': "
				" Disc needs exactly two time points.\n");
		return false;
	}

//	remember local solutions
	const local_vector_type& sol = vLocSol.solution(0);
	const local_vector_type& oldSol = vLocSol.solution(1);
	number dt = vLocSol.time(0) - vLocSol.time(1);

// 	loop BF (i.e. the bf of the higher dim object)
	for(size_t ip = 0; ip < geo.num_bf(); ++ip)
	{
	// 	get current SCVF
		const typename FV1ManifoldBoundary<TElem, dim>::BF& bf = geo.bf(ip);

		UG_LOG("SIZE Boundary face on manifold: "<<bf.volume()<<"\n");

	//	get node id
		const size_t id = bf.node_id();

	//  add to local defect
		d(_PHI_, id) += bf.volume() * (m_capacity *
						( sol(_VM_, id) - oldSol(_VM_, id) ) / dt + 36*pow(sol(_N_,id),4)*(sol(_VM_,id) + 77)
						+ 120*pow(sol(_M_,id),3)*sol(_H_,id) * (sol(_VM_, id) - 50) + 0.3*(sol(_VM_,id) - 54.4) );
	}

// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
// 	we're done
	return true;
}


template<typename TDomain>
template<typename TElem>
bool FV1TimeNeumannBoundary<TDomain>::assemble_f(local_vector_type& d)
{
// 	we're done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
FV1TimeNeumannBoundary<TDomain>::
FV1TimeNeumannBoundary()
{
//	register assemling functions
	register_all_fv1_funcs(false);
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

template<>
void FV1TimeNeumannBoundary<Domain1d>::
register_all_fv1_funcs(bool bHang)
{
	throw(UGFatalError("Not implemented."));
}

template<>
void FV1TimeNeumannBoundary<Domain2d>::
register_all_fv1_funcs(bool bHang)
{
	register_fv1_func<Edge>();
}

template<>
void FV1TimeNeumannBoundary<Domain3d>::
register_all_fv1_funcs(bool bHang)
{
	register_fv1_func<Triangle>();
	register_fv1_func<Quadrilateral>();
}


template<typename TDomain>
template<typename TElem>
void FV1TimeNeumannBoundary<TDomain>::register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	set_prep_elem_loop_fct(	id, &T::template prepare_element_loop<TElem>);
	set_prep_elem_fct(		id, &T::template prepare_element<TElem>);
	set_fsh_elem_loop_fct( 	id, &T::template finish_element_loop<TElem>);
	set_ass_JA_elem_fct(	id, &T::template assemble_JA<TElem>);
	set_ass_JM_elem_fct(	id, &T::template assemble_JM<TElem>);
	set_ass_dA_elem_fct(	id, &T::template assemble_A<TElem>);
	set_ass_dM_elem_fct(	id, &T::template assemble_M<TElem>);
	set_ass_rhs_elem_fct(	id, &T::template assemble_f<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FV1TimeNeumannBoundary<Domain1d>;
template class FV1TimeNeumannBoundary<Domain2d>;
template class FV1TimeNeumannBoundary<Domain3d>;


} // namespace ug

