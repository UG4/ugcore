/*
 * convection_diffusion_fe.cpp
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion.h"

#include "lib_discretization/spatial_discretization/disc_util/finite_element_geometry.h"
#include "common/util/provider.h"
#include "lib_discretization/local_finite_element/lagrange/lagrange.h"
#include "lib_discretization/local_finite_element/lagrange/lagrangep1.h"
#include "lib_discretization/quadrature/gauss_quad/gauss_quad.h"

namespace ug{

template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_loop_prepare_fe()
{
//	get reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	static const ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

//	request geometry
	static typename THolder::Type& geo = THolder::get();

//	prepare geometry for type and order
	if(!geo.update_local(roid, m_lfeID, m_quadOrder))
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::elem_loop_prepare_fe':"
				" Cannot update Finite Element Geometry.\n");
		return false;
	}

//	set local positions
	m_imDiffusion.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip());
	m_imVelocity.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip());
	m_imSource.template    set_local_ips<refDim>(geo.local_ips(), geo.num_ip());
	m_imReaction.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip());
	m_imMassScale.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip());

//	done
	return true;
}

template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_loop_finish_fe()
{
//	nothing to do
	return true;
}

template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_prepare_fe(TElem* elem, const local_vector_type& u)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

//	request geometry
	static typename THolder::Type& geo = THolder::get();

	if(!geo.update(elem, &m_vCornerCoords[0]))
	{
		UG_LOG("ERROR in 'ConvectionDiffusionElemDisc::prepare_element':"
				" Cannot update Finite Element Geometry.\n");
		return false;
	}

//	set global positions for rhs
	m_imDiffusion.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVelocity. set_global_ips(geo.global_ips(), geo.num_ip());
	m_imSource.   set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReaction. set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMassScale.set_global_ips(geo.global_ips(), geo.num_ip());

//	done
	return true;
}

template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_JA_fe(local_matrix_type& J, const local_vector_type& u)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

	MathVector<dim> v, Dgrad;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop trial space
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
		//	Diffusion
			if(m_imDiffusion.data_given())
				MatVecMult(Dgrad, m_imDiffusion[ip], geo.grad_global(ip, j));
			else
				VecSet(Dgrad, 0.0);

		//  Convection
			if(m_imVelocity.data_given())
				VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_imVelocity[ip]);

		//	loop test space
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	compute integrand
				number integrand = VecDot(Dgrad, geo.grad_global(ip, i));

			// 	Reaction
				if(m_imReaction.data_given())
					integrand += m_imReaction[ip] * geo.shape(ip, j) * geo.shape(ip, i);

			//	multiply by weight
				integrand *= geo.weight(ip);

			//	add to local matrix
				J(_C_, i, _C_, j) += integrand;
			}
		}
	}

//	done
	return true;
}


template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_JM_fe(local_matrix_type& J, const local_vector_type& u)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test space
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	loop trial space
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
			//	compute integrand
				number val = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);

			//	add MassScale
				if(m_imMassScale.data_given())
					val *= m_imMassScale[ip];

			//	add to local matrix
				J(_C_, i, _C_, j) += val;
			}
		}
	}

//	done
	return true;
}


template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_dA_fe(local_vector_type& d, const local_vector_type& u)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

	number integrand, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(grad_u, 0.0);
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			VecScaleAppend(grad_u, u(_C_,j), geo.grad_global(ip, j));
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

	// 	Diffusion
		if(m_imDiffusion.data_given())
			MatVecMult(Dgrad_u, m_imDiffusion[ip], grad_u);
		else
			VecSet(Dgrad_u, 0.0);

	// 	Convection
		if(m_imReaction.data_given())
			VecScaleAppend(Dgrad_u, -1*shape_u, m_imVelocity[ip]);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute integrand
			integrand = VecDot(Dgrad_u, geo.grad_global(ip, i));

		// 	add Reaction
			if(m_imReaction.data_given())
				integrand += m_imReaction[ip] * shape_u * geo.shape(ip, i);

		//	multiply by integration weight
			integrand *= geo.weight(ip);

		//	add to local defect
			d(_C_, i) += integrand;
		}
	}

	return true;
}


template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_dM_fe(local_vector_type& d, const local_vector_type& u)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

	number shape_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add MassScaling
			if(m_imMassScale.data_given())
				val *= m_imMassScale[ip];

		//	add to local defect
			d(_C_, i) +=  val;
		}
	}

//	done
	return true;
};

template<typename TDomain>
template<typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
elem_rhs_fe(local_vector_type& d)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	skip if no source present
	if(!m_imSource.data_given()) return true;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add contribution to local defect
			d(_C_, i) +=  m_imSource[ip] * geo.shape(ip, i) * geo.weight(ip);
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for all dim
template<>
void ConvectionDiffusionElemDisc<Domain1d>::
register_all_fe_funcs(int order, int quadOrder)
{
//	Edge
	register_fe_func<Edge, FlexGeomHolder<DimFEGeometry<dim, 1> > >();
}

// register for all dim
template<>
void ConvectionDiffusionElemDisc<Domain2d>::
register_all_fe_funcs(int order, int quadOrder)
{
	if(quadOrder != 2*order+1)
	{
		register_fe_func<Triangle, FlexGeomHolder<DimFEGeometry<dim, 2> > >();
		register_fe_func<Quadrilateral, FlexGeomHolder<DimFEGeometry<dim, 2> > >();
	}

//	special compiled cases

//	Triangle
	switch(order)
	{
		case 1:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1>, GaussQuadrature<ReferenceTriangle, 2> > FEGeom;
				 register_fe_func<Triangle, ClassHolder<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2>, GaussQuadrature<ReferenceTriangle, 4> > FEGeom;
				 register_fe_func<Triangle, ClassHolder<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3>, GaussQuadrature<ReferenceTriangle, 6> > FEGeom;
				 register_fe_func<Triangle, ClassHolder<FEGeom> >(); break;}
		default: register_fe_func<Triangle, FlexGeomHolder<DimFEGeometry<dim, 2> > >();  break;
	}

//	Quadrilateral
	switch(order) {
		case 1:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 1>, GaussQuadrature<ReferenceQuadrilateral, 2> > FEGeom;
				 register_fe_func<Quadrilateral, ClassHolder<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 2>, GaussQuadrature<ReferenceQuadrilateral, 6> > FEGeom;
				 register_fe_func<Quadrilateral, ClassHolder<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 3>, GaussQuadrature<ReferenceQuadrilateral, 10> > FEGeom;
				 register_fe_func<Quadrilateral, ClassHolder<FEGeom> >(); break;}
		default: register_fe_func<Quadrilateral, FlexGeomHolder<DimFEGeometry<dim, 2> > >();  break;
	}
}

// register for all dim
template<>
void ConvectionDiffusionElemDisc<Domain3d>::
register_all_fe_funcs(int order, int quadOrder)
{
	if(quadOrder != 2*order+1)
	{
		register_fe_func<Tetrahedron, FlexGeomHolder<DimFEGeometry<dim, 3> > >();
		register_fe_func<Prism, FlexGeomHolder<DimFEGeometry<dim, 3> > >();
		register_fe_func<Pyramid, FlexGeomHolder<DimFEGeometry<dim, 3> > >();
		register_fe_func<Hexahedron, FlexGeomHolder<DimFEGeometry<dim, 3> > >();
	}

//	special compiled cases

//	Tetrahedron
	switch(order)
	{
		case 1:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 1>, GaussQuadrature<ReferenceTetrahedron, 2> > FEGeom;
				 register_fe_func<Tetrahedron, ClassHolder<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 2>, GaussQuadrature<ReferenceTetrahedron, 4> > FEGeom;
				 register_fe_func<Tetrahedron, ClassHolder<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 3>, GaussQuadrature<ReferenceTetrahedron, 6> > FEGeom;
				 register_fe_func<Tetrahedron, ClassHolder<FEGeom> >(); break;}
		default: register_fe_func<Tetrahedron, FlexGeomHolder<DimFEGeometry<dim, 3> > >();  break;
	}

//	Prism
	switch(order) {
		case 1:	{typedef FEGeometry<Prism, dim, LagrangeLSFS<ReferencePrism, 1>, GaussQuadrature<ReferencePrism, 2> > FEGeom;
				 register_fe_func<Prism, ClassHolder<FEGeom> >(); break;}
		default: register_fe_func<Prism, FlexGeomHolder<DimFEGeometry<dim, 3> > >();  break;
	}

//	Pyramid
	switch(order)
	{
		default: register_fe_func<Pyramid, FlexGeomHolder<DimFEGeometry<dim, 3> > >();  break;
	}

//	Hexahedron
	switch(order)
	{
		case 1:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 1>, GaussQuadrature<ReferenceHexahedron, 2> > FEGeom;
				 register_fe_func<Hexahedron, ClassHolder<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 2>, GaussQuadrature<ReferenceHexahedron, 6> > FEGeom;
				 register_fe_func<Hexahedron, ClassHolder<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 3>, GaussQuadrature<ReferenceHexahedron, 10> > FEGeom;
				 register_fe_func<Hexahedron, ClassHolder<FEGeom> >(); break;}
		default: register_fe_func<Hexahedron, FlexGeomHolder<DimFEGeometry<dim, 3> > >();  break;
	}
}


template<typename TDomain>
template<typename TElem, typename THolder>
void ConvectionDiffusionElemDisc<TDomain>::
register_fe_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	reg_prepare_elem_loop_fct(id, &T::template elem_loop_prepare_fe<TElem, THolder>);
	reg_prepare_elem_fct(	  id, &T::template elem_prepare_fe<TElem, THolder>);
	reg_finish_elem_loop_fct( id, &T::template elem_loop_finish_fe<TElem, THolder>);
	reg_ass_JA_elem_fct(	  id, &T::template elem_JA_fe<TElem, THolder>);
	reg_ass_JM_elem_fct(	  id, &T::template elem_JM_fe<TElem, THolder>);
	reg_ass_dA_elem_fct(	  id, &T::template elem_dA_fe<TElem, THolder>);
	reg_ass_dM_elem_fct(	  id, &T::template elem_dM_fe<TElem, THolder>);
	reg_ass_rhs_elem_fct(	  id, &T::template elem_rhs_fe<TElem, THolder>);
}

} // namespace ug

