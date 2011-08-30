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
				MatVecMult(Dgrad, m_imDiffusion[ip], geo.global_grad(ip, j));
			else
				VecSet(Dgrad, 0.0);

		//  Convection
			if(m_imVelocity.data_given())
				VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_imVelocity[ip]);

		//	loop test space
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	compute integrand
				number integrand = VecDot(Dgrad, geo.global_grad(ip, i));

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
			VecScaleAppend(grad_u, u(_C_,j), geo.global_grad(ip, j));
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

	// 	Diffusion
		if(m_imDiffusion.data_given())
			MatVecMult(Dgrad_u, m_imDiffusion[ip], grad_u);
		else
			VecSet(Dgrad_u, 0.0);

	// 	Convection
		if(m_imVelocity.data_given())
			VecScaleAppend(Dgrad_u, -1*shape_u, m_imVelocity[ip]);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute integrand
			integrand = VecDot(Dgrad_u, geo.global_grad(ip, i));

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

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_def_velocity_fe(const local_vector_type& u,
                     std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add to local defect
			VecScale(vvvLinDef[ip][_C_][i], geo.global_grad(ip, i),
			         	 	 	 	 	 	 geo.weight(ip) * shape_u);
		}
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_def_diffusion_fe(const local_vector_type& u,
                      std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                      const size_t nip)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

	MathVector<dim> grad_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(grad_u, 0.0);
		for(size_t j = 0; j < geo.num_sh(); ++j)
			VecScaleAppend(grad_u, u(_C_,j), geo.global_grad(ip, j));

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t k=0; k < (size_t)dim; ++k)
				for(size_t j = 0; j < (size_t)dim; ++j)
					(vvvLinDef[ip][_C_][i])(k,j) += grad_u[j] * geo.global_grad(ip, i)[k]
												* geo.weight(ip);
		}
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_def_reaction_fe(const local_vector_type& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_def_source_fe(const local_vector_type& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add contribution to local defect
			vvvLinDef[ip][_C_][i] = geo.shape(ip, i) * geo.weight(ip);
		}
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
lin_def_mass_scale_fe(const local_vector_type& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}

//	we're done
	return true;
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
ex_concentration_fe(const local_vector_type& u,
                     const MathVector<dim> vGlobIP[],
                     const MathVector<THolder::Type::dim> vLocIP[],
                     const size_t nip,
                     number vValue[],
                     bool bDeriv,
                     std::vector<std::vector<number> > vvvDeriv[])
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::num_corners;

//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	FE ip
	if(vLocIP == geo.local_ips())
	{
	//	Loop ips
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				vValue[ip] += u(_C_, sh) * geo.shape(ip, sh);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = geo.shape(ip, sh);
		}
	}
// 	general case
	else
	{
	//	request for trial space
		try{
		const DimLocalShapeFunctionSet<dim>& rTrialSpace
			 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u(_C_, sh) * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_C_][sh] = vShape[sh];
		}

		}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
			UG_LOG("ERROR in ConvectionDiffusionElemDisc::ex_concentration_fe: "
					<< ex.get_msg() << ".\n");
			return false;
		}
	}

//	we're done
	return true;
}

template<typename TDomain>
template <typename TElem, typename THolder>
bool ConvectionDiffusionElemDisc<TDomain>::
ex_concentration_grad_fe(const local_vector_type& u,
                          const MathVector<dim> vGlobIP[],
                          const MathVector<THolder::Type::dim> vLocIP[],
                          const size_t nip,
                          MathVector<dim> vValue[],
                          bool bDeriv,
                          std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
//	request geometry
	static const typename THolder::Type& geo = THolder::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::num_corners;

//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	FE
	if(vLocIP == geo.local_ips())
	{
	//	Loop ip
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			VecSet(vValue[ip], 0.0);
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u(_C_, sh), geo.global_grad(ip, sh));

			if(bDeriv)
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = geo.global_grad(ip, sh);
		}
	}
// 	general case
	else
	{
	//	request for trial space
		try{
		const DimLocalShapeFunctionSet<dim>& rTrialSpace
			 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(&m_vCornerCoords[0]);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);
		}

		}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
			UG_LOG("ERROR in ConvectionDiffusionElemDisc::ex_concentration_grad_fe: "
					<< ex.get_msg() << ".\n");
			return false;
		}
	}

//	we're done
	return true;
};

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
		case 1:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1>, GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
				 register_fe_func<Triangle, ClassHolder<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2>, GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
				 register_fe_func<Triangle, ClassHolder<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3>, GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
				 register_fe_func<Triangle, ClassHolder<FEGeom> >(); break;}
		default: register_fe_func<Triangle, FlexGeomHolder<DimFEGeometry<dim, 2> > >();  break;
	}

//	Quadrilateral
	switch(order) {
		case 1:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 1>, GaussQuadrature<ReferenceQuadrilateral, 3> > FEGeom;
				 register_fe_func<Quadrilateral, ClassHolder<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 2>, GaussQuadrature<ReferenceQuadrilateral, 7> > FEGeom;
				 register_fe_func<Quadrilateral, ClassHolder<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 3>, GaussQuadrature<ReferenceQuadrilateral, 11> > FEGeom;
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
	static const int refDim = reference_element_traits<TElem>::dim;

	set_prep_elem_loop_fct(id, &T::template elem_loop_prepare_fe<TElem, THolder>);
	set_prep_elem_fct(	   id, &T::template elem_prepare_fe<TElem, THolder>);
	set_fsh_elem_loop_fct( id, &T::template elem_loop_finish_fe<TElem, THolder>);
	set_ass_JA_elem_fct(   id, &T::template elem_JA_fe<TElem, THolder>);
	set_ass_JM_elem_fct(   id, &T::template elem_JM_fe<TElem, THolder>);
	set_ass_dA_elem_fct(   id, &T::template elem_dA_fe<TElem, THolder>);
	set_ass_dM_elem_fct(   id, &T::template elem_dM_fe<TElem, THolder>);
	set_ass_rhs_elem_fct(  id, &T::template elem_rhs_fe<TElem, THolder>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity_fe<TElem, THolder>);
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion_fe<TElem, THolder>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction_fe<TElem, THolder>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source_fe<TElem, THolder>);
	m_imMassScale.set_fct(id, this, &T::template lin_def_mass_scale_fe<TElem, THolder>);

//	exports
	m_exConcentration.	  template set_fct<T,refDim>(id, this, &T::template ex_concentration_fe<TElem, THolder>);
	m_exConcentrationGrad.template set_fct<T,refDim>(id, this, &T::template ex_concentration_grad_fe<TElem, THolder>);
}

} // namespace ug

