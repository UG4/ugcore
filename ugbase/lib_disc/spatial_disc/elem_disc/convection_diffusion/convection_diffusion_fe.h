/*
 * convection_diffusion_fe.cpp
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion.h"

#include "lib_disc/spatial_disc/disc_util/finite_element_geometry.h"
#include "common/util/provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"
#include "lib_disc/quadrature/gauss_quad/gauss_quad.h"

namespace ug{

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
elem_loop_prepare_fe()
{
//	get reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	static const ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

//	request geometry
	static typename TGeomProvider::Type& geo = TGeomProvider::get();

//	prepare geometry for type and order
	if(!geo.update_local(roid, m_lfeID, m_quadOrder))
		UG_THROW_FATAL("ConvectionDiffusion::elem_loop_prepare_fe:"
						" Cannot update Finite Element Geometry.");

//	set local positions
	m_imDiffusion.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imVelocity.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imSource.template    set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imReactionRate.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imReaction.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imMassScale.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imMass.template 	   set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
elem_loop_finish_fe()
{}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
elem_prepare_fe(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

//	request geometry
	static typename TGeomProvider::Type& geo = TGeomProvider::get();

	if(!geo.update(elem, &m_vCornerCoords[0]))
		UG_THROW_FATAL("ConvectionDiffusion::elem_prepare_fe:"
						" Cannot update Finite Element Geometry.");

//	set global positions for rhs
	m_imDiffusion.	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVelocity. 	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imSource.   	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReactionRate.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReaction. 	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMassScale.	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMass.	  	set_global_ips(geo.global_ips(), geo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ass_JA_elem_fe(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
				if(m_imReactionRate.data_given())
					integrand += m_imReactionRate[ip] * geo.shape(ip, j) * geo.shape(ip, i);

			//	no explicit dependency on m_imReaction

			//	multiply by weight
				integrand *= geo.weight(ip);

			//	add to local matrix
				J(_C_, i, _C_, j) += integrand;
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ass_JM_elem_fe(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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

			//	no explicit dependency on m_imMass

			//	add to local matrix
				J(_C_, i, _C_, j) += val;
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ass_dA_elem_fe(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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

		// 	add Reaction Rate
			if(m_imReactionRate.data_given())
				integrand += m_imReactionRate[ip] * shape_u * geo.shape(ip, i);

		// 	add Reaction
			if(m_imReaction.data_given())
				integrand += m_imReaction[ip] * geo.shape(ip, i);

		//	multiply by integration weight
			integrand *= geo.weight(ip);

		//	add to local defect
			d(_C_, i) += integrand;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ass_dM_elem_fe(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
			number val = shape_u;

		//	add MassScaling
			if(m_imMassScale.data_given())
				val *= m_imMassScale[ip];

		//	add Maxx
			if(m_imMass.data_given())
				val += m_imMass[ip];

		//	add to local defect
			d(_C_, i) +=  val * geo.shape(ip, i) * geo.weight(ip);
		}
	}
};

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ass_rhs_elem_fe(LocalVector& d)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	skip if no source present
	if(!m_imSource.data_given()) return;

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
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_velocity_fe(const LocalVector& u,
                     std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_diffusion_fe(const LocalVector& u,
                      std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                      const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_reaction_fe(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_reaction_rate_fe(const LocalVector& u,
                         std::vector<std::vector<number> > vvvLinDef[],
                         const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
}


//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_source_fe(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_mass_scale_fe(const LocalVector& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
lin_def_mass_fe(const LocalVector& u,
                std::vector<std::vector<number> > vvvLinDef[],
                const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ex_value_fe(const LocalVector& u,
            const MathVector<dim> vGlobIP[],
            const MathVector<TGeomProvider::Type::dim> vLocIP[],
            const size_t nip,
            number vValue[],
            bool bDeriv,
            std::vector<std::vector<number> > vvvDeriv[])
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
			UG_THROW_FATAL("ConvectionDiffusion::ex_value_fe: "<< ex.get_msg());
		}
	}
}

template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::
ex_grad_fe(const LocalVector& u,
           const MathVector<dim> vGlobIP[],
           const MathVector<TGeomProvider::Type::dim> vLocIP[],
           const size_t nip,
           MathVector<dim> vValue[],
           bool bDeriv,
           std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

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
			UG_LOG("ERROR in ConvectionDiffusion::ex_grad_fe: "<< ex.get_msg());
		}
	}
};

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for all dim
template<>
void ConvectionDiffusion<Domain1d>::
register_all_fe_funcs(int order, int quadOrder)
{
//	Edge
	register_fe_func<Edge, FlexGeomProvider<DimFEGeometry<dim, 1> > >();
}

// register for all dim
template<>
void ConvectionDiffusion<Domain2d>::
register_all_fe_funcs(int order, int quadOrder)
{
	if(quadOrder != 2*order+1)
	{
		register_fe_func<Triangle, FlexGeomProvider<DimFEGeometry<dim, 2> > >();
		register_fe_func<Quadrilateral, FlexGeomProvider<DimFEGeometry<dim, 2> > >();
	}

//	special compiled cases

//	Triangle
	switch(order)
	{
		case 1:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1>, GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
				 register_fe_func<Triangle, Provider<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2>, GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
				 register_fe_func<Triangle, Provider<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3>, GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
				 register_fe_func<Triangle, Provider<FEGeom> >(); break;}
		default: register_fe_func<Triangle, FlexGeomProvider<DimFEGeometry<dim, 2> > >();  break;
	}

//	Quadrilateral
	switch(order) {
		case 1:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 1>, GaussQuadrature<ReferenceQuadrilateral, 3> > FEGeom;
				 register_fe_func<Quadrilateral, Provider<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 2>, GaussQuadrature<ReferenceQuadrilateral, 7> > FEGeom;
				 register_fe_func<Quadrilateral, Provider<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 3>, GaussQuadrature<ReferenceQuadrilateral, 11> > FEGeom;
				 register_fe_func<Quadrilateral, Provider<FEGeom> >(); break;}
		default: register_fe_func<Quadrilateral, FlexGeomProvider<DimFEGeometry<dim, 2> > >();  break;
	}
}

// register for all dim
template<>
void ConvectionDiffusion<Domain3d>::
register_all_fe_funcs(int order, int quadOrder)
{
	if(quadOrder != 2*order+1)
	{
		register_fe_func<Tetrahedron, FlexGeomProvider<DimFEGeometry<dim, 3> > >();
		register_fe_func<Prism, FlexGeomProvider<DimFEGeometry<dim, 3> > >();
		register_fe_func<Pyramid, FlexGeomProvider<DimFEGeometry<dim, 3> > >();
		register_fe_func<Hexahedron, FlexGeomProvider<DimFEGeometry<dim, 3> > >();
	}

//	special compiled cases

//	Tetrahedron
	switch(order)
	{
		case 1:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 1>, GaussQuadrature<ReferenceTetrahedron, 3> > FEGeom;
				 register_fe_func<Tetrahedron, Provider<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 2>, GaussQuadrature<ReferenceTetrahedron, 5> > FEGeom;
				 register_fe_func<Tetrahedron, Provider<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 3>, GaussQuadrature<ReferenceTetrahedron, 7> > FEGeom;
				 register_fe_func<Tetrahedron, Provider<FEGeom> >(); break;}
		default: register_fe_func<Tetrahedron, FlexGeomProvider<DimFEGeometry<dim, 3> > >();  break;
	}

//	Prism
	switch(order) {
		case 1:	{typedef FEGeometry<Prism, dim, LagrangeLSFS<ReferencePrism, 1>, GaussQuadrature<ReferencePrism, 2> > FEGeom;
				 register_fe_func<Prism, Provider<FEGeom> >(); break;}
		default: register_fe_func<Prism, FlexGeomProvider<DimFEGeometry<dim, 3> > >();  break;
	}

//	Pyramid
	switch(order)
	{
		default: register_fe_func<Pyramid, FlexGeomProvider<DimFEGeometry<dim, 3> > >();  break;
	}

//	Hexahedron
	switch(order)
	{
		case 1:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 1>, GaussQuadrature<ReferenceHexahedron, 3> > FEGeom;
				 register_fe_func<Hexahedron, Provider<FEGeom> >(); break;}
		case 2:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 2>, GaussQuadrature<ReferenceHexahedron, 7> > FEGeom;
				 register_fe_func<Hexahedron, Provider<FEGeom> >(); break;}
		case 3:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 3>, GaussQuadrature<ReferenceHexahedron, 11> > FEGeom;
				 register_fe_func<Hexahedron, Provider<FEGeom> >(); break;}
		default: register_fe_func<Hexahedron, FlexGeomProvider<DimFEGeometry<dim, 3> > >();  break;
	}
}


template <typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusion<TDomain>::register_fe_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->enable_fast_ass_elem(true);
	this->set_prep_elem_loop_fct(id, &T::template elem_loop_prepare_fe<TElem, TGeomProvider>);
	this->set_prep_elem_fct(	   id, &T::template elem_prepare_fe<TElem, TGeomProvider>);
	this->set_fsh_elem_loop_fct( id, &T::template elem_loop_finish_fe<TElem, TGeomProvider>);
	this->set_ass_JA_elem_fct(   id, &T::template ass_JA_elem_fe<TElem, TGeomProvider>);
	this->set_ass_JM_elem_fct(   id, &T::template ass_JM_elem_fe<TElem, TGeomProvider>);
	this->set_ass_dA_elem_fct(   id, &T::template ass_dA_elem_fe<TElem, TGeomProvider>);
	this->set_ass_dM_elem_fct(   id, &T::template ass_dM_elem_fe<TElem, TGeomProvider>);
	this->set_ass_rhs_elem_fct(  id, &T::template ass_rhs_elem_fe<TElem, TGeomProvider>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity_fe<TElem, TGeomProvider>);
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion_fe<TElem, TGeomProvider>);
	m_imReactionRate. set_fct(id, this, &T::template lin_def_reaction_rate_fe<TElem, TGeomProvider>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction_fe<TElem, TGeomProvider>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source_fe<TElem, TGeomProvider>);
	m_imMassScale.set_fct(id, this, &T::template lin_def_mass_scale_fe<TElem, TGeomProvider>);
	m_imMass.	  set_fct(id, this, &T::template lin_def_mass_fe<TElem, TGeomProvider>);

//	exports
	m_exConcentration->	  template set_fct<T,refDim>(id, this, &T::template ex_value_fe<TElem, TGeomProvider>);
	m_exConcentrationGrad->template set_fct<T,refDim>(id, this, &T::template ex_grad_fe<TElem, TGeomProvider>);
}

} // namespace ug

