/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__

#include "common/common.h"
#include "common/util/provider.h"
#include "common/math/ugmath_types.h"
#include "lib_disc/spatial_disc/disc_util/fv_geom_base.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// Interface for Convection shapes
/////////////////////////////////////////////////////////////////////////////

/// Interface class for upwind methods
/**
 * Upwind is a way of stabilization of discretizations of convection(-diffusion)
 * PDEs. For a detailed description of the idea, cf.
 * <ul>
 *  <li> P. Frolkovic, H. De Schepper, Numerical Modelling of Convection Dominated
 *       Transport Coupled with Density Driven Flow in Porous Media.
 *       Advances in Water Resources 24 (1): 63–72, 2000, DOI: 10.1016/S0309-1708(00)00025-7 </li>
 *  <li> K. W. Morton, D. Mayers, Numerical Solution of Partial Differential Equations:
 *       An Introduction (2nd ed.), Cambridge: Cambridge University Press, 2005, DOI: 10.1017/CBO9780511812248 </li>
 *  <li> R. J. LeVeque, Numerical Methods for Conservation Laws (2nd ed.), Springer Basel AG, 1992, DOI: 10.1007/978-3-0348-8629-1 </li>
 * </ul>
 * Convection shapes is one of possible implementations of upwind methods,
 * in which the upwind is considered as a special kind of interpolation of degrees
 * of freedom of the unknown function (for ex., concentration at the corners of
 * the element) into integration points inside the grid element. E.g. the so-called
 * 'full upwind' method for the vertex-centered FV discretizations sets the
 * concentration in the whole grid element be equal to the concentration at a
 * certain specially selected 'upwind corner' chosen according to the direction
 * of the convection velocity. Like every interpolation, this one can be introduced
 * by a set of shape functions associated with the degrees of freedom of the
 * grid function. These shape functions for an upwind method are said to be
 * its convection shapes.
 *
 * But in terms of this interface, the convection shapes are the values of the
 * convection shape functions at the integration points premultiplied by the norm
 * of the projection of the velocity to the normal to the corresponding
 * subcontrol volume faces and by the area of these faces. Thus, for a convective
 * term of the form
 * \f[ \nabla \cdot (\mathbf{v} c) \f]
 * (with \f$\mathbf{v}\f$ being the convective velocity) the convective flux
 * through a subcontrol volume face with the integration point \f$ip\f$ is
 * \f[ \sum_{s=1}^{n_s} \phi_s (ip) \, c_s, \f]
 * where \f$\phi_s (ip)\f$ are the values of the \f$n_s\f$ ``convection shapes''
 * at \f$ip\f$, whereas \f$c_s\f$ are the degrees of freedom of \f$c\f$ corresponding
 * to the shapes.
 *
 * Although the convection shapes for every upwind method are based on the
 * same idea, they depend on the geometry of the secondary grid (i.e., on the
 * finite volume geometry). For this, for every FV geometry, an own object
 * of the convection shapes class is registered.
 *
 * Note furthermore that the convection shapes may depend not only on the
 * velocity but also on the diffusion. Such methods raise up the asymptotic
 * convergence of the discretization of convection-diffusion problems. As
 * the diffusion tensor may also depend on the unknowns (via its dispersion
 * part), the derivatives of the shapes w.r.t. the components of the diffusion
 * tensor are taken into account.
 *
 * \tparam dim	the world dimension (in which the velocity is given)
 */
template <int dim>
class IConvectionShapes
{
	public:
	/// Abbreviation for own type
		typedef IConvectionShapes<dim> this_type;

	///	type of update function
		typedef bool (this_type::*UpdateFunc)
				(const FVGeometryBase* geo,
				 const MathVector<dim>* Velocity,
				 const MathMatrix<dim, dim>* Diffusion,
				 bool computeDeriv);

	public:
	///	constructor
		IConvectionShapes()
			: m_numScvf(0), m_numSh(0), m_bNonZeroDerivDiffusion(true)
		{
			m_vUpdateFunc.clear();
			m_vUpShape.clear();
			m_vUpShapeVel.clear();
			m_vUpShapeDiffusion.clear();
		}

	/// destructor
		virtual ~IConvectionShapes() {};

	///	returns number of shapes
		size_t num_sh() const {return m_numSh;}

	///	returns number of sub control volume faces
		size_t num_scvf() const {return m_numScvf;}

	/// shape value
		number operator()(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShape.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShape[scvf].size(), "Invalid index");
			return m_vUpShape[scvf][sh];
		}

	///	upwind shape for corner vel
		const MathVector<dim>& D_vel(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShapeVel.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeVel[scvf].size(), "Invalid index");
			return m_vUpShapeVel[scvf][sh];
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_deriv_diffusion() const {return m_bNonZeroDerivDiffusion;}

	///	upwind shapes for ip vel
		const MathMatrix<dim,dim>& D_diffusion(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShapeDiffusion.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeDiffusion[scvf].size(), "Invalid index");
			return m_vUpShapeDiffusion[scvf][sh];
		}

		bool update(const FVGeometryBase* geo,
					const MathVector<dim>* Velocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv)
			{return (this->*(m_vUpdateFunc[m_id])) (geo, Velocity, DiffDisp, computeDeriv);}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	resize the data arrays
		void set_sizes(size_t numScvf, size_t numSh);

	///	sets the shape ip flag
		void set_non_zero_deriv_diffusion_flag(bool flag) {m_bNonZeroDerivDiffusion = flag;}

	/// non-const access to ip velocity (i.e. interpolated velocity at ip)
		number& conv_shape(size_t scvf, size_t sh)
		{
			UG_ASSERT(scvf < m_vUpShape.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShape[scvf].size(), "Invalid index");
			return m_vUpShape[scvf][sh];
		}

	///	non-const access to upwind shapes for corner vel
		MathVector<dim>& D_vel(size_t scvf, size_t sh)
		{
			UG_ASSERT(scvf < m_vUpShapeVel.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeVel[scvf].size(), "Invalid index");
			return m_vUpShapeVel[scvf][sh];
		}

	///	non-const access to upwind shapes for ip vel
		MathMatrix<dim,dim>& conv_shape_diffusion(size_t scvf, size_t sh)
		{
			UG_ASSERT(scvf < m_vUpShapeDiffusion.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeDiffusion[scvf].size(), "Invalid index");
			return m_vUpShapeDiffusion[scvf][sh];
		}


	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	upwind shapes for corners shape functions
		std::vector<std::vector<number> > m_vUpShape;

	///	upwind shapes for corners shape functions
		std::vector<std::vector<MathVector<dim> > > m_vUpShapeVel;

	///	upwind shapes for corners shape functions
		std::vector<std::vector<MathMatrix<dim,dim> > > m_vUpShapeDiffusion;

	///	flag if ip shapes are non-zero
		bool m_bNonZeroDerivDiffusion;

	//////////////////////////
	// registering process
	//////////////////////////

	public:
	///	register a update function for a Geometry
		template <typename TFVGeom, typename TAssFunc>
		void register_update_func(TAssFunc func);

	///	set the Geometry type to use for next updates
		template <typename TFVGeom>
		bool set_geometry_type(const TFVGeom& geo);

	protected:
	///	Vector holding all update functions
		std::vector<UpdateFunc> m_vUpdateFunc;

	///	id of current geometry type
		int m_id;
};

/////////////////////////////////////////////////////////////////////////////
// Interface for Stabilization
/////////////////////////////////////////////////////////////////////////////

//	register a update function for a Geometry
template <int dim>
template <typename TFVGeom, typename TAssFunc>
void
IConvectionShapes<dim>::
register_update_func(TAssFunc func)
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	make sure that there is enough space
	if((size_t)id >= m_vUpdateFunc.size())
		m_vUpdateFunc.resize(id+1, NULL);

//	set pointer
	m_vUpdateFunc[id] = (UpdateFunc)func;
}

//	set the Geometry type to use for next updates
template <int dim>
template <typename TFVGeom>
bool
IConvectionShapes<dim>::
set_geometry_type(const TFVGeom& geo)
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	check that function exists
	if(id >= m_vUpdateFunc.size() || m_vUpdateFunc[id] == NULL)
	{
		UG_LOG("ERROR in 'IConvectionShapes::set_geometry_type':"
				" No update function registered for this Geometry.\n");
		return false;
	}

//	set current geometry
	m_id = id;

//	set sizes
	set_sizes(geo.num_scvf(), geo.num_sh());

//	we're done
	return true;
}


//	resize the data arrays
template <int dim>
void
IConvectionShapes<dim>::
set_sizes(size_t numScvf, size_t numSh)
{
//	remember sizes
	m_numScvf = numScvf;
	m_numSh = numSh;

//	adjust arrays
	m_vUpShape.resize(m_numScvf);
	m_vUpShapeVel.resize(m_numScvf);
	m_vUpShapeDiffusion.resize(m_numScvf);
	for(size_t scvf = 0; scvf < m_numScvf; ++scvf)
	{
		m_vUpShape[scvf].resize(m_numSh);
		m_vUpShapeVel[scvf].resize(m_numSh);
		m_vUpShapeDiffusion[scvf].resize(m_numSh);
	}
}



} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__ */
