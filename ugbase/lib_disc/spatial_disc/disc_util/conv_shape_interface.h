/*
 * conv_shape_interface.h
 *
 *  Created on: 08.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__

#include "common/common.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "common/util/provider.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// Interface for Convection shapes
/////////////////////////////////////////////////////////////////////////////


/**
 * The idea of convection shapes is to generalize the upwinding, such that
 * the upwind may also depend on additional parameters.
 *
 * In the case of d3f for example, one is interested in the evaluation of
 * the convective flux at an integration point
 *
 * \f[ conv_flux_{ip} := (\omega \vec{q} \cdot \vec{n})|_{ip} \f]
 *
 * In order to do so, the flux is interpolated by shape functions as:
 * 
 * \f[ conv_flux_{ip} := \sum_{sh=1}^{n_sh} \phi_{sh}^{conv}|_{ip} \cdot \omega_{sh} \f]
 *
 * The computation of the convective flux does not only depend on the geometry
 * but does also depend on the convection velocity and the Diffusion-Dispersion-
 * Tensor in order to allow a weighting between diffusive and convective flux.
 * Therefore, also the convection shape depend on those quantities and the
 * derivatives of the convection shapes w.r.t. the velocity and the Diff-Disp-Tensor
 * must be taken into account when computing an exact jacobian.
 *
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
				 const MathVector<dim>* DarcyVelocity,
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
					const MathVector<dim>* DarcyVelocity,
					const MathMatrix<dim, dim>* DiffDisp,
		            bool computeDeriv)
			{return (this->*(m_vUpdateFunc[m_id]))(	geo, DarcyVelocity,DiffDisp, computeDeriv);}

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
