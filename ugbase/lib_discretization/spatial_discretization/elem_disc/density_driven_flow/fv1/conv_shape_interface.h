/*
 * conv_shape_interface.h
 *
 *  Created on: 08.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// Interface for Convection shapes
/////////////////////////////////////////////////////////////////////////////

template <int dim, typename TAlgebra>
class IConvectionShapes
{
	public:
	/// Abbreviation for own type
		typedef IConvectionShapes<dim, TAlgebra> this_type;

	///	type of update function
		typedef bool (this_type::*UpdateFunc)
				(const FVGeometryBase* geo,
				 const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
				 const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
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
		number conv_shape(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShape.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShape[scvf].size(), "Invalid index");
			return m_vUpShape[scvf][sh];
		}

	///	upwind shape for corner vel
		const MathVector<dim>& conv_shape_vel(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShapeVel.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeVel[scvf].size(), "Invalid index");
			return m_vUpShapeVel[scvf][sh];
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_deriv_diffusion() const {return m_bNonZeroDerivDiffusion;}

	///	upwind shapes for ip vel
		const MathMatrix<dim,dim>& conv_shape_diffusion(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShapeDiffusion.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeDiffusion[scvf].size(), "Invalid index");
			return m_vUpShapeDiffusion[scvf][sh];
		}

		bool update(const FVGeometryBase* geo,
		            const DataImport<MathVector<dim>, dim, TAlgebra>& DarcyVelocity,
		            const DataImport<MathMatrix<dim, dim>, dim, TAlgebra>& Diffusion,
		            bool computeDeriv)
			{return (this->*(m_vUpdateFunc[m_id]))(	geo, DarcyVelocity,Diffusion, computeDeriv);}

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
			return m_vUpShape[scvf][sh];
		}

	///	non-const access to upwind shapes for corner vel
		MathVector<dim>& conv_shape_vel(size_t scvf, size_t sh)
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
		bool set_geometry_type();

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
template <int dim, typename TAlgebra>
template <typename TFVGeom, typename TAssFunc>
void
IConvectionShapes<dim, TAlgebra>::
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
template <int dim, typename TAlgebra>
template <typename TFVGeom>
bool
IConvectionShapes<dim, TAlgebra>::
set_geometry_type()
{
//	get unique geometry id
	size_t id = GetUniqueFVGeomID<TFVGeom>();

//	check that function exists
	if(id >= m_vUpdateFunc.size() || m_vUpdateFunc[id] == NULL)
	{
		UG_LOG("ERROR in 'INavierStokesStabilization::set_geometry_type':"
				" No update function registered for this Geometry.\n");
		return false;
	}

//	set current geometry
	m_id = id;

//	set sizes
	TFVGeom& geo = FVGeometryProvider::get_geom<TFVGeom>();
	set_sizes(geo.num_scvf(), geo.num_scv());

//	we're done
	return true;
}


//	resize the data arrays
template <int dim, typename TAlgebra>
void
IConvectionShapes<dim, TAlgebra>::
set_sizes(size_t numScvf, size_t numSh)
{
//	remember sizes
	m_numScvf = numScvf;
	m_numSh = numSh;

//	adjust arrays
	m_vUpShape.resize(m_numScvf);
	m_vUpShapeVel.resize(m_numScvf);
	m_vUpShapeVel.resize(m_numScvf);
	for(size_t scvf = 0; scvf < m_numScvf; ++scvf)
	{
		m_vUpShape[scvf].resize(m_numSh);
		m_vUpShapeVel[scvf].resize(m_numSh);
		m_vUpShapeDiffusion[scvf].resize(m_numSh);
	}
}



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__CONV_SHAPE_INTERFACE__ */
