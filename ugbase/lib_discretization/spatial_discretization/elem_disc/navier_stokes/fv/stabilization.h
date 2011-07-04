/*
 * stabilization.h
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__

#include "upwind.h"

namespace ug{

enum DIFFUSION_LENGTH
{
    NS_RAW = 0,
    NS_FIVEPOINT,
    NS_COR
};

template <int dim>
class INavierStokesStabilization
{
	public:
	/// Local vector type
		typedef LocalVector local_vector_type;

	/// Abbreviation for own type
		typedef INavierStokesStabilization<dim> this_type;

	///	type of update function
		typedef bool (this_type::*UpdateFunc)(	const FVGeometryBase* geo,
												const local_vector_type& vCornerValue,
												const DataImport<number, dim>& kinVisco,
									            const DataImport<MathVector<dim>, dim>* pSource,
												const local_vector_type* pvCornerValueOldTime, number dt);

	public:
	///	constructor
		INavierStokesStabilization()
			: m_numScvf(0), m_numSh(0), m_pUpwind(NULL)
		{
			m_vUpdateFunc.clear();
			m_vDiffLengthSqInv.clear();

		//	default setup
			set_diffusion_length("NS_RAW");
		}

	///	sets the type of diff length used for evaluation
		bool set_diffusion_length(std::string diffLength);

	///	sets the upwind method
		void set_upwind(INavierStokesUpwind<dim>& upwind)
		{
			m_pUpwind = &upwind;
			m_pConstUpwind = const_cast<const INavierStokesUpwind<dim>*>(m_pUpwind);
		}

	///	returns the upwind
		const INavierStokesUpwind<dim>* get_upwind() const
			{return m_pConstUpwind;}

	///	diff length
		number diff_length_sq_inv(size_t scvf) const
		{
			UG_ASSERT(scvf < m_vDiffLengthSqInv.size(), "Invalid index");
			return m_vDiffLengthSqInv[scvf];
		}

	/// stabilized velocity
		const MathVector<dim>& stab_vel(size_t scvf) const
		{
			UG_ASSERT(scvf < m_vStabVel.size(), "Invalid index");
			return m_vStabVel[scvf];
		}

	///	returns if stab velocity comp depends on other vel components
		bool vel_comp_connected() const {return m_bVelCompConnected;}

	/// computed stab shape for velocity. This is: The stab_vel derivative
	/// w.r.t velocity unknowns in the corner for each component
		number stab_shape_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh) const
		{
			UG_ASSERT(scvf < m_vvvvStabShapeVel.size(), "Invalid index.");
			UG_ASSERT(compOut < m_vvvvStabShapeVel[scvf].size(), "Invalid index.");
			UG_ASSERT(compIn < m_vvvvStabShapeVel[scvf][compOut].size(), "Invalid index.");
			UG_ASSERT(sh < m_vvvvStabShapeVel[scvf][compOut][compIn].size(), "Invalid index.");
			return m_vvvvStabShapeVel[scvf][compOut][compIn][sh];
		}

	///	computed stab shape for pressure.
		number stab_shape_p(size_t scvf, size_t compOut, size_t sh) const
		{
			UG_ASSERT(scvf < m_vvvvStabShapePressure.size(), "Invalid index.");
			UG_ASSERT(compOut < m_vvvvStabShapePressure[scvf].size(), "Invalid index.");
			UG_ASSERT(sh < m_vvvvStabShapePressure[scvf][compOut].size(), "Invalid index.");
			return m_vvvvStabShapePressure[scvf][compOut][sh];
		}


	///	compute values for new geometry and corner velocities
		bool update(const FVGeometryBase* geo,
		            const local_vector_type& vCornerValue,
		            const DataImport<number, dim>& kinVisco,
		            const DataImport<MathVector<dim>, dim>* pSource,
		            const local_vector_type* pvCornerValueOldTime, number dt)
			{return (this->*(m_vUpdateFunc[m_id]))(	geo, vCornerValue,
													kinVisco, pSource,
													pvCornerValueOldTime, dt);}

	/////////////////////////////////////////
	// forward methods of Upwind Velocity
	/////////////////////////////////////////

	///	Convection Length
		number conv_length(size_t scvf) const
		{
			UG_ASSERT(m_pConstUpwind != NULL, "No upwind object");
			return m_pConstUpwind->conv_length(scvf);
		}

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_ASSERT(m_pConstUpwind != NULL, "No upwind object");
			return m_pConstUpwind->upwind_shape_sh(scvf, sh);
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const
		{
			UG_ASSERT(m_pConstUpwind != NULL, "No upwind object");
			return m_pConstUpwind->non_zero_shape_ip();
		}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_ASSERT(m_pConstUpwind != NULL, "No upwind object");
			return m_pConstUpwind->upwind_shape_ip(scvf, scvf2);
		}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	resize the data arrays
		void set_sizes(size_t numScvf, size_t numSh);

	///	Upwind values
		INavierStokesUpwind<dim>* m_pUpwind;
		const INavierStokesUpwind<dim>* m_pConstUpwind;

	///	computes the diffusion length
		template <typename TFVGeom>
		bool compute_upwind(const TFVGeom& geo, const local_vector_type& vCornerValue);

	///	type of diffusion length computation
		int m_diffLengthType;

	///	vector holding diffusion Length squared and inverted
		std::vector<number> m_vDiffLengthSqInv;

	///	computes the diffusion length
		template <typename TFVGeom>
		bool compute_diff_length(const TFVGeom& geo);

	///	values of stabilized velocity at ip
		std::vector<MathVector<dim> > m_vStabVel;

	///	flag if velocity components are interconnected
		bool m_bVelCompConnected;

	///	sets the vel comp connected flag
		void set_vel_comp_connected(bool bVelCompConnected) {m_bVelCompConnected = bVelCompConnected;}

	///	stab shapes w.r.t vel
		std::vector<std::vector<std::vector<std::vector<number> > > > m_vvvvStabShapeVel;

	///	stab shapes w.r.t pressure
		std::vector<std::vector<std::vector<number> > > m_vvvvStabShapePressure;

	/// stabilized velocity
		MathVector<dim>& stab_vel(size_t scvf)
		{
			UG_ASSERT(scvf < m_vStabVel.size(), "Invalid index");
			return m_vStabVel[scvf];
		}

	/// computed stab shape for velocity. This is: The stab_vel derivative
	/// w.r.t velocity unknowns in the corner for each component
		number& stab_shape_vel(size_t scvf, size_t compOut, size_t compIn, size_t sh)
		{
			UG_ASSERT(scvf < m_vvvvStabShapeVel.size(), "Invalid index.");
			UG_ASSERT(compOut < m_vvvvStabShapeVel[scvf].size(), "Invalid index.");
			UG_ASSERT(compIn < m_vvvvStabShapeVel[scvf][compOut].size(), "Invalid index.");
			UG_ASSERT(sh < m_vvvvStabShapeVel[scvf][compOut][compIn].size(), "Invalid index.");
			return m_vvvvStabShapeVel[scvf][compOut][compIn][sh];
		}

	///	computed stab shape for pressure.
		number& stab_shape_p(size_t scvf, size_t compOut, size_t sh)
		{
			UG_ASSERT(scvf < m_vvvvStabShapePressure.size(), "Invalid index.");
			UG_ASSERT(compOut < m_vvvvStabShapePressure[scvf].size(), "Invalid index.");
			UG_ASSERT(sh < m_vvvvStabShapePressure[scvf][compOut].size(), "Invalid index.");
			return m_vvvvStabShapePressure[scvf][compOut][sh];
		}

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
// FIELDS
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesFIELDSStabilization
	: public INavierStokesStabilization<TDim>
{
	public:
	///	Base class
		typedef INavierStokesStabilization<TDim> base_type;

	///	This class
		typedef NavierStokesFIELDSStabilization<TDim> this_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::register_update_func;
		using base_type::diff_length_sq_inv;
		using base_type::stab_shape_vel;
		using base_type::stab_shape_p;
		using base_type::stab_vel;
		using base_type::set_vel_comp_connected;

	//	functions from upwind
		using base_type::conv_length;
		using base_type::upwind_shape_sh;
		using base_type::non_zero_shape_ip;
		using base_type::upwind_shape_ip;

	public:
	///	constructor
		NavierStokesFIELDSStabilization()
		{
		//	vel comp not interconnected
			set_vel_comp_connected(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo,
		            const local_vector_type& vCornerValue,
		            const DataImport<number, dim>& kinVisco,
		            const DataImport<MathVector<dim>, dim>* pSource,
		            const local_vector_type* pvCornerValueOldTime, number dt);

	private:
		void register_func(Int2Type<1>)
		{register_func<Edge>();}

		void register_func(Int2Type<2>)
		{	register_func(Int2Type<1>());
			register_func<Triangle>();
			register_func<Quadrilateral>();}

		void register_func(Int2Type<3>)
		{	register_func(Int2Type<2>());
			register_func<Tetrahedron>();
			register_func<Pyramid>();
			register_func<Prism>();
			register_func<Hexahedron>();}

		template <typename TElem>
		void register_func()
		{
			typedef FV1Geometry<TElem, dim> TGeom;
			typedef bool (this_type::*TFunc)(const TGeom* geo,
											 const local_vector_type& vCornerValue,
											 const DataImport<number, dim>& kinVisco,
											 const DataImport<MathVector<dim>, dim>* pSource,
											 const local_vector_type* pvCornerValueOldTime, number dt);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};


} // end namespace ug

// include implementation
#include "stabilization_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__STABILIZATION__ */
