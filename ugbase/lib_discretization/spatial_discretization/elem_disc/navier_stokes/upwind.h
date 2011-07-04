/*
 * upwind.h
 *
 *  Created on: 10.03.2011
 *      Author: andreasvogel
 */

#ifndef NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__
#define NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__

#include <vector>

#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/disc_util/finite_volume_geometry.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////
// Interface for Upwinds
/////////////////////////////////////////////////////////////////////////////


template <int dim>
class INavierStokesUpwind
{
	public:
	/// Local vector type
		typedef LocalVector local_vector_type;

	/// Abbreviation for own type
		typedef INavierStokesUpwind<dim> this_type;

	///	type of update function
		typedef bool (this_type::*UpdateFunc)(	const FVGeometryBase* obj,
												const local_vector_type& vCornerVels);

	public:
	///	constructor
		INavierStokesUpwind()
			: m_pCornerValue(NULL), m_numScvf(0), m_numSh(0), m_bNonZeroShapeIp(true)
		{
			m_vUpdateFunc.clear();
			m_vConvLength.clear();
			m_vIPVel.clear();
			m_vUpShapeSh.clear();
			m_vUpShapeIp.clear();
		}

	///	returns number of shapes
		size_t num_sh() const {return m_numSh;}

	///	returns number of sub control volume faces
		size_t num_scvf() const {return m_numScvf;}

	///	Convection Length
		number conv_length(size_t scvf) const
		{
			UG_ASSERT(scvf < m_vConvLength.size(), "Invalid index");
			return m_vConvLength[scvf];
		}

	///	ip velocity (i.e. interpolated velocity at ip)
		const MathVector<dim>& ip_vel(size_t scvf) const
		{
			UG_ASSERT(scvf < m_vIPVel.size(), "Invalid index");
			return m_vIPVel[scvf];
		}

	/// returns the upwind velocity
		MathVector<dim> upwind_vel(size_t scvf) const;

	///	upwind shape for corner vel
		number upwind_shape_sh(size_t scvf, size_t sh) const
		{
			UG_ASSERT(scvf < m_vUpShapeSh.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeSh[scvf].size(), "Invalid index");
			return m_vUpShapeSh[scvf][sh];
		}

	///	returns if upwind shape w.r.t. ip vel is non-zero
		bool non_zero_shape_ip() const {return m_bNonZeroShapeIp;}

	///	upwind shapes for ip vel
		number upwind_shape_ip(size_t scvf, size_t scvf2) const
		{
			UG_ASSERT(scvf < m_vUpShapeIp.size(), "Invalid index");
			UG_ASSERT(scvf2 < m_vUpShapeIp[scvf].size(), "Invalid index");
			return m_vUpShapeIp[scvf][scvf2];
		}

	///	compute values for new geometry and corner velocities
		bool update(const FVGeometryBase* geo, const local_vector_type& vCornerValue)
			{	m_pCornerValue = &vCornerValue;
				return (this->*(m_vUpdateFunc[m_id]))(geo, vCornerValue);}

	//////////////////////////
	// internal handling
	//////////////////////////

	protected:
	///	resize the data arrays
		void set_sizes(size_t numScvf, size_t numSh);

	///	sets the shape ip flag
		void set_shape_ip_flag(bool flag) {m_bNonZeroShapeIp = flag;}

	/// non-const access to ip velocity (i.e. interpolated velocity at ip)
		MathVector<dim>& ip_vel(size_t scvf)
		{
			UG_ASSERT(scvf < m_vIPVel.size(), "Invalid index");
			return m_vIPVel[scvf];
		}

	///	non-const access to upwind shapes for corner vel
		number& upwind_shape_sh(size_t scvf, size_t sh)
		{
			UG_ASSERT(scvf < m_vUpShapeSh.size(), "Invalid index");
			UG_ASSERT(sh < m_vUpShapeSh[scvf].size(), "Invalid index");
			return m_vUpShapeSh[scvf][sh];
		}

	///	non-const access to upwind shapes for ip vel
		number& upwind_shape_ip(size_t scvf, size_t scvf2)
		{
			UG_ASSERT(scvf < m_vUpShapeIp.size(), "Invalid index");
			UG_ASSERT(scvf2 < m_vUpShapeIp[scvf].size(), "Invalid index");
			return m_vUpShapeIp[scvf][scvf2];
		}

	///	non-const access to Convection Length
		number& conv_length(size_t scvf)
		{
			UG_ASSERT(scvf < m_vConvLength.size(), "Invalid index");
			return m_vConvLength[scvf];
		}

	///	pointer to currently used values
		const local_vector_type* m_pCornerValue;

	///	interpolated value at ip
		std::vector<MathVector<dim> > m_vIPVel;

	///	number of current scvf
		size_t m_numScvf;

	///	number of current shape functions (usually in corners)
		size_t m_numSh;

	///	convection length
		std::vector<number> m_vConvLength;

	///	upwind shapes for corners shape functions
		std::vector<std::vector<number> > m_vUpShapeSh;

	///	flag if ip shapes are non-zero
		bool m_bNonZeroShapeIp;

	///	upwind shapes for ip vels
		std::vector<std::vector<number> > m_vUpShapeIp;

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
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesNoUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesNoUpwind<TDim> this_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::upwind_shape_sh;
		using base_type::upwind_shape_ip;
		using base_type::conv_length;
		using base_type::ip_vel;
		using base_type::register_update_func;

	public:
	///	constructor
		NavierStokesNoUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue);

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
			typedef bool (this_type::*TFunc)(const TGeom* geo, const local_vector_type& vCornerValue);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesFullUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesFullUpwind<TDim> this_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::upwind_shape_sh;
		using base_type::upwind_shape_ip;
		using base_type::conv_length;
		using base_type::ip_vel;
		using base_type::register_update_func;

	public:
	///	constructor
		NavierStokesFullUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue);

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
			typedef bool (this_type::*TFunc)(const TGeom* geo, const local_vector_type& vCornerValue);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesSkewedUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesSkewedUpwind<TDim> this_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::upwind_shape_sh;
		using base_type::upwind_shape_ip;
		using base_type::conv_length;
		using base_type::ip_vel;
		using base_type::register_update_func;

	public:
	///	constructor
		NavierStokesSkewedUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue);

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
			typedef bool (this_type::*TFunc)(const TGeom* geo, const local_vector_type& vCornerValue);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

/////////////////////////////////////////////////////////////////////////////
// Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesLinearProfileSkewedUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesLinearProfileSkewedUpwind<TDim> this_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::upwind_shape_sh;
		using base_type::upwind_shape_ip;
		using base_type::conv_length;
		using base_type::ip_vel;
		using base_type::register_update_func;

	public:
	///	constructor
		NavierStokesLinearProfileSkewedUpwind()
		{
		//	shapes for ip vels are zero (no dependency between ip shapes, only to corners)
		//	Note: during resize, values are initialized to zero, thus values
		//		  for shapes depending on ip values are always correct (i.e. zero)
			set_shape_ip_flag(false);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue);

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
			typedef bool (this_type::*TFunc)(const TGeom* geo, const local_vector_type& vCornerValue);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};


/////////////////////////////////////////////////////////////////////////////
// Positive Upwind
/////////////////////////////////////////////////////////////////////////////

template <int TDim>
class NavierStokesPositiveUpwind
	: public INavierStokesUpwind<TDim>
{
	public:
	///	Base class
		typedef INavierStokesUpwind<TDim> base_type;

	///	This class
		typedef NavierStokesPositiveUpwind<TDim> this_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Dimension
		static const int dim = TDim;

	protected:
	//	explicitly forward some function
		using base_type::set_shape_ip_flag;
		using base_type::upwind_shape_sh;
		using base_type::upwind_shape_ip;
		using base_type::conv_length;
		using base_type::ip_vel;
		using base_type::register_update_func;

	public:
	///	constructor
		NavierStokesPositiveUpwind()
		{
		//	shapes for ip vels are non-zero (dependency between ip shapes)
			set_shape_ip_flag(true);

		//	register evaluation function
			register_func(Int2Type<dim>());
		}

	///	update of values for FV1Geometry
		template <typename TElem>
		bool update(const FV1Geometry<TElem, dim>* geo, const local_vector_type& vCornerValue);

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
			typedef bool (this_type::*TFunc)(const TGeom* geo, const local_vector_type& vCornerValue);

			this->template register_update_func<TGeom, TFunc>(&this_type::template update<TElem>);
		}
};

} // end namespace ug

// include implementation
#include "upwind_impl.h"

#endif /* NEW_STABILIZATION_IMPL_H___H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__UPWIND__ */
