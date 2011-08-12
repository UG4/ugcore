/*
 * constant_equation_impl.h
 *
 *  Created on: 13.05.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONSTANT_EQUATION__FV1__CONSTANT_EQUATION__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONSTANT_EQUATION__FV1__CONSTANT_EQUATION__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/ip_data/data_import_export.h"

namespace ug{

/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Volume Element Discretization for a Constant Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the equation:
 * The Equation has the form
 * \f[
 * 	\partial_t m - \nabla \left( \vec{v} \right) = f
 * \f]
 * with
 * <ul>
 * <li>	<b>no</b> unknown solution
 * <li>	\f$ m \equiv m(\vec{x},t) \f$ is the Mass Scaling Term
 * <li>	\f$ v \equiv \vec{v}(\vec{x},t) \f$ is the Velocity Field
 * <li>	\f$ f \equiv f(\vec{x},t) \f$ is a Source Term
 * </ul>
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain>
class FVConstantEquationElemDisc : public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	Own type
		typedef FVConstantEquationElemDisc<TDomain> this_type;

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
		FVConstantEquationElemDisc();

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
		void set_velocity(IPData<MathVector<dim>, dim>& user) {m_imVelocity.set_data(user);}

	///	sets the source / sink term
	/**
	 * This method sets the source/sink value. A zero value is assumed as
	 * default.
	 */
		void set_source(IPData<number, dim>& user)	{m_imSource.set_data(user);}

	///	sets mass scale
	/**
	 * This method sets the mass scale value. A value of 1.0 is assumed as
	 * default.
	 */
		void set_mass_scale(IPData<number, dim>& user)	{m_imMassScale.set_data(user);}

	public:
	///	number of functions used
		virtual size_t num_fct(){return 1;}

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
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool prepare_element_loop();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool prepare_element(TElem* elem, const local_vector_type& u);

	///	finishes the loop over all elements
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u);

	///	assembles the stiffness part of the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool assemble_A(local_vector_type& d, const local_vector_type& u);

	///	assembles the mass part of the local defect
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool assemble_M(local_vector_type& d, const local_vector_type& u);

	///	assembles the local right hand side
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool assemble_f(local_vector_type& d);

	protected:
	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool lin_defect_velocity(const local_vector_type& u);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool lin_defect_source(const local_vector_type& u);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool lin_defect_mass_scale(const local_vector_type& u);

	private:
	///	Corner Coordinates
		const position_type* m_vCornerCoords;

	///	abbreviation for the local solution
		static const size_t _C_ = 0;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim > m_imVelocity;

	///	Data import for the right-hand side
		DataImport<number, dim> m_imSource;

	///	Data import for the mass scale
		DataImport<number, dim> m_imMassScale;

	public:
	///	returns the export of the concentration
		IPData<number, dim>& get_concentration() {return m_exConcentration;}

	///	returns the export of gradient of the concentration
		IPData<MathVector<dim>, dim>& get_concentration_grad() {return m_exConcentrationGrad;}

	protected:
	///	computes the concentration
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool ex_concentration(const local_vector_type& u,
							 const MathVector<dim>* vGlobIP,
							 const MathVector<TFVGeom<TElem, dim>::dim>* vLocIP,
							 const size_t nip,
							 number* vValue,
							 bool bDeriv,
							 std::vector<std::vector<std::vector<number> > >& vvvDeriv);

	///	computes the gradient of the concentration
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		bool ex_concentration_grad(const local_vector_type& u,
								  const MathVector<dim>* vGlobIP,
								  const MathVector<TFVGeom<TElem, dim>::dim>* vLocIP,
								  const size_t nip,
								  MathVector<dim>* vValue,
								  bool bDeriv,
								  std::vector<std::vector<std::vector<MathVector<dim> > > >& vvvDeriv);

	///	Export for the concentration
		DataExport<number, dim> m_exConcentration;

	///	Export for the gradient of concentration
		DataExport<MathVector<dim>, dim> m_exConcentrationGrad;

	private:
		void register_all_fv1_funcs(bool bHang);

		template <template <class Elem, int WorldDim> class TFVGeom>
		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fv1_func<TElem, TFVGeom>();}
		};

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void register_fv1_func();

};

/// @}

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONSTANT_EQUATION__FV1__CONSTANT_EQUATION__*/
