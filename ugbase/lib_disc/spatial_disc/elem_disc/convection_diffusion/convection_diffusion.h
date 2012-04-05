/*
 * convection_diffusion.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/ip_data/data_import_export.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"
#include "lib_disc/spatial_disc/disc_util/finite_element_geometry.h"

namespace ug{

/// \ingroup lib_disc_elem_disc
/// @{

/// Discretization for the Convection-Diffusion Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation.
 * The Equation has the form
 * \f[
 * 	\partial_t (m1*c + m2) - \nabla \left( D \nabla c - \vec{v} c \right) +
 * 		r1 \cdot c + r2 = f
 * \f]
 * with
 * <ul>
 * <li>	\f$ c \f$ is the unknown solution
 * <li>	\f$ m1 \equiv m(\vec{x},t) \f$ is the Mass Scaling Term
 * <li>	\f$ m2 \equiv m(\vec{x},t) \f$ is the Mass Term
 * <li>	\f$ D \equiv D(\vec{x},t) \f$ is the Diffusion Tensor
 * <li>	\f$ v \equiv \vec{v}(\vec{x},t) \f$ is the Velocity Field
 * <li>	\f$ r1 \equiv r(\vec{x},t) \f$ is the Reaction Rate
 * <li>	\f$ r2 \equiv r(\vec{x},t) \f$ is a Reaction Term
 * <li>	\f$ f \equiv f(\vec{x},t) \f$ is a Source Term
 * </ul>
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain>
class ConvectionDiffusion
: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	Own type
		typedef ConvectionDiffusion<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	///	Constructor
		ConvectionDiffusion(const char* functions, const char* subsets);

	///	set the upwind method
	/**
	 * This method sets the upwind method used to upwind the convection.
	 *
	 * \param	shapes		upwind method
	 */
		void set_upwind(SmartPtr<IConvectionShapes<dim> > shapes);

	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
	///	\{
		void set_diffusion(SmartPtr<IPData<MathMatrix<dim, dim>, dim> > user);
		void set_diffusion(number val);
#ifdef UG_FOR_LUA
		void set_diffusion(const char* fctName);
#endif
	///	\}

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_velocity(SmartPtr<IPData<MathVector<dim>, dim> > user);
		void set_velocity(number vel_x);
		void set_velocity(number vel_x, number vel_y);
		void set_velocity(number vel_x, number vel_y, number vel_z);
#ifdef UG_FOR_LUA
		void set_velocity(const char* fctName);
#endif
	/// \}

	///	sets the reaction rate
	/**
	 * This method sets the Reaction Rate. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction_rate(SmartPtr<IPData<number, dim> > user);
		void set_reaction_rate(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate(const char* fctName);
#endif
	///	\}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction(SmartPtr<IPData<number, dim> > user);
		void set_reaction(number val);
#ifdef UG_FOR_LUA
		void set_reaction(const char* fctName);
#endif
	///	\}

	///	sets the source / sink term
	/**
	 * This method sets the source/sink value. A zero value is assumed as
	 * default.
	 */
	///	\{
		void set_source(SmartPtr<IPData<number, dim> > user);
		void set_source(number val);
#ifdef UG_FOR_LUA
		void set_source(const char* fctName);
#endif
	///	\}

	///	sets mass scale
	/**
	 * This method sets the mass scale value. A value of 1.0 is assumed as
	 * default.
	 */
	///	\{
		void set_mass_scale(SmartPtr<IPData<number, dim> > user);
		void set_mass_scale(number val);
#ifdef UG_FOR_LUA
		void set_mass_scale(const char* fctName);
#endif
	///	\}

	///	sets mass
	/**
	 * This method sets the mass value. A value of 0.0 is assumed as
	 * default.
	 */
	///	\{
		void set_mass(SmartPtr<IPData<number, dim> > user);
		void set_mass(number val);
#ifdef UG_FOR_LUA
		void set_mass(const char* fctName);
#endif
	///	\}

	public:
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID);

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular);

	///	returns if hanging nodes are needed
		virtual bool use_hanging() const;

	///	sets the disc scheme
		void set_disc_scheme(const char* c_scheme);

	///	sets the quad order for fe / fv
		void set_quad_order(size_t order)
			{m_quadOrder = order; m_quadOrderSCV = order;
			 m_quadOrderSCVF = order;m_bQuadOrderUserDef = true;}

	///	sets the quad order for scv of fv
		void set_quad_order_scv(size_t order)
			{m_quadOrderSCV = order; m_bQuadOrderUserDef = true;}

	///	sets the quad order for fe / fv
		void set_quad_order_scvf(size_t order)
			{m_quadOrderSCVF = order; m_bQuadOrderUserDef = true;}

	protected:
	///	sets the requested assembling routines
		void set_ass_funcs();

	///	current type of disc scheme
		std::string m_discScheme;

	///	current order of disc scheme
		int m_order;

	///	current shape function set
		LFEID m_lfeID;

	///	current integration order
		bool m_bQuadOrderUserDef;
		int m_quadOrder;
		int m_quadOrderSCV;
		int m_quadOrderSCVF;

	///	current regular grid flag
		bool m_bNonRegularGrid;

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

		/////////////////////////////////////
		//	Finite Volume assemblings (FV1)
		/////////////////////////////////////

	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void elem_loop_prepare_fv1();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void elem_prepare_fv1(TElem* elem, const LocalVector& u);

	///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void elem_loop_finish_fv1();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void ass_JA_elem_fv1(LocalMatrix& J, const LocalVector& u);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void ass_JM_elem_fv1(LocalMatrix& J, const LocalVector& u);

	///	assembles the stiffness part of the local defect
		template <typename TElem, typename TFVGeom>
		void ass_dA_elem_fv1(LocalVector& d, const LocalVector& u);

	///	assembles the mass part of the local defect
		template <typename TElem, typename TFVGeom>
		void ass_dM_elem_fv1(LocalVector& d, const LocalVector& u);

	///	assembles the local right hand side
		template <typename TElem, typename TFVGeom>
		void ass_rhs_elem_fv1(LocalVector& d);

		/////////////////////////////////////
		//	Finite Volume assemblings (FVHO)
		/////////////////////////////////////

	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void elem_loop_prepare_fvho();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void elem_prepare_fvho(TElem* elem, const LocalVector& u);

	///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void elem_loop_finish_fvho();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void ass_JA_elem_fvho(LocalMatrix& J, const LocalVector& u);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void ass_JM_elem_fvho(LocalMatrix& J, const LocalVector& u);

	///	assembles the stiffness part of the local defect
		template <typename TElem, typename TFVGeom>
		void ass_dA_elem_fvho(LocalVector& d, const LocalVector& u);

	///	assembles the mass part of the local defect
		template <typename TElem, typename TFVGeom>
		void ass_dM_elem_fvho(LocalVector& d, const LocalVector& u);

	///	assembles the local right hand side
		template <typename TElem, typename TFVGeom>
		void ass_rhs_elem_fvho(LocalVector& d);

		/////////////////////////////////////
		//	Finite Element assemblings
		/////////////////////////////////////

		template<typename TElem, typename TGeomProvider>
		void elem_loop_prepare_fe();

		template<typename TElem, typename TGeomProvider>
		void elem_prepare_fe(TElem* elem, const LocalVector& u);

		template<typename TElem, typename TGeomProvider>
		void elem_loop_finish_fe();

		template<typename TElem, typename TGeomProvider>
		void ass_JA_elem_fe(LocalMatrix& J, const LocalVector& u);

		template<typename TElem, typename TGeomProvider>
		void ass_JM_elem_fe(LocalMatrix& J, const LocalVector& u);

		template<typename TElem, typename TGeomProvider>
		void ass_dA_elem_fe(LocalVector& d, const LocalVector& u);

		template<typename TElem, typename TGeomProvider>
		void ass_dM_elem_fe(LocalVector& d, const LocalVector& u);

		template<typename TElem, typename TGeomProvider>
		void ass_rhs_elem_fe(LocalVector& d);

	protected:
	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_velocity_fv1(const LocalVector& u,
		                          std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                          const size_t nip);

	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_diffusion_fv1(const LocalVector& u,
		                           std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFVGeom>
		void lin_def_reaction_fv1(const LocalVector& u,
		                          std::vector<std::vector<number> > vvvLinDef[],
		                          const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFVGeom>
		void lin_def_reaction_rate_fv1(const LocalVector& u,
		                               std::vector<std::vector<number> > vvvLinDef[],
		                               const size_t nip);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem, typename TFVGeom>
		void lin_def_source_fv1(const LocalVector& u,
		                        std::vector<std::vector<number> > vvvLinDef[],
		                        const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFVGeom>
		void lin_def_mass_scale_fv1(const LocalVector& u,
		                            std::vector<std::vector<number> > vvvLinDef[],
		                            const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFVGeom>
		void lin_def_mass_fv1(const LocalVector& u,
		                      std::vector<std::vector<number> > vvvLinDef[],
		                      const size_t nip);

	protected:
	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_velocity_fvho(const LocalVector& u,
		                          std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                          const size_t nip);

	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_diffusion_fvho(const LocalVector& u,
		                           std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the reaction rate
		template <typename TElem, typename TFVGeom>
		void lin_def_reaction_rate_fvho(const LocalVector& u,
		                                std::vector<std::vector<number> > vvvLinDef[],
		                                const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFVGeom>
		void lin_def_reaction_fvho(const LocalVector& u,
								  std::vector<std::vector<number> > vvvLinDef[],
								  const size_t nip);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem, typename TFVGeom>
		void lin_def_source_fvho(const LocalVector& u,
		                        std::vector<std::vector<number> > vvvLinDef[],
		                        const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFVGeom>
		void lin_def_mass_scale_fvho(const LocalVector& u,
		                            std::vector<std::vector<number> > vvvLinDef[],
		                            const size_t nip);

	///	computes the linearized defect w.r.t to the mass term
		template <typename TElem, typename TFVGeom>
		void lin_def_mass_fvho(const LocalVector& u,
		                       std::vector<std::vector<number> > vvvLinDef[],
		                       const size_t nip);

	protected:
	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TGeomProvider>
		void lin_def_velocity_fe(const LocalVector& u,
		                          std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                          const size_t nip);

	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TGeomProvider>
		void lin_def_diffusion_fe(const LocalVector& u,
		                           std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the reaction rate
		template <typename TElem, typename TGeomProvider>
		void lin_def_reaction_rate_fe(const LocalVector& u,
		                              std::vector<std::vector<number> > vvvLinDef[],
		                              const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TGeomProvider>
		void lin_def_reaction_fe(const LocalVector& u,
								  std::vector<std::vector<number> > vvvLinDef[],
								  const size_t nip);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem, typename TGeomProvider>
		void lin_def_source_fe(const LocalVector& u,
		                        std::vector<std::vector<number> > vvvLinDef[],
		                        const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TGeomProvider>
		void lin_def_mass_scale_fe(const LocalVector& u,
		                            std::vector<std::vector<number> > vvvLinDef[],
		                            const size_t nip);

	///	computes the linearized defect w.r.t to the mass term
		template <typename TElem, typename TGeomProvider>
		void lin_def_mass_fe(const LocalVector& u,
		                     std::vector<std::vector<number> > vvvLinDef[],
		                     const size_t nip);

	private:
	///	Corner Coordinates
		const position_type* m_vCornerCoords;

	///	abbreviation for the local solution
		static const size_t _C_ = 0;

	/// method to compute the upwind shapes
		SmartPtr<IConvectionShapes<dim> > m_spConvShape;

	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim > m_imVelocity;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReactionRate;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReaction;

	///	Data import for the right-hand side
		DataImport<number, dim> m_imSource;

	///	Data import for the mass scale
		DataImport<number, dim> m_imMassScale;

	///	Data import for the mass scale
		DataImport<number, dim> m_imMass;

	public:
		typedef SmartPtr<IPData<number, dim> > NumberExport;
		typedef SmartPtr<IPData<MathVector<dim>, dim> > GradExport;

	///	returns the export of the value of associated unknown function
		SmartPtr<IPData<number, dim> > value();

	///	returns the export of the gradient of associated unknown function
		SmartPtr<IPData<MathVector<dim>, dim> > gradient();

	protected:
		typedef IConvectionShapes<dim> conv_shape_type;

	///	returns the updated convection shapes
		const IConvectionShapes<dim>& get_updated_conv_shapes(const FVGeometryBase& geo);

	///	computes the concentration
		template <typename TElem, typename TFVGeom>
		void ex_value_fv1(const LocalVector& u,
		                  const MathVector<dim> vGlobIP[],
		                  const MathVector<TFVGeom::dim> vLocIP[],
		                  const size_t nip,
		                  number vValue[],
		                  bool bDeriv,
		                  std::vector<std::vector<number> > vvvDeriv[]);

	///	computes the gradient of the concentration
		template <typename TElem, typename TFVGeom>
		void ex_grad_fv1(const LocalVector& u,
		                 const MathVector<dim> vGlobIP[],
		                 const MathVector<TFVGeom::dim> vLocIP[],
		                 const size_t nip,
		                 MathVector<dim> vValue[],
		                 bool bDeriv,
		                 std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	///	computes the concentration
		template <typename TElem, typename TGeomProvider>
		void ex_value_fvho(const LocalVector& u,
		                   const MathVector<dim> vGlobIP[],
		                   const MathVector<TGeomProvider::Type::dim> vLocIP[],
		                   const size_t nip,
		                   number vValue[],
		                   bool bDeriv,
		                   std::vector<std::vector<number> > vvvDeriv[]);

	///	computes the gradient of the concentration
		template <typename TElem, typename TGeomProvider>
		void ex_grad_fvho(const LocalVector& u,
		                  const MathVector<dim> vGlobIP[],
		                  const MathVector<TGeomProvider::Type::dim> vLocIP[],
		                  const size_t nip,
		                  MathVector<dim> vValue[],
		                  bool bDeriv,
		                  std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	///	computes the concentration
		template <typename TElem, typename TGeomProvider>
		void ex_value_fe(const LocalVector& u,
		                 const MathVector<dim> vGlobIP[],
		                 const MathVector<TGeomProvider::Type::dim> vLocIP[],
		                 const size_t nip,
		                 number vValue[],
		                 bool bDeriv,
		                 std::vector<std::vector<number> > vvvDeriv[]);

	///	computes the gradient of the concentration
		template <typename TElem, typename TGeomProvider>
		void ex_grad_fe(const LocalVector& u,
		                const MathVector<dim> vGlobIP[],
		                const MathVector<TGeomProvider::Type::dim> vLocIP[],
		                const size_t nip,
		                MathVector<dim> vValue[],
		                bool bDeriv,
		                std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exConcentration;

	///	Export for the gradient of concentration
		SmartPtr<DataExport<MathVector<dim>, dim> > m_exConcentrationGrad;

	protected:
	// 	FV1 Assemblings
		void register_all_fv1_funcs(bool bHang);

		template <typename TElem, typename TFVGeom>
		void register_fv1_func();


	// 	FVHO Assemblings
		void register_all_fvho_funcs(int order, int quadOrderSCV, int quadOrderSCVF);

		template<typename TElem, typename TGeomProvider>
		void register_fvho_func();


	// 	FE Assemblings
		void register_all_fe_funcs(int order, int quadOrder);

		template <typename TElem, typename TGeomProvider>
		void register_fe_func();

	//	helper class holding a geometry
		template<typename TGeom>
		struct FlexGeomProvider
		{
			typedef TGeom Type;
			static inline TGeom& get(){static TGeom inst; return inst;}
		};
};

/// @}

} // end namespace ug


#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__*/
