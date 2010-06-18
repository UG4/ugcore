
#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_ASSEMBLE__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_ASSEMBLE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/domain_discretization/disc_helper/fvgeom.h"

#include "lib_discretization/domain_discretization/plug_in_disc/plug_in_element_disc_interface.h"
#include "lib_discretization/domain_discretization/disc_coupling/element_data.h"

namespace ug{

template<typename TDomain, typename TAlgebra, typename TElem >
class DensityDrivenFlow : public DataExportingClass<MathVector<TDomain::dim>, MathVector<TDomain::dim>, TAlgebra>{
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	public:
		// forward constants and types

		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = TDomain::dim;
		static const int ref_dim = ref_elem_type::dim;

		// position type
		typedef typename TDomain::position_type position_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef typename algebra_type::matrix_type::local_matrix_type local_matrix_type;

		// local vector type
		typedef typename algebra_type::vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename algebra_type::vector_type::local_index_type local_index_type;

	protected:
		typedef void (*Pososity_fct)(number&);
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);
		typedef void (*Mol_Diff_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Permeability_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Gravity_fct)(MathVector<dim>&);

	public:
		DensityDrivenFlow(TDomain& domain, number upwind_amount, Pososity_fct Porosity, Viscosity_fct Viscosity, Density_fct Density, D_Density_fct D_Density,
				Mol_Diff_Tensor_fct Mol_Diff, Permeability_Tensor_fct Permeability_Tensor, Gravity_fct Gravity, DataClassExportPossibility<MathVector<dim>, MathVector<dim>,TAlgebra>& Darcy_Velocity_export)
		: m_domain(domain), m_upwind_amount(upwind_amount),
			m_Porosity(Porosity), m_Viscosity(Viscosity), m_Density(Density), m_D_Density(D_Density),
			m_Mol_Diff_Tensor(Mol_Diff), m_Permeability_Tensor(Permeability_Tensor), m_Gravity(Gravity), m_Darcy_Velocity_export(Darcy_Velocity_export)
		{};

	public:

		// total number of shape functions on elements of type 'TElem'
		inline size_t num_sh(){return 2*ref_elem_type::num_corners;};

		// number of shape functions on elements of type 'TElem' for the 'i'-th fundamental function
		inline size_t num_sh(size_t loc_fct){return ref_elem_type::num_corners;};

		// prepares the loop. Must be called, before prepare_element can be used
		inline IPlugInReturn prepare_element_loop();

		// prepares the element. Must be called before assemble_element_XXX can be used. Must be called after prepare_element_loop
		inline IPlugInReturn prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		inline IPlugInReturn assemble_element_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_JABnd(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_ABnd(local_vector_type& d, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_f(local_vector_type& d, number time=0.0);

		inline IPlugInReturn finish_element_loop();

	private:
		// domain
		TDomain& m_domain;

		// position access
		ReferenceMapping<ref_elem_type, TDomain::dim> m_mapping;
		typename TDomain::position_type m_corners[ref_elem_type::num_corners];
		typename TDomain::position_accessor_type m_aaPos;

		// Finite Volume Element Geometry
		FVElementGeometry<TElem, dim>* m_geo;

		// amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		number m_upwind_amount;

		// User functions
		Pososity_fct m_Porosity;
		Viscosity_fct m_Viscosity;
		Density_fct m_Density;
		D_Density_fct m_D_Density;
		Mol_Diff_Tensor_fct m_Mol_Diff_Tensor;
		Permeability_Tensor_fct m_Permeability_Tensor;
		Gravity_fct m_Gravity;

		DataClassExportPossibility<MathVector<dim>, MathVector<dim>, TAlgebra>& m_Darcy_Velocity_export;
		// constant values, read in once
		number m_porosity;

	private:
		// local constants for readability (local function 0 == _C_, local function 1 == _P_)
		static const size_t _C_ = 0;
		static const size_t _P_ = 1;

		// local access to local solution vector splitted for each component
		// TODO: Implement access to u, J, d as member functions, avoid defines !!!

	private:
		void export1(std::vector<MathVector<dim> >& val, std::vector<std::vector<MathVector<dim> > >& deriv, const std::vector<MathVector<ref_dim> >& pos, const local_vector_type& u, bool compute_derivatives)
		{
		#define u(fct, i)    ( (u)[ref_elem_type::num_corners*(fct) + (i)])

			const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
					LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);
			const size_t num_co = ref_elem_type::num_corners;


			number c;
			MathVector<TDomain::dim> grad_p, grad_p_local;
			MathMatrix<TDomain::dim, ref_elem_type::dim> JTInv;
			number shape;
			MathVector<TDomain::dim> grad_local;

			VecSet(grad_p, 0.0);
			m_mapping.update(m_corners);
			for(size_t i = 0; i < val.size(); ++i)
			{
				// compute c and grad_p
				VecSet(grad_p_local, 0.0);
				c = 0.0;

				for(size_t co = 0; co < num_co; ++co)
				{
					if(TrialSpace.evaluate(co, pos[i], shape) == false) UG_ASSERT(0, "");
					c += u(_C_, co) * shape;

					if(TrialSpace.evaluate_grad(co, pos[i], grad_local) == false)  UG_ASSERT(0, "");
					if(m_mapping.jacobian_transposed_inverse(pos[i], JTInv) == false) UG_ASSERT(0, "");
					VecScaleAppend(grad_p_local, u(_P_,co), grad_local);

					MatVecMult(grad_p, JTInv, grad_p_local);
				}

				compute_ip_Darcy_velocity(val[i], c, grad_p);
			}

			if(compute_derivatives == true)
			{

			}
		#undef u
		}



	public:
		/* HELP FUNCTIONS, TODO: Make a nicer solution */
		void compute_ip_Darcy_velocity(MathVector<dim>& Darcy_vel, number c_ip, const MathVector<dim>& grad_p_ip)
		{
			number s, viscosity_ip;
			MathVector<dim> vel;
			MathMatrix<dim, dim> K;

			m_Density(s, c_ip);
			m_Viscosity(viscosity_ip, c_ip);
			m_Gravity(vel);
			m_Permeability_Tensor(K);
			VecScale(vel, vel, s);
			VecSubtract(vel, vel, grad_p_ip);
			MatVecMult(Darcy_vel, K, vel);
			VecScale(Darcy_vel, Darcy_vel, 1./viscosity_ip);
		};


		void compute_D_ip_Darcy_velocity(	const SubControlVolumeFace<TElem, dim>& scvf,
											MathVector<dim>& Darcy_vel, MathVector<dim> D_Darcy_vel_c[], MathVector<dim> D_Darcy_vel_p[],
											number c_ip, const MathVector<dim>& grad_p_ip)
		{
			const int num_co = ref_elem_type::num_corners;
			number s, mu_ip;
			MathVector<dim> vel, gravity;
			MathVector<dim> D_vel_c[num_co], D_vel_p[num_co];
			MathMatrix<dim, dim> K;
			const SD_Values<TElem, dim>& sdv = scvf.sdv();

			m_Density(s, c_ip);
			m_Gravity(gravity);
			m_Viscosity(mu_ip, c_ip);
			m_Permeability_Tensor(K);
			VecScale(vel, gravity, s);
			VecSubtract(vel, vel, grad_p_ip);
			MatVecMult(Darcy_vel, K, vel);
			VecScale(Darcy_vel, Darcy_vel, 1./mu_ip);

			m_D_Density(s, c_ip);
			for(int co = 0; co < num_co; ++co)
			{
				VecScale(D_vel_c[co], gravity, s*sdv.shape(co));
				VecScale(D_vel_p[co], sdv.grad_global(co), -1);
				MatVecMult(D_Darcy_vel_c[co], K, D_vel_c[co]);
				MatVecMult(D_Darcy_vel_p[co], K, D_vel_p[co]);

				VecScale(D_Darcy_vel_c[co],D_Darcy_vel_c[co],1./mu_ip);
				VecScale(D_Darcy_vel_p[co],D_Darcy_vel_p[co],1./mu_ip);
			}

			// D_Viscosity == 0 !!!!
		};


};


template <typename TDomain, typename TAlgebra>
class DensityDrivenFlowPlugIn : public IPlugInElementDiscretization<TAlgebra>, public DataExportingClass<MathVector<TDomain::dim>, MathVector<TDomain::dim>,TAlgebra>{

	public:
		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = domain_type::dim;

		// position type
		typedef typename domain_type::position_type position_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef typename algebra_type::matrix_type::local_matrix_type local_matrix_type;

		// local vector tyoe
		typedef typename algebra_type::vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename algebra_type::vector_type::local_index_type local_index_type;

	protected:
		typedef void (*Pososity_fct)(number&);
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);
		typedef void (*Mol_Diff_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Permeability_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Gravity_fct)(MathVector<dim>&);

	public:
		DensityDrivenFlowPlugIn(size_t c_fct, size_t p_fct, TDomain& domain, number upwind_amount,
				Pososity_fct Porosity, Viscosity_fct Viscosity, Density_fct Density, D_Density_fct D_Density,
				Mol_Diff_Tensor_fct Mol_Diff, Permeability_Tensor_fct Permeability_Tensor, Gravity_fct Gravity) :
			m_c_fct(c_fct), m_p_fct(p_fct),
			m_Darcy_velocity_export("Darcy velocity"),
			m_ImpTriangle(domain, upwind_amount, Porosity, Viscosity, Density, D_Density, Mol_Diff, Permeability_Tensor, Gravity, m_Darcy_velocity_export),
			m_ImpQuadrilateral(domain, upwind_amount, Porosity, Viscosity, Density, D_Density, Mol_Diff, Permeability_Tensor, Gravity, m_Darcy_velocity_export),
			m_ImpTetrahedron(domain, upwind_amount, Porosity, Viscosity, Density, D_Density, Mol_Diff, Permeability_Tensor, Gravity, m_Darcy_velocity_export)
			{};

	public:
		/* GENERAL INFORMATIONS */
		// number of fundamental functions required for this assembling
		inline size_t num_fct(){return 2;}

		// local shape function set required for the 'i'-th fundamental function
		inline LocalShapeFunctionSetID local_shape_function_set(size_t i)
		{
			UG_ASSERT(i < num_fct(), "Accessing fundamental function, that is not contained in this assembling.");
			return LSFS_LAGRANGEP1;
		}

		// global number of fundamental function of local fundamental function
		inline size_t fct(size_t i)
		{
			UG_ASSERT(i < 2, "D3F has only two component.");
			switch(i)
			{
			case 0: return m_c_fct;
			case 1: return m_p_fct;
			}
			return (size_t)-1;
		}

	protected:
		// number of fundamental function, where this assembling works
		size_t m_c_fct;
		size_t m_p_fct;

	public:
		bool register_exports(DataContainer& Cont)
		{
			if(Cont.register_item(m_Darcy_velocity_export) != true)
			{
				UG_ASSERT(0, "Must work.");
				return false;
			}
			return true;
		}

		bool unregister_exports(DataContainer& Cont)
		{
			if(Cont.register_item(m_Darcy_velocity_export) != true)
			{
				UG_ASSERT(0, "Must work.");
				return false;
			}
			return true;
		}

		bool register_imports(DataContainer& Cont)
		{
			return true;
		}

		bool unregister_imports(DataContainer& Cont)
		{
			return true;
		}

	protected:
		DataClassExportPossibility<MathVector<dim>, MathVector<dim>,TAlgebra> m_Darcy_velocity_export;

	/* ELEMENT WISE ASSEMBLNG */
	public:
		// support assembling on triangles
		inline size_t num_sh(Triangle* elem)
		{ return m_ImpTriangle.num_sh();};

		inline size_t num_sh(Triangle* elem, size_t loc_fct)
		{ return m_ImpTriangle.num_sh(loc_fct);};

		inline IPlugInReturn prepare_element_loop(Triangle* elem)
		{ return m_ImpTriangle.prepare_element_loop(); };

		inline IPlugInReturn prepare_element(Triangle* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpTriangle.prepare_element(elem, u, glob_ind); };

		inline IPlugInReturn assemble_element_JA(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_JA(J, u, time); };

		inline IPlugInReturn assemble_element_JM(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_JM(J, u, time); };

		inline IPlugInReturn assemble_element_A(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_A(d, u, time); };

		inline IPlugInReturn assemble_element_M(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_M(d, u, time); };

		inline IPlugInReturn assemble_element_f(Triangle* elem, local_vector_type& d, number time=0.0)
		{ return m_ImpTriangle.assemble_element_f(d, time); };

		inline IPlugInReturn finish_element_loop(Triangle* elem)
		{ return m_ImpTriangle.finish_element_loop(); };

	protected:
		DensityDrivenFlow<domain_type, algebra_type, Triangle> m_ImpTriangle;


	public:
		// support assembling on triangles
		inline size_t num_sh(Quadrilateral* elem)
		{ return m_ImpQuadrilateral.num_sh();};

		inline size_t num_sh(Quadrilateral* elem, size_t loc_fct)
		{ return m_ImpQuadrilateral.num_sh(loc_fct);};

		inline IPlugInReturn prepare_element_loop(Quadrilateral* elem)
		{ return m_ImpQuadrilateral.prepare_element_loop(); };

		inline IPlugInReturn prepare_element(Quadrilateral* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpQuadrilateral.prepare_element(elem, u, glob_ind); };

		inline IPlugInReturn prepare_element_loop(Quadrilateral* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpQuadrilateral.prepare_element_loop(u, glob_ind); };

		inline IPlugInReturn assemble_element_JA(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_JA(J, u, time); };

		inline IPlugInReturn assemble_element_JM(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_JM(J, u, time); };

		inline IPlugInReturn assemble_element_A(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_A(d, u, time); };

		inline IPlugInReturn assemble_element_M(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_M(d, u, time); };

		inline IPlugInReturn assemble_element_f(Quadrilateral* elem, local_vector_type& d, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_f(d, time); };

		inline IPlugInReturn finish_element_loop(Quadrilateral* elem)
		{ return m_ImpQuadrilateral.finish_element_loop(); };

	protected:
		DensityDrivenFlow<domain_type, algebra_type, Quadrilateral> m_ImpQuadrilateral;

	public:
		// support assembling on triangles
		inline size_t num_sh(Tetrahedron* elem)
		{ return m_ImpTetrahedron.num_sh();};

		inline size_t num_sh(Tetrahedron* elem, size_t loc_fct)
		{ return m_ImpTetrahedron.num_sh(loc_fct);};

		inline IPlugInReturn prepare_element_loop(Tetrahedron* elem)
		{ return m_ImpTetrahedron.prepare_element_loop(); };

		inline IPlugInReturn prepare_element(Tetrahedron* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpTetrahedron.prepare_element(elem, u, glob_ind); };

		inline IPlugInReturn prepare_element_loop(Tetrahedron* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpTetrahedron.prepare_element_loop(u, glob_ind); };

		inline IPlugInReturn assemble_element_JA(Tetrahedron* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpTetrahedron.assemble_element_JA(J, u, time); };

		inline IPlugInReturn assemble_element_JM(Tetrahedron* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpTetrahedron.assemble_element_JM(J, u, time); };

		inline IPlugInReturn assemble_element_A(Tetrahedron* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpTetrahedron.assemble_element_A(d, u, time); };

		inline IPlugInReturn assemble_element_M(Tetrahedron* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpTetrahedron.assemble_element_M(d, u, time); };

		inline IPlugInReturn assemble_element_f(Tetrahedron* elem, local_vector_type& d, number time=0.0)
		{ return m_ImpTetrahedron.assemble_element_f(d, time); };

		inline IPlugInReturn finish_element_loop(Tetrahedron* elem)
		{ return m_ImpTetrahedron.finish_element_loop(); };

	protected:
		DensityDrivenFlow<domain_type, algebra_type, Tetrahedron> m_ImpTetrahedron;

};

}

#include "density_driven_flow_assemble_impl.h"

#endif /*__H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW_ASSEMBLE__*/
