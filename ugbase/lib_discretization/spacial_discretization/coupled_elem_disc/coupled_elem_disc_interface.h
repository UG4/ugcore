/*
 * coupled_elem_disc_interface.h
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_INTERFACE__

#include <vector>
#include <string>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

#include "lib_discretization/spacial_discretization/domain_discretization_interface.h"
#include "lib_discretization/spacial_discretization/disc_coupling/element_data.h"

namespace ug {

// base class
template <typename TDomain, typename TAlgebra>
class ISystemDomainDiscretization {

	public:
		// domain type
		typedef TDomain domain_type;

		// position type
		typedef typename domain_type::position_type position_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef typename algebra_type::matrix_type::local_matrix_type local_matrix_type;

		// local vector type
		typedef typename algebra_type::vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename algebra_type::vector_type::local_index_type local_index_type;

	/* GENERAL INFORMATIONS */
	public:
		// number of fundamental functions required for this assembling
		virtual size_t num_fct() = 0;

		// local shape function set required for the 'i'-th fundamental function
		virtual LocalShapeFunctionSetID local_shape_function_set(size_t i) = 0;

		// get global fundamental function number from local on
		virtual size_t fct(size_t i) = 0;

	public:
		virtual size_t num_imports() = 0;

		virtual DataImportItem* import(size_t i) = 0;

		virtual bool register_exports(DataContainer& Cont) = 0;

		virtual bool unregister_exports(DataContainer& Cont) = 0;

		virtual bool register_imports(DataContainer& Cont) = 0;

		virtual bool unregister_imports(DataContainer& Cont) = 0;

	/* ELEMENT WISE ASSEMBLNG */
	public:
		// support assembling on triangles
		virtual size_t num_sh(Triangle* elem) = 0;

		virtual size_t num_sh(Triangle* elem, size_t fct) = 0;

		virtual IPlugInReturn prepare_element_loop(Triangle* elem) = 0;

		virtual IPlugInReturn prepare_element(Triangle* elem, const local_vector_type& u, const local_index_type& glob_ind) = 0;

		virtual IPlugInReturn assemble_element_JA(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_JM(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_A(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_M(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_f(Triangle* elem, local_vector_type& d, number time=0.0) = 0;

		virtual IPlugInReturn finish_element_loop(Triangle* elem) = 0;

	public:
		// support assembling on quadrilateral
		virtual size_t num_sh(Quadrilateral* elem) = 0;

		virtual size_t num_sh(Quadrilateral* elem, size_t fct) = 0;

		virtual IPlugInReturn prepare_element_loop(Quadrilateral* elem) = 0;

		virtual IPlugInReturn prepare_element(Quadrilateral* elem, const local_vector_type& u, const local_index_type& glob_ind) = 0;

		virtual IPlugInReturn assemble_element_JA(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_JM(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_A(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_M(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_f(Quadrilateral* elem, local_vector_type& d, number time=0.0) = 0;

		virtual IPlugInReturn finish_element_loop(Quadrilateral* elem) = 0;

		virtual ~ISystemDomainDiscretization()
		{};

};


template <typename TDiscreteFunction>
class CoupledSystemDomainDiscretization : public IDimensionDomainDiscretization<TDiscreteFunction> {
	protected:
		// forward constants and typenames

		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		// type of domain
		typedef typename TDiscreteFunction::domain_type domain_type;

		// type of grid used
		typedef typename domain_type::grid_type grid_type;

		// type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

		// type of dof manager
		typedef typename TDiscreteFunction::dof_manager_type dof_manager_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of local matrix
		typedef typename matrix_type::local_matrix_type local_matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		// type of algebra vector
		typedef typename vector_type::local_vector_type local_vector_type;

		// type of multi_index used in algebra
		typedef typename matrix_type::index_type index_type;

		// type of local index container
		typedef typename matrix_type::local_index_type local_index_type;

		// Virtual base class for system domains
		typedef ISystemDomainDiscretization<domain_type, algebra_type> system_type;

	public:
		CoupledSystemDomainDiscretization(std::string name)
			: m_name(name)
		{};

		std::string name() const {return m_name;}

		inline void add_system_discretization(system_type& sys)
		{
			m_systems.push_back(&sys);
			sys.register_exports(m_ElemDataContainer);
			sys.register_imports(m_ElemDataContainer);
		};

		size_t num_system_discretizations() {return m_systems.size();}

		bool print_export_possibilities(){return m_ElemDataContainer.print_export_possibilities(DCI_LINKS);}
		bool print_exports(){return m_ElemDataContainer.print_exports(DCI_LINKS);}
		bool print_imports(){return m_ElemDataContainer.print_imports(DCI_LINKS);};
		bool print_linkage(){return m_ElemDataContainer.print_linkage();};
		bool link(size_t exp, size_t imp)
		{
			m_ElemDataContainer.create_export(exp);
			return m_ElemDataContainer.link(exp, imp);
		}

	protected:
		// list of systems coupled by this discretization
		std::vector<system_type*> m_systems;

		// Containes the Data, that will be linked together
		DataContainer m_ElemDataContainer;

		// name of discretiztaion
		std::string m_name;

	public:
		// TODO: NEEDED ?
		virtual size_t num_fct() const
		{
			size_t num = 0;
			for(size_t sys = 0; sys < m_systems.size(); ++sys)
			{
				num += m_systems[sys]->num_fct();
			}
			return num;
		}

};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_INTERFACE__ */
