/*
 * problem.h
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION__

#include <vector>
#include <string>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

#include "lib_discretization/domain_discretization/domain_discretization_interface.h"
#include "lib_discretization/domain_discretization/disc_coupling/element_data.h"

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

	public:
		// support assembling on tetrahedron
		virtual size_t num_sh(Tetrahedron* elem) = 0;

		virtual size_t num_sh(Tetrahedron* elem, size_t fct) = 0;

		virtual IPlugInReturn prepare_element_loop(Tetrahedron* elem) = 0;

		virtual IPlugInReturn prepare_element(Tetrahedron* elem, const local_vector_type& u, const local_index_type& glob_ind) = 0;

		virtual IPlugInReturn assemble_element_JA(Tetrahedron* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_JM(Tetrahedron* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_A(Tetrahedron* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_M(Tetrahedron* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) = 0;

		virtual IPlugInReturn assemble_element_f(Tetrahedron* elem, local_vector_type& d, number time=0.0) = 0;

		virtual IPlugInReturn finish_element_loop(Tetrahedron* elem) = 0;


		virtual ~ISystemDomainDiscretization()
		{};

};

// implement all virtual functions by Elem Disc
template <typename TDomain, typename TAlgebra, template<typename TDomain, typename TAlgebra > class TElemDisc>
class SystemDomainDiscretizationPlugIn : public ISystemDomainDiscretization<TDomain, TAlgebra>{

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

	public:
		SystemDomainDiscretizationPlugIn(TElemDisc<domain_type, algebra_type>& elemDisc) :
			  m_elemDisc(elemDisc)
			{};

	/* GENERAL INFORMATIONS */
	public:
		// number of fundamental functions required for this assembling
		virtual size_t num_fct(){return m_elemDisc.num_fct();}

		// local shape function set required for the 'i'-th fundamental function
		virtual LocalShapeFunctionSetID local_shape_function_set(size_t i){return m_elemDisc.local_shape_function_set(i);}

		// get global fundamental function number from local on
		virtual size_t fct(size_t i) {return m_elemDisc.fct(i);}

	public:
		virtual bool register_exports(DataContainer& Cont) {return m_elemDisc.register_exports(Cont);}

		virtual bool unregister_exports(DataContainer& Cont) {return m_elemDisc.unregister_exports(Cont);}

		virtual bool register_imports(DataContainer& Cont) {return m_elemDisc.register_imports(Cont);}

		virtual bool unregister_imports(DataContainer& Cont) {return m_elemDisc.unregister_imports(Cont);}

	/* ELEMENT WISE ASSEMBLNG */
	public:
		// support assembling on triangles
		virtual size_t num_sh(Triangle* elem){return m_elemDisc.num_sh(elem);}

		virtual size_t num_sh(Triangle* elem, size_t fct){return m_elemDisc.num_sh(elem, fct);}

		virtual IPlugInReturn prepare_element_loop(Triangle* elem){return m_elemDisc.prepare_element_loop(elem);}

		virtual IPlugInReturn prepare_element(Triangle* elem, const local_vector_type& u, const local_index_type& glob_ind){return m_elemDisc.prepare_element(elem, u, glob_ind);}

		virtual IPlugInReturn assemble_element_JA(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_JA(elem, J, u, time);}

		virtual IPlugInReturn assemble_element_JM(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_JM(elem, J, u, time);}

		virtual IPlugInReturn assemble_element_A(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) {return m_elemDisc.assemble_element_A(elem, d, u, time);}

		virtual IPlugInReturn assemble_element_M(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_M(elem, d, u, time);}

		virtual IPlugInReturn assemble_element_f(Triangle* elem, local_vector_type& d, number time=0.0) {return m_elemDisc.assemble_element_f(elem, d, time);}

		virtual IPlugInReturn finish_element_loop(Triangle* elem){return m_elemDisc.finish_element_loop(elem);}

	public:
		// support assembling on quadrilateral
		virtual size_t num_sh(Quadrilateral* elem){return m_elemDisc.num_sh(elem);}

		virtual size_t num_sh(Quadrilateral* elem, size_t fct){return m_elemDisc.num_sh(elem, fct);}

		virtual IPlugInReturn prepare_element_loop(Quadrilateral* elem){return m_elemDisc.prepare_element_loop(elem);}

		virtual IPlugInReturn prepare_element(Quadrilateral* elem, const local_vector_type& u, const local_index_type& glob_ind){return m_elemDisc.prepare_element(elem, u, glob_ind);}

		virtual IPlugInReturn assemble_element_JA(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_JA(elem, J, u, time);}

		virtual IPlugInReturn assemble_element_JM(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_JM(elem, J, u, time);}

		virtual IPlugInReturn assemble_element_A(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) {return m_elemDisc.assemble_element_A(elem, d, u, time);}

		virtual IPlugInReturn assemble_element_M(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_M(elem, d, u, time);}

		virtual IPlugInReturn assemble_element_f(Quadrilateral* elem, local_vector_type& d, number time=0.0) {return m_elemDisc.assemble_element_f(elem, d, time);}

		virtual IPlugInReturn finish_element_loop(Quadrilateral* elem){return m_elemDisc.finish_element_loop(elem);}

	public:
		// support assembling on quadrilateral
		virtual size_t num_sh(Tetrahedron* elem){return m_elemDisc.num_sh(elem);}

		virtual size_t num_sh(Tetrahedron* elem, size_t fct){return m_elemDisc.num_sh(elem, fct);}

		virtual IPlugInReturn prepare_element_loop(Tetrahedron* elem){return m_elemDisc.prepare_element_loop(elem);}

		virtual IPlugInReturn prepare_element(Tetrahedron* elem, const local_vector_type& u, const local_index_type& glob_ind){return m_elemDisc.prepare_element(elem, u, glob_ind);}

		virtual IPlugInReturn assemble_element_JA(Tetrahedron* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_JA(elem, J, u, time);}

		virtual IPlugInReturn assemble_element_JM(Tetrahedron* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_JM(elem, J, u, time);}

		virtual IPlugInReturn assemble_element_A(Tetrahedron* elem, local_vector_type& d, const local_vector_type& u, number time=0.0) {return m_elemDisc.assemble_element_A(elem, d, u, time);}

		virtual IPlugInReturn assemble_element_M(Tetrahedron* elem, local_vector_type& d, const local_vector_type& u, number time=0.0){return m_elemDisc.assemble_element_M(elem, d, u, time);}

		virtual IPlugInReturn assemble_element_f(Tetrahedron* elem, local_vector_type& d, number time=0.0) {return m_elemDisc.assemble_element_f(elem, d, time);}

		virtual IPlugInReturn finish_element_loop(Tetrahedron* elem){return m_elemDisc.finish_element_loop(elem);}

	protected:
		TElemDisc<domain_type, algebra_type>& m_elemDisc;

};


template <typename TDiscreteFunction>
class CoupledSystemDomainDiscretization : public IDomainDiscretization<TDiscreteFunction> {
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
		CoupledSystemDomainDiscretization(std::string name, DirichletBNDValues<discrete_function_type>& bndDisc, domain_type& domain)
		: m_name(name),
		  m_dirichletDisc(bndDisc), m_domain(domain)
		{
			m_elem_subset.resize(4);
			for(std::size_t i = 0; i < 3; ++i)
				m_elem_subset[i].clear();
		};

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
		bool link(int exp, int imp){return m_ElemDataContainer.link(exp, imp);};

	protected:
		// list of systems coupled by this discretization
		std::vector<system_type*> m_systems;

		// Containes the Data, that will be linked together
		DataContainer m_ElemDataContainer;

		// name of discretiztaion
		std::string m_name;

	protected:
		// TODO: These implementations of the interface are nearly the same as those for PlugInDomainDiscretization
		// Maybe we could find a common framework to avoid repetitive code

		// Assemble routines for time independent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u);
		IAssembleReturn assemble_solution(discrete_function_type& u);

		// Assemble routines for time dependent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_solution(discrete_function_type& u, number time);

	public:
		size_t num_elem_subsets(size_t d) {return m_elem_subset[d].size();}
		int elem_subset(size_t d, size_t i) {return m_elem_subset[d][i];}
		bool add_elem_assemble_subset(size_t d, int s)
		{
			std::vector<int>::iterator iter = find(m_elem_subset[d].begin(), m_elem_subset[d].end(), s);
			if(iter == m_elem_subset[d].end())
				m_elem_subset[d].push_back(s);
			return true;
		}

		virtual size_t num_fct() const
		{
			size_t num = 0;
			for(size_t sys = 0; sys < m_systems.size(); ++sys)
			{
				num += m_systems[sys]->num_fct();
			}
			return num;
		}

	protected:
		std::vector<std::vector<int> > m_elem_subset; // 3 = max dimensions

	protected:
		DirichletBNDValues<discrete_function_type>& m_dirichletDisc;

		domain_type& m_domain;

		virtual bool is_dirichlet(int si, size_t fct) {	return m_dirichletDisc.is_dirichlet(si, fct);}


	protected:
		template <typename TElem>
		bool assemble_jacobian_defect	(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											matrix_type& J, vector_type& d, const discrete_function_type& u,
											number time, number s_m, number s_a);
		template <typename TElem>
		bool assemble_jacobian			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											matrix_type& J, const discrete_function_type& u,
											number time, number s_m, number s_a);
		template <typename TElem>
		bool assemble_defect			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											vector_type& d, const discrete_function_type& u,
											number time, number s_m, number s_a);

		template <typename TElem>
		bool assemble_linear			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											matrix_type& mat, vector_type& rhs, const discrete_function_type& u);

		bool clear_dirichlet_jacobian_defect(	geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, size_t sys, size_t loc_fct, int si, matrix_type& J, vector_type& d, const discrete_function_type& u, number time = 0.0);
		bool clear_dirichlet_jacobian(			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, size_t sys, size_t loc_fct, int si, matrix_type& J, const discrete_function_type& u, number time = 0.0);
		bool clear_dirichlet_defect(			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, size_t sys, size_t loc_fct, int si, vector_type& d, const discrete_function_type& u, number time = 0.0);
		bool set_dirichlet_solution( 			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, size_t sys, size_t loc_fct, int si, vector_type& x, const discrete_function_type& u, number time = 0.0);
		bool set_dirichlet_linear(				geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, size_t sys, size_t loc_fct, int si, matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time = 0.0);
};

} // end namespace ug

#include "coupled_system_domain_discretization_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION__ */
