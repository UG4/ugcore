/*
 * coupled_system.h
 *
 *  Created on: 09.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_SYSTEM__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_SYSTEM__

#include "./disc_coupling/element_data.h"

namespace ug {

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type>
class CoupledSystem {
	protected:
		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		typedef TAlgebra algebra_type;
	public:
		CoupledSystem(std::string name)
			: m_name(name)
		{};

		std::string name() const {return m_name;}

		inline bool add_system_discretization(ICoupledElemDisc<algebra_type>& sys)
		{
			m_vSystem.push_back(&sys);
			if(!sys.register_exports(m_ElemDataContainer))
				{UG_LOG("CoupledSystem::add_system_discretization: Cannot register exports of system"); return false;}
			if(!sys.register_imports(m_ElemDataContainer))
				{UG_LOG("CoupledSystem::add_system_discretization: Cannot register imports of system"); return false;}
			if(!set_sys_ids())
				{UG_LOG("CoupledSystem::add_system_discretization: Cannot assign system ids"); return false;}
			return true;
		};

		size_t num_sys() {return m_vSystem.size();}
		ICoupledElemDisc<algebra_type>& sys(size_t i) {return *m_vSystem[i];}
		size_t sys_fct(size_t sys, size_t loc_fct)
		{
			size_t fct = 0;
			for(size_t s = 0; s < sys; ++s)
			{
				fct += m_vSystem[s]->num_fct();
			}
			return fct + loc_fct;
		}

		bool set_sys_ids()
		{
			for(size_t id = 0; id < num_sys(); ++id)
			{
				if(!m_vSystem[id]->set_sys_id(id)) return false;
			}
			return true;
		}


		DataContainer& get_elem_data_container() {return m_ElemDataContainer;}

		bool print_export_possibilities(){return m_ElemDataContainer.print_export_possibilities(DCI_LINKS);}
		bool print_exports(){return m_ElemDataContainer.print_exports(DCI_LINKS);}
		bool print_imports(){return m_ElemDataContainer.print_imports(DCI_LINKS);};
		bool print_linkage(){return m_ElemDataContainer.print_linkage();};
		bool link(size_t exp, size_t imp)
		{
			if(!m_ElemDataContainer.create_export(exp))
				{UG_LOG("CoupledSystem::link: Cannot create export"); return false;}

			if(!m_ElemDataContainer.link(exp, imp))
				{UG_LOG("CoupledSystem::link: Cannot link created export to import"); return false;}

			return true;
		}

		virtual ~CoupledSystem(){}
	protected:
		// list of systems coupled by this discretization
		std::vector<ICoupledElemDisc<algebra_type>*> m_vSystem;

		// Containes the Data, that will be linked together
		DataContainer m_ElemDataContainer;

		// name of discretiztaion
		std::string m_name;

	public:
		// TODO: NEEDED ?
		virtual size_t num_fct() const
		{
			size_t num = 0;
			for(size_t sys = 0; sys < m_vSystem.size(); ++sys)
			{
				num += m_vSystem[sys]->num_fct();
			}
			return num;
		}

};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_SYSTEM__ */
