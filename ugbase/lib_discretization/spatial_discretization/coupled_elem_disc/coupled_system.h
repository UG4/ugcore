/*
 * coupled_system.h
 *
 *  Created on: 09.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_SYSTEM__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_SYSTEM__

#include "./elem_data/element_data.h"
#include "lib_discretization/common/function_group.h"

namespace ug {

template <typename TAlgebra>
class CoupledSystem {
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	public:
		/// Constructor
		CoupledSystem(std::string name)
			: m_name(name)
		{};

		/// name of coupled system
		std::string name() const {return m_name;}


		///////////////////////////////////
		// Administration of systems
		///////////////////////////////////

		/// add a system
		inline bool add_system(ICoupledElemDisc<algebra_type>& sys)
		{
			m_vSystem.push_back(&sys);
			if(!sys.register_exports(m_ElemDataContainer))
				{UG_LOG("CoupledSystem::add_system_discretization: Cannot register exports of system"); return false;}
			if(!sys.register_imports(m_ElemDataContainer))
				{UG_LOG("CoupledSystem::add_system_discretization: Cannot register imports of system"); return false;}
			if(!set_sys_ids())
				{UG_LOG("CoupledSystem::add_system_discretization: Cannot assign system ids"); return false;}

			update_sub_function_maps();
			update_function_group();
			return true;
		};

		/// number of systems coupled
		size_t num_system() const {return m_vSystem.size();}

		/// access to a system
		ICoupledElemDisc<algebra_type>& system(size_t sys)
		{
			UG_ASSERT(sys < num_system(), "Invalid system id.");
			return *m_vSystem[sys];
		}

		/// get SubFunctionMap for a system
		const SubFunctionMap& sub_function_map(size_t sys) const
		{
			UG_ASSERT(sys < num_system(), "Invalid system id.");
			return m_vSubFunctionMap[sys];
		}

		/// set system id to all systems
		bool set_sys_ids()
		{
			for(size_t id = 0; id < num_system(); ++id)
			{
				if(!system(id).set_sys_id(id))
					{UG_LOG("CoupledSystem::set_sys_ids: Cannot set sys_id for system " << id << ".\n"); return false;}
			}
			return true;
		}

		const FunctionGroup& get_function_group() const {return m_FunctionGroup;}

		/// number of functions this coupled discretization requires
		size_t num_fct() const
		{
			/// TODO: Think about, if this value should be cached
			size_t num = 0;
			for(size_t sys = 0; sys < num_system(); ++sys)
			{
				num += m_vSystem[sys]->num_fct();
			}
			return num;
		}

		/// Destructor
		virtual ~CoupledSystem(){}

		///////////////////////////////////
		// Data Coupling
		///////////////////////////////////

		/// get data container
		DataContainer& get_elem_data_container() {return m_ElemDataContainer;}

		/// print export possibilities
		bool print_export_possibilities(){return m_ElemDataContainer.print_export_possibilities(DCI_LINKS);}

		/// print existing exports
		bool print_exports(){return m_ElemDataContainer.print_exports(DCI_LINKS);}

		/// print existing imports
		bool print_imports(){return m_ElemDataContainer.print_imports(DCI_LINKS);};

		/// print current linkage between import and export
		bool print_linkage(){return m_ElemDataContainer.print_linkage();};

		/// link export to import
		bool link(size_t exp, size_t imp)
		{
			if(!m_ElemDataContainer.link(exp, imp))
				{UG_LOG("CoupledSystem::link: Cannot link created export to import.\n"); return false;}

			return true;
		}

		/// create export
		bool create_export(size_t exp)
		{
			if(m_ElemDataContainer.create_export(exp) == NULL)
				{UG_LOG("CoupledSystem::link: Cannot create export.\n"); return false;}
			return true;
		}

	protected:
		void update_sub_function_maps()
		{
			m_vSubFunctionMap.resize(num_system());

			size_t fct = 0;
			for(size_t sys = 0; sys < num_system(); ++sys)
			{
				m_vSubFunctionMap[sys].clear();
				for(size_t loc_fct = 0; loc_fct < system(sys).num_fct(); ++loc_fct)
				{
					m_vSubFunctionMap[sys].select(fct);
					fct++;
				}
			}
		}

		void update_function_group()
		{
			m_FunctionGroup.clear();

			for(size_t sys = 0; sys < num_system(); ++sys)
			{
				const FunctionGroup& sysFunctionGroup = system(sys).get_function_group();

				// todo: checks if all underlying patterns are equal
				// todo: check if all subsets are equal

				for(size_t loc_fct = 0; loc_fct < sysFunctionGroup.num_fct(); ++loc_fct)
				{
					// todo: check that each function only once added (Does it make sense otherwise?)
					m_FunctionGroup.add_function(sysFunctionGroup[loc_fct]);
				}
			}
		}

	protected:
		// function group
		FunctionGroup m_FunctionGroup;

		// vector of systems coupled by this discretization
		std::vector<ICoupledElemDisc<algebra_type>*> m_vSystem;

		/// vector of subfunction maps: Mapping loc_fct_id of system to index in Function Group
		std::vector<SubFunctionMap> m_vSubFunctionMap;

		// Contains the Data, that will be linked together
		DataContainer m_ElemDataContainer;

		// name of discretization
		std::string m_name;
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_SYSTEM__ */
