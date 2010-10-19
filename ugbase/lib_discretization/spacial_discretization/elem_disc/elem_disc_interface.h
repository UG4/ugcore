/*
 * elem_disc_interface.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__

// extern headers
#include <vector>

// intern headers
#include "lib_discretization/assemble.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{

enum IElemDiscNeed {
	IEDN_NONE = 0,
	IEDN_DEFECT = 1 << 0,
	IEDN_JACOBIAN = 1 << 1,
	IEDN_LINEAR = 1 << 2,
	IEDN_STIFFNESS = 1 << 3,
	IEDN_MASS = 1 << 4
};

template <typename TAlgebra>
class IElemDisc{
	protected:
		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef LocalMatrix<typename TAlgebra::matrix_type::entry_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::entry_type> local_vector_type;

		// local index type
		typedef LocalIndices local_index_type;

	public:
		IElemDisc() : m_pPattern(NULL) {}

		void set_pattern(const FunctionPattern& pattern)
		{
			m_pPattern = &pattern;
			m_FunctionGroup.set_function_pattern(pattern);
			m_SubsetGroup.set_subset_handler(*pattern.get_subset_handler());
		}

		virtual bool set_functions(const char* functions);
		virtual bool set_functions(const FunctionGroup& funcGroup);

		virtual bool set_subsets(const char* subsets);
		virtual bool set_subsets(const SubsetGroup& subsetGroup);

		const FunctionGroup& get_function_group() const {return m_FunctionGroup;}
		const SubsetGroup& get_subset_group() const {return m_SubsetGroup;}

	protected:

		const FunctionPattern* m_pPattern;
		FunctionGroup m_FunctionGroup;
		SubsetGroup m_SubsetGroup;

	public:
		// number of functions this discretization handles
		virtual size_t num_fct() = 0;

		// shape function set of the functions handled by this discretization
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct) = 0;

		// returns, if hanging nodes should be considered in elem disc as well (default is false)
		virtual bool use_hanging() const {return false;}

		// sets the geometric object type
		// ATTENTION: type must be set, before other public functions can be called
		bool set_geometric_object_type(int id, IElemDiscNeed need);

		// preparing and finishing of loop
		bool prepare_element_loop()						{return (this->*(m_vPrepareElementLoopFunc[m_id]))();}
		bool prepare_element(GeometricObject* obj, const local_vector_type& u, const local_index_type& glob_ind)
														{return (this->*(m_vPrepareElementFunc[m_id]))(obj, u, glob_ind);}
		bool finish_element_loop()						{return (this->*(m_vFinishElementLoopFunc[m_id]))();}

		// assembling of Jacobian
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0) 	{return (this->*(m_vAssembleJAFunc[m_id]))(J, u, time);}
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0)		{return (this->*(m_vAssembleJMFunc[m_id]))(J, u, time);}

		// assembling of defect
		bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0)		{return (this->*(m_vAssembleAFunc[m_id]))(d, u, time);}
		bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0)		{return (this->*(m_vAssembleMFunc[m_id]))(d, u, time);}

		// assembling of right hand side for linear case
		bool assemble_f(local_vector_type& d, number time=0.0) 									{return (this->*(m_vAssembleFFunc[m_id]))(d, time);}

		// virtual destructor
		virtual ~IElemDisc() {}

	protected:
		// register the functions
		template <typename TAssFunc> void register_prepare_element_loop_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_prepare_element_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_finish_element_loop_function(int id, TAssFunc func);

		template <typename TAssFunc> void register_assemble_JA_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_JM_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_A_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_M_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_f_function(int id, TAssFunc func);

	protected:
		// checks if the needed functions are registered for the id type
		bool function_registered(int id, IElemDiscNeed need);

		// checks if the functions are present
		bool prepare_element_loop_function_registered(int id);
		bool prepare_element_function_registered(int id);
		bool finish_element_loop_function_registered(int id);

		// checks if the functions are present
		bool assemble_JA_function_registered(int id);
		bool assemble_JM_function_registered(int id);
		bool assemble_A_function_registered(int id);
		bool assemble_M_function_registered(int id);
		bool assemble_f_function_registered(int id);

	private:
		// types of loop function pointers
		typedef bool (IElemDisc<TAlgebra>::*PrepareElementLoopFunc)();
		typedef bool (IElemDisc<TAlgebra>::*PrepareElementFunc)(GeometricObject* obj, const local_vector_type& u, const local_index_type& glob_ind);
		typedef bool (IElemDisc<TAlgebra>::*FinishElementLoopFunc)();

		// types of Jacobian assemble functions
		typedef bool (IElemDisc<TAlgebra>::*AssembleJAFunc)(local_matrix_type& J, const local_vector_type& u, number time);
		typedef bool (IElemDisc<TAlgebra>::*AssembleJMFunc)(local_matrix_type& J, const local_vector_type& u, number time);

		// types of Defect assemble functions
		typedef bool (IElemDisc<TAlgebra>::*AssembleAFunc)(local_vector_type& d, const local_vector_type& u, number time);
		typedef bool (IElemDisc<TAlgebra>::*AssembleMFunc)(local_vector_type& d, const local_vector_type& u, number time);

		// types of right hand side assemble functions
		typedef bool (IElemDisc<TAlgebra>::*AssembleFFunc)(local_vector_type& d, number time);

	private:
		// loop function pointers
		std::vector<PrepareElementLoopFunc> m_vPrepareElementLoopFunc;
		std::vector<PrepareElementFunc> 	m_vPrepareElementFunc;
		std::vector<FinishElementLoopFunc> 	m_vFinishElementLoopFunc;

		// Jacobian function pointers
		std::vector<AssembleJAFunc> 	m_vAssembleJAFunc;
		std::vector<AssembleJMFunc> 	m_vAssembleJMFunc;

		// Defect function pointers
		std::vector<AssembleAFunc> 	m_vAssembleAFunc;
		std::vector<AssembleMFunc> 	m_vAssembleMFunc;

		// Rhs function pointers
		std::vector<AssembleFFunc> 	m_vAssembleFFunc;

	protected:
		// current Geometric Object
		int m_id;

};

} // end namespace ug

#include "elem_disc_interface_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__ */
