/*
 * post_process_interface.h
 *
 *  Created on: 04.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__POST_PROCESS_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__POST_PROCESS_INTERFACE__

// extern headers
#include <vector>

// intern headers
#include "lib_discretization/common/local_algebra.h"

namespace ug{

enum IDiscPostProcessNeed {
	IEDN_NONE = 0,
	IEDN_DEFECT = 1 << 0,
	IEDN_JACOBIAN = 1 << 1,
	IEDN_LINEAR = 1 << 2
};

template <typename TAlgebra>
class IDiscPostProcess{
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
		// sets the geometric object type
		// ATTENTION: type must be set, before other public functions can be called
		bool set_geometric_object_type(int id, IDiscPostProcessNeed need);

		// preparing and finishing of loop
		bool prepare_element_loop()						{return (this->*(m_vPrepareElementLoopFunc[m_id]))();}
		bool prepare_element(GeometricObject* obj, const local_vector_type& u, const local_index_type& glob_ind)
														{return (this->*(m_vPrepareElementFunc[m_id]))(obj, u, glob_ind);}
		bool finish_element_loop()						{return (this->*(m_vFinishElementLoopFunc[m_id]))();}

		// post processing  of Jacobian
		bool post_process_J(local_matrix_type& J, const local_vector_type& u, number time=0.0) 	{return (this->*(m_vPostProcessJFunc[m_id]))(J, u, time);}

		// post processing  of defect
		bool post_process_d(local_vector_type& d, const local_vector_type& u, number time=0.0)	{return (this->*(m_vPostProcessDFunc[m_id]))(d, u, time);}

		// post processing of right hand side for linear case
		bool post_process_f(local_vector_type& d, number time=0.0) 								{return (this->*(m_vPostProcessFFunc[m_id]))(d, time);}

		// virtual destructor
		virtual ~IDiscPostProcess() {}

	protected:
		// register the functions
		template <typename TAssFunc> void register_prepare_element_loop_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_prepare_element_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_finish_element_loop_function(int id, TAssFunc func);

		template <typename TAssFunc> void register_post_process_J_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_post_process_d_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_post_process_f_function(int id, TAssFunc func);

	protected:
		// checks if the needed functions are registered for the id type
		bool function_registered(int id, IDiscPostProcessNeed need);

		// checks if the functions are present
		bool prepare_element_loop_function_registered(int id);
		bool prepare_element_function_registered(int id);
		bool finish_element_loop_function_registered(int id);

		// checks if the functions are present
		bool post_process_J_function_registered(int id);
		bool post_process_d_function_registered(int id);
		bool post_process_f_function_registered(int id);

	private:
		// types of loop function pointers
		typedef bool (IDiscPostProcess<TAlgebra>::*PrepareElementLoopFunc)();
		typedef bool (IDiscPostProcess<TAlgebra>::*PrepareElementFunc)(GeometricObject* obj, const local_vector_type& u, const local_index_type& glob_ind);
		typedef bool (IDiscPostProcess<TAlgebra>::*FinishElementLoopFunc)();

		// types of Jacobian assemble functions
		typedef bool (IDiscPostProcess<TAlgebra>::*PostProcessJFunc)(local_matrix_type& J, const local_vector_type& u, number time);

		// types of Defect assemble functions
		typedef bool (IDiscPostProcess<TAlgebra>::*PostProcessDFunc)(local_vector_type& d, const local_vector_type& u, number time);

		// types of right hand side assemble functions
		typedef bool (IDiscPostProcess<TAlgebra>::*PostProcessFFunc)(local_vector_type& d, number time);

	private:
		// loop function pointers
		std::vector<PrepareElementLoopFunc> m_vPrepareElementLoopFunc;
		std::vector<PrepareElementFunc> 	m_vPrepareElementFunc;
		std::vector<FinishElementLoopFunc> 	m_vFinishElementLoopFunc;

		// Jacobian function pointers
		std::vector<PostProcessJFunc> 	m_vPostProcessJFunc;

		// Defect function pointers
		std::vector<PostProcessDFunc> 	m_vPostProcessDFunc;

		// Rhs function pointers
		std::vector<PostProcessFFunc> 	m_vPostProcessFFunc;

	protected:
		// current Geometric Object
		int m_id;

};

} // end namespace ug

#include "post_process_interface_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__POST_PROCESS_INTERFACE__ */
