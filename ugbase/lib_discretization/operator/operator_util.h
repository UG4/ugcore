/*
 * operator_util.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__
#define __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__

#include "lib_discretization/time_discretization/time_discretization_interface.h"
#include "lib_discretization/io/vtkoutput.h"
#include "common/profiler/profiler.h"

namespace ug{

template <typename TGridFunction>
bool ApplyLinearSolver(	ILinearOperator<typename TGridFunction::vector_type, typename TGridFunction::vector_type>& A,
						TGridFunction& u, TGridFunction& b,
						ILinearOperatorInverse<typename TGridFunction::vector_type, typename TGridFunction::vector_type>& solver)
{
// step 1: Prepare: Assemble matrix
	PROFILE_BEGIN(assembleLinearMatrix);
	if(!A.init())
		{UG_LOG("ApplyLinearSolver: Cannot init Operator.\n"); return false;}
	PROFILE_END();

// step 2: Init Linear Inverse Operator
	PROFILE_BEGIN(initLinearSolver);
	if(!solver.init(A))
		{UG_LOG("ApplyLinearSolver: Cannot init Inverse operator.\n"); return false;}
	PROFILE_END();

// step 4: Apply Operator
	PROFILE_BEGIN(applyLinearSolver);
	if(!solver.apply_return_defect(u,b))
		{UG_LOG("ApplyLinearSolver: Cannot apply Inverse operator.\n"); return false;}
	PROFILE_END();

	return true;
}

template <typename TGridFunction>
bool
PerformTimeStep(IOperatorInverse<typename TGridFunction::vector_type, typename TGridFunction::vector_type>& newton,
				TGridFunction& u,
				ITimeDiscretization<typename TGridFunction::dof_distribution_type, typename TGridFunction::algebra_type>& timestep,
				size_t timesteps, size_t step,
				number time, number dt,
				VTKOutput<TGridFunction>& out, const char* outName)
{
//	declare Vector type
	typedef typename TGridFunction::algebra_type::vector_type vector_type;

//  create deques for old solutions and old timesteps
	std::deque<vector_type*> u_old;
	std::deque<number> time_old;
	u_old.resize(timestep.num_prev_steps());
	time_old.resize(timestep.num_prev_steps());

//  set start time
	time_old[0] = time;

//  get help vectors
	TGridFunction* uOldFunc = &u.clone();
	u_old[0] = uOldFunc;

//  end timestep
	const size_t end_timestep = timesteps + step -1;

//  loop over time steps
	for(; step <= end_timestep; ++step)
	{
		UG_LOG("++++++ TIMESTEP " << step << " BEGIN ++++++\n");

		//	prepare time step
		timestep.prepare_step(u_old, time_old, dt);

		// prepare newton solver
		if(!newton.prepare(u))
			{UG_LOG("Cannot prepare Time step " << step << ". Aborting.\n"); return false;}

		//	execute newton solver
		if(!newton.apply(u))
			{UG_LOG("Time step " << step << " did not converge. Aborting.\n"); return false;}

		// update time
		time = time_old[0] + dt;
		time_old.pop_back();
		time_old.push_front(time);

		// update old solutions
		vector_type* current_u = u_old.back();
		u_old.pop_back();
		*current_u = u;
		u_old.push_front(current_u);

		// plot solution to file
		out.print(outName, u, step, time);
		UG_LOG("++++++ TIMESTEP " << step << "  END ++++++\n");
	}

	// free help vectors
	delete uOldFunc;

	return true;
}


} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__ */
