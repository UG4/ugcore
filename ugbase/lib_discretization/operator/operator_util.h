/*
 * operator.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__
#define __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__

#include "lib_discretization/time_discretization/time_discretization_interface.h"
#include "lib_discretization/io/vtkoutput.h"

namespace ug{

template <typename TGridFunction>
bool ApplyLinearSolver(	ILinearOperator<TGridFunction, TGridFunction>& A,
						TGridFunction& u, TGridFunction& b,
						ILinearizedOperatorInverse<TGridFunction, TGridFunction>& solver)
{
	// step 1: Prepare operator, compute b and set dirichlet values in u
	PROFILE_BEGIN(CreateOperator);
	if(!A.prepare(u,b))
		{UG_LOG("ApplyLinearSolver: Cannot init Operator.\n"); return false;}
	PROFILE_END();

	// step 2: Init Linear Inverse Operator
	if(!solver.init(A))
		{UG_LOG("ApplyLinearSolver: Cannot init Inverse operator.\n"); return false;}

	// step 3: Prepare Linear Inverse Operator
	PROFILE_BEGIN(CreateInverseOperators);
	if(!solver.prepare(u, b, u))
		{UG_LOG("ApplyLinearSolver: Cannot prepare Inverse operator.\n"); return false;}
	PROFILE_END();

	// step 4: Apply Operator
	PROFILE_BEGIN(InvertOperator);
	if(!solver.apply(b,u))
		{UG_LOG("ApplyLinearSolver: Cannot apply Inverse operator.\n"); return false;}
	PROFILE_END();

	return true;
}

template <typename TGridFunction, typename TLevelFunction>
bool
PerformTimeStep(IOperatorInverse<TGridFunction, TGridFunction>& newton,
				TGridFunction& u,
				ITimeDiscretization<TGridFunction>& surface_timestep,
				ITimeDiscretization<TLevelFunction>& level_timestep,
				size_t timesteps, number dt,
				VTKOutput<TGridFunction>& out, const char* outName)
{
//  create deques for old solutions and old timesteps
	std::deque<TGridFunction*> u_old;
	std::deque<number> time_old;
	u_old.resize(surface_timestep.num_prev_steps());
	time_old.resize(surface_timestep.num_prev_steps());

//  set start time
	time_old[0] = 0.0;

//  get help vectors
	u_old[0] = &u.clone();

//  loop over time steps
	for(size_t step = 1; step <= timesteps; ++step)
	{
		UG_LOG("++++++ TIMESTEP " << step << " BEGIN ++++++\n");

		//	prepare time step
		surface_timestep.prepare_step(u_old, time_old, dt);
		// TODO: This is wrong, should be u_old_level
		level_timestep.prepare_step(u_old, time_old, dt);

		// preapre newton solver
		if(!newton.prepare(u))
			{UG_LOG("Cannot prepare Time step " << step << ". Aborting.\n"); return false;}

		//	execute newton solver
		if(!newton.apply(u))
			{UG_LOG("Time step " << step << " did not converge. Aborting.\n"); return true;}

		// update time
		number current_time = time_old[0] + dt;
		time_old.pop_back();
		time_old.push_front(current_time);

		// update old solutions
		TGridFunction* current_u = u_old.back();
		u_old.pop_back();
		*current_u = u;
		u_old.push_front(current_u);

		// plot solution to file
		out.print(outName, u, step, current_time);
		UG_LOG("++++++ TIMESTEP " << step << " END ++++++\n");
	}

	// free help vectors
	delete u_old[0];

	return true;
}


} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__ */
