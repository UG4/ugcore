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
bool
PerformTimeStep(IOperatorInverse<typename TGridFunction::vector_type, typename TGridFunction::vector_type>& newton,
				TGridFunction& u,
				ITimeDiscretization<typename TGridFunction::dof_distribution_type::implementation_type, typename TGridFunction::algebra_type>& timestep,
				size_t timesteps, size_t step,
				number time, number dt,
				VTKOutput<TGridFunction>& out, const char* outName, bool bDoOutput)
{
//	declare Vector type
	typedef typename TGridFunction::algebra_type::vector_type vector_type;

//  create deques for old solutions and old timesteps
	PreviousSolutions<vector_type> prevSol;

//	clone current solution
	TGridFunction* uOldFunc = &u.clone();

//	store current time as old solution
	prevSol.push(*uOldFunc, time);

//  end timestep
	const size_t end_timestep = timesteps + step -1;

//  loop over time steps
	for(; step <= end_timestep; ++step)
	{
		UG_LOG("++++++ TIMESTEP " << step << " BEGIN ++++++\n");

		//	prepare time step
		timestep.prepare_step(prevSol, dt);

		// prepare newton solver
		if(!newton.prepare(u))
			{UG_LOG("Cannot prepare Time step " << step << ". Aborting.\n"); return false;}

		//	execute newton solver
		if(!newton.apply(u))
			{UG_LOG("Time step " << step << " did not converge. Aborting.\n"); return false;}

		// update time
		time = prevSol.time(0) + dt;

		// plot solution to file
		if(bDoOutput) out.print(outName, u, step, time);

		// get oldest solution
		vector_type& oldestSol = prevSol.oldest_solution();

		// copy values into oldest solution
		oldestSol = u;

		// push oldest solutions with new values to front
		prevSol.push_discard_oldest(oldestSol, time);

		UG_LOG("++++++ TIMESTEP " << step << "  END ++++++\n");
	}

	// free help vectors
	delete uOldFunc;

	return true;
}


} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__OPERATOR__OPERATOR_UTIL__ */
