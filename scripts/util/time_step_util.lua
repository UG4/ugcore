-- Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
-- Author: Andreas Vogel
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.

--[[!
\file time_step_util.lua
\defgroup scripts_util_timestep Timestep Utility
\ingroup scripts_util
\{
]]--

util = util or {}

--!
--! @param domainDisc
--! @param timeScheme theta, impleuler, expleuer, crank-nicolson, alexander, fracstep, or bdf
--! @param orderOrTheta theat is timeScheme=theta
function util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)

	local timeDisc = nil;

	if timeScheme:lower() == "theta" then 
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(orderOrTheta) -- 1.0 is implicit euler

	elseif timeScheme:lower() == "impleuler" then 
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(1.0) -- 1.0 is implicit euler

	elseif timeScheme:lower() == "expleuler" then 
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(0.0) -- 0.0 is explicit euler

	elseif timeScheme:lower() == "crank-nicolson" then 
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(0.5) -- 0.5 is crank-nicolson
	
	elseif timeScheme:lower() == "alexander" then 
	timeDisc = ThetaTimeStep(domainDisc, "Alexander")
	
	elseif timeScheme:lower() == "fracstep" then 
	timeDisc = ThetaTimeStep(domainDisc, "FracStep")
	
	elseif timeScheme:lower() == "bdf" then
	timeDisc = BDF(domainDisc)

	elseif timeScheme:lower() == "sdirk" then
	timeDisc = SDIRK(domainDisc, orderOrTheta)
	
	else
		print("Time scheme '"..timeScheme.."' not found. Supported Schemes:")
		print("Theta, ImplEuler, ExplEuler, Crank-Nicolson, Alexander, FracStep, BDF, SDIRK")
		print("Aborting")
		exit();
	end

	return timeDisc
end


function util.PrintUsageOfSolveTimeProblem()
	print("Usage of ")
	print("util.SolveNonlinearTimeProblem(u, domainDisc, solver, out, filename")
	print("                               timeScheme, orderOrTheta, startTime, endTime,")
	print("                               maxStepSize, minStepSize, reductionFactor)")
	print("util.SolveLinearTimeProblem(   u, domainDisc, solver, out, filename")
	print("                               timeScheme, orderOrTheta, startTime, endTime,")
	print("                               maxStepSize, minStepSize, reductionFactor)")
	print("Parameters:")
	print("------")
	print("u                -- [in] GridFunction with Startvalues, [out] Solution")
	print("domainDisc       -- Domain Discretization ")
	print("solver           -- Linear or Nonlinear Solver")
	print("out              -- a VTKOutput (pass nil for no output)")
	print("filename         -- filename for output")
	print("timeScheme       -- Name of time step scheme: ")
	print("                    Theta, ImplEuler, ExplEuler, Crank-Nicolson,")
	print("                    Alexander, FracStep, BDF")
	print("orderOrTheta	    -- theta param if 'Theta', order if 'BDF'")
	print("startTime        -- start time point")
	print("endTime          -- end time point")
	print("maxStepSize      -- maximal step sized used")	
	print("minStepSize      -- (optional) minimal step sized used")	
	print("reductionFactor  -- (optional) factor by which the step size is ")
	print("                    reduced, if the problem was not solved. ")
	print("                    Iterated until minStepSize is reached.")
	print("bFinishTimeStep  -- (optional) boolean if finish_timestep should be")
	print("                    called or not.")
	print("useCheckpointing	-- (optional) if true, use checkpointing.")
	print("postProcess		-- (optional) if passed, will be called after every time step.")
	print("startTSNo		-- (optional) time step number of the initial condition (normally 0).")
	print("endTSNo			-- (optional) if passed, stop after the time step with this number.")
end

--! Time stepping with a fixed step size. Returns number of time steps done and the last time.
--! @param u 				[in] GridFunction with Startvalues, [out] Solution"
--! @param domainDisc		Domain Discretization
--! @param newtonSolver		Nonlinear Solver
--! @param out				a VTKOutput (pass nil for no output)
--! @param filename			filename for output
--! @param timeScheme   	Name of time step scheme:
--!     	            	Theta, ImplEuler, ExplEuler, Crank-Nicolson
--! 	                	Alexander, FracStep, BDF
--! @param orderOrTheta		theta param if 'Theta', order if 'BDF'
--! @param startTime		start time point
--! @param endTime			end time point")
--! @param maxStepSize		maximal step sized used
--! @param minStepSize 		(optional) minimal step sized used	
--! @param reductionFactor 	(optional) factor by which the step size is
--! 						reduced, if the problem was not solved.
--! 						Iterated until minStepSize is reached.
--! @param bFinishTimeStep 	(optional) boolean if finish_timestep should be 
--! 						called or not.
--! @param useCheckpointing (optional) if true, use checkpointing.
--! @param postProcess		(optional) if passed, will be called after every time step.
--! @param startTSNo		(optional) time step number of the initial condition (normally 0).
--! @param endTSNo			(optional) if passed, stop after the time step with this number.
function util.SolveNonlinearTimeProblem(
	u,
	domainDisc,
	newtonSolver,
	out,
	filename,
	timeScheme,
	orderOrTheta,
	startTime,
	endTime,
	maxStepSize,
	minStepSize,
	reductionFactor,
	bFinishTimeStep,
	useCheckpointing,
	postProcess,
	startTSNo,
	endTSNo)
	
	if u == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: No grid function for the solution specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if domainDisc == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: No domain discretization specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if newtonSolver == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: No nonlin. solver specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if timeScheme == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: No time scheme specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if startTime == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: Start time not specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if endTime == nil and endTSNo == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: End time or number of steps not specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if maxStepSize == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: No max. time step specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end
	
	-- create time disc
	local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)
	
	
	-- print newtonSolver setup	
	print("SolveNonlinearTimeProblem, Newton Solver setup:")
	print(newtonSolver:config_string())

	-- start
	local time = startTime
	local step = startTSNo or 0
	local nlsteps = 0;
	
	
	if useCheckpointing then
		--- Read Checkpoint if necessary
		if util.HasParamOption("-restartWithCheckpoint") then
			cp = util.ReadCheckpoint(u)
			if cp ~= nil then
				time = cp.myData.time
				step = cp.myData.step
			end
		end
	end
	
	-- write start solution
	print(">> Writing start values")
	if type(out) == "function" then out(u, step, time) end
	if type(out) == "userdata" then out:print(filename, u, step, time) end

	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(u:clone(), time)

	-- update newtonSolver	
	newtonSolver:init(AssembledOperator(timeDisc, u:grid_level()))

	-- store old solution (used for reinit in multistep)
	local uOld
	if timeDisc:num_stages() > 1 then uOld = u:clone() end			

	-- set order for bdf to 1 (initially)
	if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end

	-- bound on t-stepper from machine precision (conservative)
	relPrecisionBound = 1e-12

	-- stop if size of remaining t-domain (relative to `maxStepSize`) lies below `relPrecisionBound`
	while (endTime == nil or ((time < endTime) and ((endTime-time)/maxStepSize > relPrecisionBound))) and ((endTSNo == nil) or (step < endTSNo)) do
		step = step+1
		print("++++++ TIMESTEP " .. step .. " BEGIN (current time: " .. time .. ") ++++++");

		-- initial t-step size
		local currdt = maxStepSize
		-- adjust in case of over-estimation
		if time+currdt > endTime then currdt = endTime-time end
		-- adjust if size of remaining t-domain (relative to `maxStepSize`) lies below `relPrecisionBound`
		if ((endTime-(time+currdt))/maxStepSize < relPrecisionBound) then currdt = endTime-time end
		
		-- try time step
		local bSuccess = false;	
		while bSuccess == false do
			TerminateAbortedRun()

			print("++++++ Time step size: "..currdt);

			-- get old solution if multistage
			if timeDisc:num_stages() > 1 then
				VecScaleAssign(u, 1.0, uOld)
			end			

			for stage = 1, timeDisc:num_stages() do
				if timeDisc:num_stages() > 1 then
					print("      +++ STAGE " .. stage .. " BEGIN ++++++")
				end
						
				timeDisc:set_stage(stage)
			
				-- setup time Disc for old solutions and timestep size
				timeDisc:prepare_step(solTimeSeries, currdt)
				
				-- prepare newton solver
				if newtonSolver:prepare(u) == false then 
					print ("\n++++++ Newton solver failed."); exit();
				end 
				
				-- apply newton solver
				if newtonSolver:apply(u) == false then 
					currdt = currdt * reductionFactor;
					write("\n++++++ Newton solver failed. "); 
					write("Trying decreased stepsize " .. currdt .. ".\n");
				else
					bSuccess = true; 
				end
		
				-- check valid step size			
				if(bSuccess == false and currdt < minStepSize) then
					write("++++++ Time Step size "..currdt.." below minimal step ")
					write("size "..minStepSize..". Cannot solve problem. Aborting.");
					test.require(false, "Time Solver failed.")
				end
				
				-- start over again if failed
				if bSuccess == false then break end
				
				-- update new time
				time = timeDisc:future_time()
				nlsteps = nlsteps + newtonSolver:num_newton_steps() 	 
				--total_linsolver_calls()
			
				-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
				if timeScheme:lower() == "bdf" and step < orderOrTheta then
					print("++++++ BDF: Increasing order to "..step+1)
					timeDisc:set_order(step+1)
					solTimeSeries:push(u:clone(), time)
				else 
					oldestSol = solTimeSeries:oldest()
					VecScaleAssign(oldestSol, 1.0, u)
					solTimeSeries:push_discard_oldest(oldestSol, time)
				end
				
				if not (bFinishTimeStep == nil) and bFinishTimeStep then 
					timeDisc:finish_step_elem(solTimeSeries, u:grid_level()) 
				end
				
				if timeDisc:num_stages() > 1 then
					print("      +++ STAGE " .. stage .. " END   ++++++")
				end
			end
		end

		-- save this solution if multistage
		if timeDisc:num_stages() > 1 then
			uOld = u
		end			
		
		-- call post process
		if type(postProcess) == "function" then postProcess(u, step, time, currdt) end
		
		-- plot solution
		if type(out) == "function" then out(u, step, time) end
		if type(out) == "userdata" then out:print(filename, u, step, time) end
	
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
		
		if useCheckpointing then			
			----------------------------------------------------------
			--- Write Checkpoint if necessary
			util.WriteCheckpointIntervallic(u, time, {time=time, step=step, endTime=endTime})
			----------------------------------------------------------
		end
	end
	
	if useCheckpointing and  timeDisc:num_stages() > 1 then
		ug_warning("WARNING: Checkpointing won't work at the moment with timeDisc:num_stages() > 1")
	end 
	
	if type(out) == "userdata" then out:write_time_pvd(filename, u) end
	
	return step, time
end


--!
--! @param u 				[in] GridFunction with Startvalues, [out] Solution"
--! @param domainDisc 		Domain Discretization
--! @param solver       	Linear or Nonlinear Solver
--! @param out				a VTKOutput (pass nil for no output)
--! @param filename			filename for output
--! @param timeScheme   	Name of time step scheme:
--!     	            	Theta, ImplEuler, ExplEuler, Crank-Nicolson
--! 	                	Alexander, FracStep, BDF
--! @param orderOrTheta		theta param if 'Theta', order if 'BDF'
--! @param startTime		start time point
--! @param endTime			end time point")
--! @param maxStepSize		maximal step sized used
--! @param minStepSize 		(optinal) minimal step sized used	
--! @param reductionFactor 	(optinal) factor by which the step size is
--! 						reduced, if the problem was not solved.
--! 						Iterated until minStepSize is reached.
--! @param useCheckpointing (optional) if true, use checkpointing.
function util.SolveLinearTimeProblem(
	u,
	domainDisc,
	linSolver,
	out,
	filename,
	timeScheme,
	orderOrTheta,
	startTime,
	endTime,
	maxStepSize,
	minStepSize,
	reductionFactor,
	useCheckpointing)

	if u == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No grid function for the solution specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if domainDisc == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No domain discretization specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if linSolver == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No lin. solver specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if timeScheme == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No time scheme specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if startTime == nil or endTime == nil then
		print("SolveLinearTimeProblem: Illegal parameters: Start or end time not specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if maxStepSize == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No max. time step specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	-- check parameters
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end

	-- create time disc
	local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)
	
	if timeScheme:lower() == "alexander" or 
	   timeScheme:lower() == "fracstep" then	
	   	print("Only Theta, ImplEuler, ExplEuler, Crank-Nicolson, BDF supported.")
	   	print("For Multistage use Nonlinear time loop.")
		print("Aborting")
		exit()
	end
	
	
	-- print newtonSolver setup	
	print("SolveLinearTimeProblem, Linear Solver setup:")
	print(linSolver:config_string())
	
	-- start
	local time = startTime
	local step = 0
	
	if useCheckpointing then
		--- Read Checkpoint if necessary
		if util.HasParamOption("-restartWithCheckpoint")  then
			cp = util.ReadCheckpoint(u)
			if cp ~= nil then
				time = cp.myData.time
				step = cp.myData.step
			end
		end
	end	
	
	-- write start solution
	print(">> Writing start values")
	if not (out==nil) then out:print(filename, u, step, time) end
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(u:clone(), time)
	local gl = u:grid_level()

	-- matrix and vectors
	local A = AssembledLinearOperator(timeDisc, gl)
	local b = u:clone()

	-- set order for bdf to 1 (initially)
	if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end

	local assembled_dt = nil
	
	while time < endTime do
		step = step + 1
		print("++++++ TIMESTEP "..step.." BEGIN (current time: " .. time .. ") ++++++");
	
		-- initial time step size
		-- assure, that not reaching beyond end of interval and care for round-off
		local currdt = maxStepSize
		if time+currdt > endTime then currdt = endTime - time end
		if ((endTime - (time+currdt))/currdt) < 1e-8 then currdt = endTime - time end
		
		-- try time step
		local bSuccess = false;	
		while bSuccess == false do
			TerminateAbortedRun()
			print("++++++ Time step size: "..currdt);

			-- reassemble matrix if necessary
			if not(currdt == assembled_dt) then 
				print("++++++ Assembling Matrix/Rhs for step size "..currdt); 
				timeDisc:prepare_step(solTimeSeries, currdt)
				timeDisc:assemble_linear(A, b, gl)
				linSolver:init(A, u)
				assembled_dt = currdt
			else
				timeDisc:prepare_step(solTimeSeries, currdt)
				timeDisc:assemble_rhs(b, gl)
			end
			
			-- apply linear solver
			if linSolver:apply(u,b) == false then 
				currdt = currdt * reductionFactor;
				write("\n++++++ Linear solver failed. "); 
				write("Trying decreased stepsize " .. currdt .. ".\n");
			else
				bSuccess = true; 
			end
	
			-- check valid step size			
			if(bSuccess == false and currdt < minStepSize) then
				write("++++++ Time Step size "..currdt.." below minimal step ")
				write("size "..minStepSize..". Cannot solve problem. Aborting.");
				test.require(false, "Time Solver failed.")
			end
		end
		
		-- update new time
		time = solTimeSeries:time(0) + currdt
		
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
		if timeScheme:lower() == "bdf" and step < orderOrTheta then
			print("++++++ BDF: Increasing order to "..step+1)
			timeDisc:set_order(step+1)
			solTimeSeries:push(u:clone(), time)
		else 
			oldestSol = solTimeSeries:oldest()
			VecScaleAssign(oldestSol, 1.0, u)
			solTimeSeries:push_discard_oldest(oldestSol, time)
		end

		-- plot solution
		if not (out==nil) then out:print(filename, u, step, time) end
		--SaveVectorForConnectionViewer(u, filename.."_t"..step..".vec")
			
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
		
		if useCheckpointing then
			----------------------------------------------------------
			--- Write Checkpoint if necessary
			util.WriteCheckpointIntervallic(u, time, {time=time, step=step, endTime=endTime})
			----------------------------------------------------------
		end
	end
	
	if not (out==nil) then out:write_time_pvd(filename, u) end
	
	if useCheckpointing and  timeDisc:num_stages() > 1 then
		ug_warning("WARNING: Checkpointing won't work at the moment with timeDisc:num_stages() > 1")
	end 
end

--! Time stepping with the adaptive step size. Returns number of time steps done and the last time.
--! @param u 				[in] GridFunction with Startvalues, [out] Solution"
--! @param domainDisc 		Domain Discretization
--! @param newtonSolver		Nonlinear Solver
--! @param out				a VTKOutput (pass nil for no output)
--! @param filename			filename for output
--! @param startTime		start time point
--! @param endTime			end time point
--! @param dt				time step
--! @param maxStepSize		maximal step sized used
--! @param minStepSize 		(optional) minimal step sized used
--! @param adaptiveStepInfo adaptive stepping controls
--! @param bFinishTimeStep 	(optional) boolean if finish_timestep should be 
--! 						called or not.
--! @param postProcess		(optional) if passed, will be called after every time step.
--! @param startTSNo		(optional) time step number of the initial condition (normally 0).
--! @param endTSNo			(optional) if passed, stop after the time step with this number.
function util.SolveNonlinearProblemAdaptiveTimestep(
	u,
	domainDisc,
	newtonSolver,
	out,
	filename,
	startTime,
	endTime,
	dt,
	minStepSize,
	maxStepSize,
	adaptiveStepInfo, 
	--reductionFactor,
	--tol,
	bFinishTimeStep,
	postProcess,
	startTSNo,
	endTSNo)

	local doControl = true
	local doExtrapolation = true
	local timeScheme = "ImplEuler"
	local orderOrTheta = 1.0

	if u == nil then
		print("SolveNonlinearProblemAdaptiveTimestep: Illegal parameters: No grid function for the solution specified.")
		exit()
	end

	if domainDisc == nil then
		print("SolveNonlinearProblemAdaptiveTimestep: Illegal parameters: No domain discretization specified.")
		exit()
	end

	if newtonSolver == nil then
		print("SolveNonlinearProblemAdaptiveTimestep: Illegal parameters: No nonlin. solver specified.")
		exit()
	end

	if startTime == nil then
		print("SolveNonlinearProblemAdaptiveTimestep: Illegal parameters: Start time not specified.")
		exit()
	end

	if endTime == nil and endTSNo == nil then
		print("SolveNonlinearProblemAdaptiveTimestep: Illegal parameters: End time or number of steps not specified.")
		exit()
	end
 
	if maxStepSize == nil then
		print("SolveNonlinearProblemAdaptiveTimestep: Illegal parameters: No max. time step specified.")
		exit()
	end

	-- read adaptive stuff
	local tol = adaptiveStepInfo["TOLERANCE"]
	local red = adaptiveStepInfo["REDUCTION"]
	local inc_fac = adaptiveStepInfo["INCREASE"]
	local safety_fac = adaptiveStepInfo["SAFETY"]
	local errorEst = adaptiveStepInfo["ESTIMATOR"]

	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end
	
	if red == nil then red = 0.5 end   -- reduction of time step
	if inc_fac == nil then inc_fac = 1.5 end   -- increase of time step
	
	if errorEst == nil then errorEst = Norm2Estimator() end
	if tol == nil then tol = 1e-3 end
	if safety_fac == nil then safety_fac = 0.8 end   -- safety factor
	
	
	-- create time disc
	local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)
	
	-- print newtonSolver setup	
	print("SolveNonlinearTimeProblem, Newton Solver setup:")
	print(newtonSolver:config_string())
		
	
	-- start
	local time = startTime
	local step = startTSNo or 0
	local nlsteps = 0;
	
	-- write start solution
	print(">> Writing start values")
	if type(out) == "function" then out(u, step, time) end
	if type(out) == "userdata" then out:print(filename, u, step, time) end
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	local old = u:clone()	
	solTimeSeries:push(old, time)

	-- update newtonSolver	
	local nlop = AssembledOperator(timeDisc, u:grid_level())
	newtonSolver:init(nlop)
	
	
	-- extra options for adaptive scheme
	
	local timeDisc2
	local solTimeSeries2
	local nlop2
	
	local timex
	local u2, old2, aux
	
	
	if (doControl) then 
		 
		 timeDisc2 = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)	
		 nlop2 = AssembledOperator(timeDisc2, u:grid_level())	
		
		 old2 = old:clone()  -- second solution
		 u2 = old:clone()  -- second solution
		 aux = old:clone() 
		
		 -- time series
		solTimeSeries2 = SolutionTimeSeries()
		solTimeSeries2:push(old2, time)
		if (doExtrapolation) then	
			-- Aitken-Neville-type time	extrapolation
			timex = AitkenNevilleTimex({1,2}, errorEst)
		end
	end		
			
	-- set order for bdf to 1 (initially)
	if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end
	local currdt = dt
		
	while ((endTime == nil) or (time < endTime)) and ((endTSNo == nil) or (step < endTSNo)) do
		step = step + 1
		print("++++++ TIMESTEP "..step.." BEGIN (current time: " .. time .. ") "..endTime.."++++++");
	
		-- initial time step size
		-- assure, that not reaching beyond end of interval
		currdt = math.min(currdt, endTime-time);
		
		-- try time step
		local bSuccess = false;	
		while bSuccess == false do
			TerminateAbortedRun()
			print("++++++ Time step size: "..currdt);

			-- setup time Disc for old solutions and timestep size
			timeDisc:prepare_step(solTimeSeries, currdt)
			newtonSolver:init(nlop)
			
			-- prepare newton solver
			if newtonSolver:prepare(u) == false then 
				print ("\n++++++ Newton solver failed."); exit();
			end 
				
			-- apply newton solver
			if newtonSolver:apply(u) == false then 
				currdt = currdt * red;
				write("\n++++++ Newton solver failed. "); 
				write("Trying decreased stepsize " .. currdt .. ".\n");
			else
				bSuccess = true; 
			end
		
			-- time step control
			-- compute solution u2
			if (bSuccess and doControl) then 
			
				newtonSolver:init(nlop2)
				----------------------------
				-- step 1/2
				print("Control 1/2:");
		
				-- time2 = time-tau/2;                        			-- intermediate time step 

				--VecScaleAdd2(u2, 1.0-0.5*currdt, old, 0.5*currdt, u);   -- w/ linear interpolation  (first guess)
				timeDisc2:prepare_step(solTimeSeries2, 0.5*currdt)
				if newtonSolver:prepare(u2) == false then
					print ("Newton solver failed at step "..step.."+1/2.");
					exit();
				end 
				if newtonSolver:apply(u2) == false then 
					print ("Newton solver failed at step "..step.."+1/2.");
				end 
		
				----------------------------
				-- step 2/2	
				print("Control 2/2:");
		
				-- push back solution
				time2 = time + 0.5*currdt
				tmp2 = solTimeSeries2:oldest()                  -- solution at time t 
				VecScaleAssign(tmp2, 1.0, u2)                   -- is removed and replaced by u2(t+tau/2)
				solTimeSeries2:push_discard_oldest(tmp2, time2) -- re-insert as new solution 
																			-- (now old2 is discarded)
		
				timeDisc2:prepare_step(solTimeSeries2, 0.5*currdt)
				if newtonSolver:prepare(u2) == false then
					print ("Newton solver failed at step "..step.."+2/2.");
					exit();
				end 
				if newtonSolver:apply(u2) == false then
				 	print ("Newton solver failed at step "..step.."+2/2."); 
				 	bSuccess = false;
				end 
		
				if (doExtrapolation) then
						
					-- extrapolation (more accurate)
					timex:set_solution(u, 0)
					timex:set_solution(u2, 1)
					timex:apply()
					eps = timex:get_error_estimate()
		
				else
					-- no extrapolation (less accurate)
					VecScaleAdd2(aux, 1.0, u, -1.0, u2);
					local l2_fine_est = L2Norm(aux, "c", 2);
					eps = 2.0*l2_fine_est;
					l2_ex_err = "---";
					l2_ex_est = "---";
				end
		
				-------------------------
				-- Adaptive step control
				-------------------------
				local lambda = math.pow(safety_fac*tol/eps, 0.5) -- (eps<=tol) implies (tol/eps >= 1) 
				dtEst = lambda*currdt  
				print("dtEst= "..dtEst..", eps="..eps..", tol = " ..tol..", fac = "..lambda)
		
				-- determine potential new step size
				currdt = math.min(dtEst, inc_fac*currdt, maxStepSize)

				if (eps <= tol) then 
					-- accept
					print ("ACCEPTING solution, dtnew="..currdt);
					bAcceptStep = true;
					bSuccess =true;
		 			
				else
	    			-- discard
	    			print ("DISCARDING solution, dtnew="..currdt);
	    			bAcceptStep = false;
	    			bSuccess = false;
	    			
	    			-- reset initial guess (important for non-linear coarse problem)
	    			VecScaleAssign(u, 1.0, old)
	    			
	    			-- reset solTimeSeries2
	    			utmp = solTimeSeries2:oldest()
	   				VecScaleAssign(utmp, 1.0, old)
	   				solTimeSeries2:push_discard_oldest(utmp, time)
				end
			end -- if (bSuccess and doControl)
		
		
			-- check valid step size			
			if(bSuccess == false and currdt < minStepSize) then
				write("++++++ Time Step size "..currdt.." below minimal step ")
				write("size "..minStepSize..". Cannot solve problem. Aborting.");
				test.require(false, "Time Solver failed.")
			end
				
			-- start over again if failed
			if bSuccess == false then break end
				
			-- update new time
			time = timeDisc:future_time()
			nlsteps = nlsteps + newtonSolver:num_newton_steps() 	 
			
			-- push oldest solutions with new values to front, oldest sol pointer is poped from end		
			oldestSol = solTimeSeries:oldest()
			VecScaleAssign(oldestSol, 1.0, u2)
			solTimeSeries:push_discard_oldest(oldestSol, time)
			
			if (doControl) then
				-- do the same for second 
				oldestSol = solTimeSeries2:oldest()
				VecScaleAssign(oldestSol, 1.0, u2)
				solTimeSeries2:push_discard_oldest(oldestSol, time)
			end
				
			if not (bFinishTimeStep == nil) and bFinishTimeStep then 
				timeDisc:finish_step_elem(solTimeSeries, u:grid_level()) 
			end
				
		end	-- while bSuccess == false
		
		-- call post process
		if type(postProcess) == "function" then postProcess(u, step, time, currdt) end
		
		
		-- plot solution
		if type(out) == "function" then out(u, step, time) end
		if type(out) == "userdata" then out:print(filename, u, step, time) end
	--	print("Integral("..time..")="..Integral(massLinker, u, "INNER2"))
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ", nlsteps: "..nlsteps..") ++++++");
	end
	
	if type(out) == "userdata" then out:write_time_pvd(filename, u) end
	
	return step, time
end


function util.SolveNonlinearProblemAdaptiveLimex(
	u,
	domainDisc,
	newtonSolver,
	out,
	filename,
	startTime,
	endTime,
	dt,
	minStepSize,
	maxStepSize,
	adaptiveStepInfo,
	--reductionFactor,
	--tol,
	bFinishTimeStep,
	postProcess,
	startTSNo,
	endTSNo)


	-- read adaptive stuff
	local tol = adaptiveStepInfo["TOLERANCE"]
	local red = adaptiveStepInfo["REDUCTION"]
	local inc_fac = adaptiveStepInfo["INCREASE"]
	local safety_fac = adaptiveStepInfo["SAFETY"]
	local errorEst = adaptiveStepInfo["ESTIMATOR"]
  
	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end

	if red == nil then red = 0.5 end   -- reduction of time step
	if inc_fac == nil then inc_fac = 1.5 end   -- increase of time step
  
	if errorEst == nil then 
		print "WARNING: Error estimator not set. Default is euclidean norm! "
		errorEst = Norm2ErrorEst() 
	end
	if tol == nil then 
		tol = 1e-3 
		print ("WARNING: Using default tolerance "..tol)
	end
	if safety_fac == nil then safety_fac = 0.8 end   -- safety factor


	local doControl = true
	local doExtrapolation = true
	local timeScheme = "ImplEuler"
	local orderOrTheta = 1.0
	if u == nil or domainDisc == nil or newtonSolver == nil or timeScheme == nil
		or startTime == nil or endTime == nil or maxStepSize == nil then
		print("Wrong usage found. Please specify parameters as below:")
		
		if (u == nil) then print ("Did not find u!"); end;
		if (domainDisc == nil) then print ("Did not find domainDisc!"); end;
		if (timeScheme == nil) then print ("Did not find timeScheme!"); end;
		
		if (startTime == nil) then print ("Did not find endTime!"); end;
    if (endTime == nil) then print ("Did not find endTime!"); end;
    if (maxStepSize == nil) then print ("Did not find maxStepSize!"); end;
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end


	print ("maxStepSize ="..maxStepSize)
	print ("minStepSize ="..minStepSize)

	print ("startTime ="..startTime)
	print ("endTime ="..endTime)
	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end
	
	-- create time disc
	local timeDisc = LinearImplicitEuler(domainDisc)
	
	-- print newtonSolver setup	
	print("SolveNonlinearTimeProblem, Newton Solver setup:")
	print(newtonSolver:config_string())
			
	-- start
	local time = startTime
	local step = startTSNo or 0
	local nlsteps = 0;
	
	-- write start solution
	print(">> Writing start values")
	if type(out) == "function" then out(u, step, time) end
	if type(out) == "userdata" then out:print(filename, u, step, time) end
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	local old = u:clone()	
	solTimeSeries:push(old, time)

	-- update newtonSolver	
	local nlop = AssembledOperator(timeDisc, u:grid_level())
	newtonSolver:init(nlop)
	
	-- extra options for adaptive scheme
	local timeDisc2
	local solTimeSeries2
	local nlop2
	
	local timex
	local u2, old2, aux
	
	
	if (doControl) then 
		 
		 timeDisc2 = LinearImplicitEuler(domainDisc)	
		 nlop2 = AssembledOperator(timeDisc2, u:grid_level())	
		
		 old2 = old:clone()  -- second solution
		 u2 = old:clone()  -- second solution
		 aux = old:clone() 
		
		 -- time series
		solTimeSeries2 = SolutionTimeSeries()
		solTimeSeries2:push(old2, time)
		if (doExtrapolation) then	
			-- Aitken-Neville-type time	extrapolation
			timex = AitkenNevilleTimex({1,2}, errorEst)
		end
	end		
			
	-- set order for bdf to 1 (initially)
	-- if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end
	local currdt = dt
		
	while ((endTime == nil) or (time < endTime)) and ((endTSNo == nil) or (step < endTSNo)) do
		step = step + 1
		print("++++++ TIMESTEP "..step.." BEGIN (current time: " .. time .. ") "..endTime.."++++++");
	
		-- initial time step size
		-- assure, that not reaching beyond end of interval
		currdt = math.min(currdt, endTime-time);
		
		-- try time step
		local bSuccess = false;	
		while bSuccess == false do
			TerminateAbortedRun()
			print("++++++ Time step size: "..currdt);

			-- setup time Disc for old solutions and timestep size
			timeDisc:prepare_step(solTimeSeries, currdt)
			newtonSolver:init(nlop)
			
			-- prepare newton solver
			if newtonSolver:prepare(u) == false then 
				print ("\n++++++ Newton solver failed."); exit();
			end 
				
			-- apply newton solver
			if newtonSolver:apply(u) == false then 
				currdt = currdt * reductionFactor;
				write("\n++++++ Newton solver failed. "); 
				write("Trying decreased stepsize " .. currdt .. ".\n");
			else
				bSuccess = true; 
			end
		
			-- time step control
			-- compute solution u2
			if (bSuccess and doControl) then 
			
				newtonSolver:init(nlop2)
				----------------------------
				-- step 1/2
				print("Control 1/2:");
		
				-- time2 = time-tau/2;                        			-- intermediate time step 
				--VecScaleAdd2(u2, 1.0-0.5*currdt, old, 0.5*currdt, u);   -- w/ linear interpolation  (first guess)

				timeDisc2:prepare_step(solTimeSeries2, 0.5*currdt)
				if newtonSolver:prepare(u2) == false then print ("Newton solver failed at step "..step.."+1/2."); exit(); end 
				if newtonSolver:apply(u2) == false then 
					print ("Newton solver failed at step "..step.."+1/2.");
				end 
		
				----------------------------
				-- step 2/2	
				print("Control 2/2:");
		
				-- push back solution
				tmp2 = solTimeSeries2:oldest()                     			-- solution at time t 
				VecScaleAssign(tmp2, 1.0, u2)                       		-- is removed and replaced by u2(t+tau/2)
				solTimeSeries2:push_discard_oldest(tmp2, time + 0.5*currdt) -- re-insert as new solution 
																			-- (now old2 is discarded)
		
				timeDisc2:prepare_step(solTimeSeries2, 0.5*currdt)
				if newtonSolver:prepare(u2) == false then print ("Newton solver failed at step "..step.."+2/2."); exit(); end 
				if newtonSolver:apply(u2) == false then
				 	print ("Newton solver failed at step "..step.."+2/2."); 
				 	bSuccess = false;
				end 
					
			
				
				
				if (doExtrapolation) then		
					-- extrapolation (more accurate)
					timex:set_solution(u, 0)
					timex:set_solution(u2, 1)
					timex:apply()
					eps = timex:get_error_estimate()
					--out:print("ExSol.vtu", u2, step, time) 
					
				else
					-- no extrapolation (less accurate)
					VecScaleAdd2(aux, 1.0, u, -1.0, u2);
					local l2_fine_est = L2Norm(aux, "c", 2);
					eps = 2.0*l2_fine_est;
					l2_ex_err = "---";
					l2_ex_est = "---";
				end
		
				--vtk = VTKOutput()
				--vtk:select_nodal(GridFunctionNumberData(u, "c"), "CNodal")
				--vtk:select_nodal(GridFunctionNumberData(u, "p"), "PNodal")
        
				--out:print("CoarseLimexSol.vtu", u, step, time) 
				--out:print("FineLimexSol.vtu", u2, step, time) 
		
				--print("TIME_ERROR (t="..time+tau..", |u|="..l2normB.. ") :\t" .. l2_coarse_err .. "\t"..l2_fine_err .. "\t"..l2_ex_err .. "\|  ---\t"..l2_fine_est .. "\t"..l2_ex_est.. "\t"..eps)
				--file:write(time+tau.."\t" ..l2normB.. "\t" .. l2_coarse_err .. "\t"..l2_fine_err .. "\t"..l2_ex_err .. "\t"..l2_fine_est .. "\t"..l2_ex_est.. "\t"..eps.."\n")
				--file:flush()
		
				-------------------------
				-- Adaptive step control
				-------------------------
					
				--dtEst = math.pow(0.9*tol/eps, 0.5)*currdt  -- (eps<=tol) implies (tol/eps >= 1) 
				--print("dtEst= "..dtEst..", eps="..eps..", tol = " ..tol..", fac = "..math.pow(0.9*tol/eps, 0.5))
				local lambda = math.pow(safety_fac*tol/eps, 0.5) -- (eps<=tol) implies (tol/eps >= 1) 
				dtEst = lambda*currdt  
				print("dtEst= "..dtEst..", eps="..eps..", tol = " ..tol..", fac = "..lambda)
		
				-- determine potential new step size
				currdt = math.min(dtEst, inc_fac*currdt, maxStepSize)

				if (eps <= tol) then 
					-- accept
					print ("ACCEPTING solution, dtnew="..currdt);
					bAcceptStep = true;
					bSuccess =true;
		 			
				else
	    			-- discard
	    			print ("DISCARDING solution, dtnew="..currdt);
	    			bAcceptStep = false;
	    			bSuccess = false;
	    			
	    			-- reset solTimeSeries2
	    			utmp = solTimeSeries2:oldest()
	   				VecScaleAssign(utmp, 1.0, old)
	   				solTimeSeries2:push_discard_oldest(utmp, time)
				end
			end
		
		
			-- check valid step size			
			if(bSuccess == false and currdt < minStepSize) then
				write("++++++ Time Step size "..currdt.." below minimal step ")
				print("size "..minStepSize..". Cannot solve problem. Aborting.");
				test.require(false, "Time Solver failed.")
			end
				
			-- start over again if failed
			if bSuccess == false then break end
				
			-- update new time
			time = timeDisc:future_time()
			nlsteps = nlsteps + newtonSolver:num_newton_steps() 	 
			
			-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
			oldestSol = solTimeSeries:oldest()
			VecScaleAssign(oldestSol, 1.0, u2)
			solTimeSeries:push_discard_oldest(oldestSol, time)
			
			if (doControl) then
				-- do the same for second 
				oldestSol = solTimeSeries2:oldest()
				VecScaleAssign(oldestSol, 1.0, u2)
				solTimeSeries2:push_discard_oldest(oldestSol, time)
			end
			
				
			if not (bFinishTimeStep == nil) and bFinishTimeStep then 
				timeDisc:finish_step_elem(solTimeSeries, u:grid_level()) 
			end
				
		end		
		
		-- call post process
		if type(postProcess) == "function" then postProcess(u, step, time, currdt) end
		
		-- plot solution
		if type(out) == "function" then out(u, step, time) end
		if type(out) == "userdata" then out:print(filename, u, step, time) end

		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ", nlsteps: "..nlsteps..") ++++++");
	end
	
	if type(out) == "userdata" then out:write_time_pvd(filename, u) end
	
	return step, time
end




--[[!
\}
]]--
