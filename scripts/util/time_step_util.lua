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

--! Parses a callback object and attachs the corresponding callbacks to a TimeIntegratorObserver object.
--! @param timeIntegratorSubject a TimeIntegratorObserver object to attach callbacks to
--! @luaobject the object by the user, to be parsed.
--!	a)	If this is a function then it is called after
--!		solving the non-linear problem in EVERY STAGE
--!		of the time-stepping scheme;
--!	b)	if this is a table, it can contain 4 optional functions:
--!		prepareTimeStep to call before the time step,
--!		preProcess to call before the non-linear solver in EVERY STAGE,
--!		postProcess as in a),
--!		finalizeTimeStep to call after the time step,
--!		rewindTimeStep to call if some of the computations of the time step failed,
--!		Arguments of the functions are: (u, step, time, dt)
--!		u, time: old before the solver, new after it
--!	c)	This can also be a list of C++ objects which inherit ITimeIntegratorObserver
--!		or Lua Callbacks created by util.LuaCallbackHelper:create()
function util.ParseTimeIntegratorCallbacks(timeIntegratorSubject, luaobject)
	-- helper function 
	function attachObservers(attachFct, element, donttraverse)
		donttraverse = donttraverse or false						

		-- user specified a function
		if type(element) == "function" then
			local cb = util.LuaCallbackHelper:create(element)
			attachFct(timeIntegratorSubject, cb.CPPCallback)

		-- user specified a class created by util.LuaCallbackHelper.create()
		elseif type(element) == "table" and element.CPPCallback ~= nil then	
			attachFct(timeIntegratorSubject, element.CPPCallback)

		-- user specified an cpp class, lets hope it implements ITimeIntegratorObserver
		elseif type(element) == "userdata" then
			attachFct(timeIntegratorSubject, element)

		-- smells like an array
		elseif type(element) == "table" and element[1] ~= nil and donttraverse == false then
			for i, child in ipairs(element) do		
				attachObservers(attachFct, child, true)
			end
		end
	end

	attachObservers(timeIntegratorSubject.attach_postprocess_observer, luaobject)		

	if luaobject.preProcess ~= nil then
		attachObservers(timeIntegratorSubject.attach_preprocess_observer, luaobject.preProcess)
	end
	if luaobject.postProcess ~= nil then
		attachObservers(timeIntegratorSubject.attach_postprocess_observer, luaobject.postProcess)
	end
	if luaobject.prepareTimeStep ~= nil then
		attachObservers(timeIntegratorSubject.attach_init_observer, luaobject.prepareTimeStep)
	end
	if luaobject.finalizeTimeStep ~= nil then
		attachObservers(timeIntegratorSubject.attach_finalize_observer, luaobject.finalizeTimeStep)
	end
	if luaobject.rewindTimeStep ~= nil then
		attachObservers(timeIntegratorSubject.attach_rewind_observer, luaobject.rewindTimeStep)
	end
	if luaobject.startProcess ~= nil then
		attachObservers(timeIntegratorSubject.attach_start_observer, luaobject.startProcess)
	end
	if luaobject.endProcess ~= nil then
		attachObservers(timeIntegratorSubject.attach_end_observer, luaobject.endProcess)
	end
end

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

-- split from SolveNonlinearTimeProblem below
function SolveNonlinearTimeProblemParams(
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
	endTSNo,
	newtonLineSearchFallbacks,
	additionalFinishedConditions)

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

	if endTime == nil and endTSNo == nil and (additionalFinishedConditions == nil or #additionalFinishedConditions == 0)  then
		print("SolveNonlinearTimeProblem: Illegal parameters: End time, number of steps or other finished conditions not specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if maxStepSize == nil then
		print("SolveNonlinearTimeProblem: Illegal parameters: No max. time step specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end
end

--! Time stepping with a fixed step size. Returns number of time steps done and the last time.
--! @param u 				[in] GridFunction with Startvalues, [out] Solution"
--! @param domainDisc		Domain Discretization
--! @param newtonSolver		Nonlinear Solver
--! @param out				userdata or function? a VTKOutput (pass nil for no output)
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
--! @param postProcess		(optional) if passed, this can be either a function or a table:
--!									a)	If this is a function then it is called after
--!										solving the non-linear problem in EVERY STAGE
--!										of the time-stepping scheme;
--!									b)	if this is a table, it can contain 4 optional functions:
--!										prepareTimeStep to call before the time step,
--!										preProcess to call before the non-linear solver in EVERY STAGE,
--!										postProcess as in a),
--!										finalizeTimeStep to call after the time step,
--!										rewindTimeStep to call if some of the computations of the time step failed,
--!										Arguments of the functions are: (u, step, time, dt)
--!										u, time: old before the solver, new after it
--!									c)	This can also be a list of C++ objects which inherit ITimeIntegratorObserver
--!										or Lua Callbacks created by util.LuaCallbackHelper:create()
--! @param startTSNo		(optional) time step number of the initial condition (normally 0).
--! @param endTSNo			(optional) if passed, stop after the time step with this number.
--! @param newtonLineSearchFallbacks	(optional) Sequence of line-search objects.
--!										Each time the newton-solver fails in a time-step,
--!										the next line-search object from the sequence is used.
--! @param additionalFinishedConditions (optional) More objects of type IFinishedConditions to add to the FinishedConditions
--!										defined by endTSNo or endTime
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
	endTSNo,
	newtonLineSearchFallbacks,
	additionalFinishedConditions)

	SolveNonlinearTimeProblemParams(
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
		endTSNo,
		newtonLineSearchFallbacks,
		additionalFinishedConditions)
	
	-- attach the given callbacks to a TimeIntegratorSubject
	-- if this class is transcribed into c++, just inherit from TimeIntegratorSubject and attach the callbacks
	local callbackDispatcher = TimeIntegratorSubject()
	cplusplus = true;

	if cplusplus then
		-- parsed below
	elseif postProcess ~= nil then
		util.ParseTimeIntegratorCallbacks(callbackDispatcher, postProcess)
	end
		
	-- bound on t-stepper from machine precision (conservative)
	relPrecisionBound = 1e-12

	-- create the finished conditions
	local finishedTester = FinishedTester()
	if endTSNo ~= nil then
		maxStepsCondition = MaxStepsFinishedCondition(endTSNo)
		finishedTester:add_condition(maxStepsCondition)
	end
	if endTime ~= nil then
		temporalCondition = TemporalFinishedCondition(endTime,maxStepSize,relPrecisionBound)
		finishedTester:add_condition(temporalCondition)
	end
	if additionalFinishedConditions ~= nil then
		for i, fc in ipairs(additionalFinishedConditions) do	
			finishedTester:add_condition(fc)
		end
	end

	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end
	if reductionFactor >= 1 then
		print("SolveNonlinearTimeProblem: Illegal parameters: reductionFactor must be < 1.")
		exit()
	end
	
	-- create time disc -- BUG: requires luacpp
	local timeDisc = createTimeDisc(domainDisc, timeScheme, orderOrTheta)

	
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
	print(">> Writing start values using "..type(out))
	if type(out) == "function" then out(u, step, time)
	elseif type(out) == "userdata" then out:print(filename, u, step, time)
	end

	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(u:clone(), time)

	-- update newtonSolver	
	newtonSolver:init(AssembledOperator(timeDisc, u:grid_level()))

	-- store old solution (used for reinit in multistep)
	local uOld
	if minStepSize <= maxStepSize * reductionFactor or newtonLineSearchFallbacks ~= nil then
		uOld = u:clone()
	end
	-- TODO: This can be optimized because the "old" solution is stored in
	-- solTimeSeries. Note that 'u' keeps not (always) the new solution but
	-- an iterate of the Newton method. If the Newton method fails, one should
	-- reset this iterate from the old solution: The Newton method needs a good
	-- initial approximation and its convergence depends on it strongly!

	-- set order for bdf to 1 (initially)
	if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end

	-- bound on t-stepper from machine precision (conservative)
	relPrecisionBound = 1e-12

	-- not needed in c++ version
	local defaultLineSearch = newtonSolver:line_search()

	local last_dt = currdt

	if cplusplus then -- c++ version
		print("NonLinear c++rework...")
		RequiredPlugins({"Luacpp"})

		loop = SolveNonlinearTimeProblemOuterLoop()
		if postProcess ~= nil then
			util.ParseTimeIntegratorCallbacks(loop, postProcess)
		end
		loop:setFinishedTester(finishedTester)
		loop:setNewtonSolver(newtonSolver)
		loop:setTimeDisc(timeDisc) -- move "createTimeDisc(domainDisc, timeScheme, orderOrTheta)" to constructor?
		                           -- also, call NewtonSolver:init from there?
		loop:setStep(step)
		loop:setOrderOrTheta(orderOrTheta)
		loop:setReductionFactor(reductionFactor)
		loop:setMinStepSize(minStepSize)
		loop:setMaxStepSize(maxStepSize)
		loop:setTimeScheme(timeScheme)
		loop:setBFinishTimeStep(bFinishTimeStep); -- what if nil? does this work?

		loop:setEndTime(endTime)
		loop:setRelPrecisionBound(relPrecisionBound)
		if util.debug_writer ~= nil then
			loop:setDebugWriterContext(util.debug_writer)
		end
		loop:attachU(u)

		function MyLuaCallback(step, time, currdt, dummy)
			--- eek, u is shared between callback and solver.
			-- but still passed as argument?
			out(u, step, time, currdt)
			return 0.
		end

		ocb = LuaCallbackObserver() -- is an ITimeIntegratorObserver...?
		-- print("ocb type ", type(ocb)) -- is "userdata"
		-- print("out type ", type(out)) -- is "userdata"
		-- ocb can be function or userdata?
		--
		-- does not work.
		-- thing = TimeIntegratorObserver()
		-- util.ParseTimeIntegratorCallbacks(thing, out)

		if(type(out) == "function") then
			-- print("setting out function", out)
			ocb:set_callback("MyLuaCallback")
		else
			-- print("incomplete. not a function", ocb);
			exit();
		end

		loop:setOutput(ocb)
--		loop:setCallback(callbackDispatcher)
		if minStepSize <= maxStepSize * reductionFactor or newtonLineSearchFallbacks ~= nil then
			loop:storeU()
		end
		loop:do_it()
		step = loop:get_step()
		time = loop:get_time()
		-- last_dt = loop:getLastDt()
		-- loop:getU(u) no. shared.
	else -- old version
		print("legacy SolveNonlinearTimeProblem")

	callbackDispatcher:notify_start(u, step, time, maxStepSize)

	while not finishedTester:is_finished(time, step) do
		step = step+1
		print("++++++ TIMESTEP " .. step .. " BEGIN (current time: " .. time .. ") ++++++");
		
		local solver_call = 0;
		
		-- initial t-step size
		local currdt = maxStepSize
		if endTime ~= nil then
			-- adjust in case of over-estimation
			if time+currdt > endTime then currdt = endTime-time end
			-- adjust if size of remaining t-domain (relative to `maxStepSize`) lies below `relPrecisionBound`
			if ((endTime-(time+currdt))/maxStepSize < relPrecisionBound) then currdt = endTime-time end
		end
		
		-- try time step
		local bSuccess = false;	
		while bSuccess == false do
			TerminateAbortedRun()

			print("++++++ Time step size: "..currdt);

			local newtonSuccess = false
			local newtonTry = 1
			newtonSolver:set_line_search(defaultLineSearch)

			-- try to solve the non-linear problem with different line-search strategies
			while newtonSuccess == false do
				-- get old solution if the restart with a smaller time step is possible
				if uOld ~= nil then
					VecAssign(u, uOld)
				end
				
				local pp_res -- result of pre- and postprocess routines
				
				pp_res = true
				if util.debug_writer ~= nil then
					util.debug_writer:enter_section ("TIMESTEP-"..step.."-Prepare-SolverCall-"..solver_call)
				end
				local pp_res = callbackDispatcher:notify_init_step(u, step, time, currdt)
				if util.debug_writer ~= nil then
					util.debug_writer:leave_section ()
				end
				if pp_res == false then
					print("\n++++++ Preparation of the time step failed.")
					break
				end
				
				for stage = 1, timeDisc:num_stages() do
					if timeDisc:num_stages() > 1 then
						print("      +++ STAGE " .. stage .. " BEGIN ++++++")
					end
					
					-- call pre process
					if util.debug_writer ~= nil then
						util.debug_writer:enter_section ("TIMESTEP-"..step.."-PreProcess-SolverCall-"..solver_call)
					end
					pp_res = callbackDispatcher:notify_preprocess_step(u, step, time, currdt)
					if util.debug_writer ~= nil then
						util.debug_writer:leave_section ()
					end
					if pp_res == false then -- i.e. not nil, not something else, but "false"!
						print("\n++++++ PreProcess failed.")
						newtonSuccess = false
						break
					end
					
						
					print("set_stage")
					timeDisc:set_stage(stage)
				
					-- setup time Disc for old solutions and timestep size
					timeDisc:prepare_step(solTimeSeries, currdt)
					
					-- enter debug section (if the debug_writer is specified)
					if util.debug_writer ~= nil then
						util.debug_writer:enter_section ("TIMESTEP-"..step.."-SolverCall-"..solver_call)
					end
			
					-- prepare newton solver
					if newtonSolver:prepare(u) == false then 
						print ("\n++++++ Newton solver failed."); exit();
					end 
					
					-- apply newton solver
					newtonSuccess = newtonSolver:apply(u)
						
					-- exit debug section (if the debug_writer is specified)
					if util.debug_writer ~= nil then
						util.debug_writer:leave_section ()
					end
			
					-- start over again if failed
					if newtonSuccess == false then break end
					
					-- call post process
					if util.debug_writer ~= nil then
						util.debug_writer:enter_section ("TIMESTEP-"..step.."-PostProcess-SolverCall-"..solver_call)
					end
					pp_res = callbackDispatcher:notify_postprocess_step(u, step, timeDisc:future_time(), currdt)
					if util.debug_writer ~= nil then
						util.debug_writer:leave_section ()
					end
					if pp_res == false then -- i.e. not nil, not something else, but "false"!
						print("\n++++++ PostProcess failed.")
						newtonSuccess = false
						break
					end

					--total_linsolver_calls()
				
					-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
					if timeScheme:lower() == "bdf" and step < orderOrTheta then
						print("++++++ BDF: Increasing order to "..step+1)
						timeDisc:set_order(step+1)
						solTimeSeries:push(u:clone(), timeDisc:future_time())
					else 
						local oldestSol = solTimeSeries:oldest()
						VecAssign(oldestSol, u)
						solTimeSeries:push_discard_oldest(oldestSol, timeDisc:future_time())
					end
					
					if not (bFinishTimeStep == nil) and bFinishTimeStep then 
						timeDisc:finish_step_elem(solTimeSeries, u:grid_level()) 
					end
					
					if timeDisc:num_stages() > 1 then
						print("      +++ STAGE " .. stage .. " END   ++++++")
					end
				end -- loop over the stages

				if newtonSuccess then
					if util.debug_writer ~= nil then
						util.debug_writer:enter_section ("TIMESTEP-"..step.."-Finalize-SolverCall-"..solver_call)
					end
					last_dt = currdt
					pp_res = callbackDispatcher:notify_finalize_step(u, step, time, currdt)
					if util.debug_writer ~= nil then
						util.debug_writer:leave_section ()
					end
					if pp_res == false then -- i.e. not nil, not something else, but "false"!
						write("\n++++++ Finalization of the time step failed.")
						newtonSuccess = false
						break
					end
				end
				
				if newtonSuccess == false and newtonLineSearchFallbacks ~= nil then
					if newtonLineSearchFallbacks[newtonTry] == nil or newtonSolver:last_num_newton_steps() == 0 then
						write("\n++++++ Adaptive Newton failed.")
						break
					else
						newtonSolver:set_line_search(newtonLineSearchFallbacks[newtonTry])
						write("Restarting Newton method with line search fallback " .. newtonTry .. ".\n");
						newtonTry = newtonTry + 1
					end
				else
					break
				end
			end -- while newtonSuccess == false

			-- call post process
			if newtonSuccess == false then
				currdt = currdt * reductionFactor;
				write("\n++++++ Newton method failed. ");
				write("Trying decreased stepsize " .. currdt .. ".\n");					
				if(bSuccess == false and currdt < minStepSize) then
					write("++++++ Time Step size "..currdt.." below minimal step ")
					write("size "..minStepSize..". Cannot solve problem. Aborting.");
					test.require(false, "Time Solver failed.")
				end
				bSuccess = false
				if rewindTimeStep ~= nil  then
					if util.debug_writer ~= nil then
						util.debug_writer:enter_section ("TIMESTEP-"..step.."-Rewind-SolverCall-"..solver_call)
					end
					callbackDispatcher:notify_rewind_step(u, step, time, currdt)
					if util.debug_writer ~= nil then
						util.debug_writer:leave_section ()
					end
				end
			else
				-- update new time
				time = timeDisc:future_time()
				nlsteps = nlsteps + newtonSolver:num_newton_steps() 	 
				bSuccess = true
			end
		
		end -- while bSuccess == false

		-- save this solution if the restard with a smaller time step is possible
		if uOld ~= nil then
			VecAssign (uOld, u)
		end			
		
		-- plot solution
		print("plot type ", type(out));
		if type(out) == "function" then
			out(u, step, time)
		elseif type(out) == "userdata" then
			out:print(filename, u, step, time)
		end
	
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
		
		if useCheckpointing then			
			----------------------------------------------------------
			--- Write Checkpoint if necessary
			util.WriteCheckpointIntervallic(u, time, {time=time, step=step, endTime=endTime})
			----------------------------------------------------------
		end
	end -- outer loop

	callbackDispatcher:notify_end(u, step, time, last_dt)

	end -- lua script version
	
	if useCheckpointing and  timeDisc:num_stages() > 1 then
		ug_warning("WARNING: Checkpointing won't work at the moment with timeDisc:num_stages() > 1")
	end 
	
	if type(out) == "userdata" then out:write_time_pvd(filename, u) end
	
	return step, time
end -- SolveNonlinearTimeProblem

function SolveLinearTimeProblemParams(
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
	useCheckpointing,
	postProcess,
	startTSNo,
	endTSNo)
	if u == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No grid function for the solution specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	local reassemble = false
	if domainDisc == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No domain discretization specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	elseif type(domainDisc) == "table" then
		reassemble = domainDisc.reassemble
		domainDisc = domainDisc.domainDisc
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

	if startTime == nil then
		print("SolveLinearTimeProblem: Illegal parameters: Start time not specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if endTime == nil and endTSNo == nil then
		print("SolveLinearTimeProblem: Illegal parameters: End time or number of steps not specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	if maxStepSize == nil then
		print("SolveLinearTimeProblem: Illegal parameters: No max. time step specified.")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end
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
--! @param postProcess		(optional) if passed, this can be either a function or a table:
--!									a)	If this is a function then it is called after
--!										solving the linear problem in every time step
--!										of the time-stepping scheme;
--!									b)	if this is a table, it can contain 2 optional functions:
--!										preProcess to call before the time step,
--!										postProcess as in a),
--!										Arguments of the functions are: (u, step, time, dt)
--!										u, time: old before the solver, new after it;
--!										there may be also values for the return codes:
--!										retValAtOK to return if the process ends normally
--!										retValAtMinStepSize to return if the min. dt is overcome
--!										retValAtSolver to return if the lin. solver failed
--! @param startTSNo		(optional) time step number of the initial condition (normally 0).
--! @param endTSNo			(optional) if passed, stop after the time step with this number.
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
	useCheckpointing,
	postProcess,
	startTSNo,
	endTSNo)

	SolveLinearTimeProblemParams(
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
	   useCheckpointing,
	   postProcess,
	   startTSNo,
	   endTSNo)
	
	local preProcess = nil
	local retValAtOK = nil
	local retValAtMinStepSize = nil
	local retValAtSolver = nil
	if postProcess ~= nil then
		if type(postProcess) ~= "function" then
			if type(postProcess) ~= "table" then
				print("SolveLinearTimeProblem: Illegal parameters: postProcess must be a function.") -- table??
				exit()
			end
			retValAtOK = postProcess.retValAtOK
			retValAtMinStepSize = postProcess.retValAtMinStepSize
			retValAtSolver = postProcess.retValAtSolver
			preProcess = postProcess.preProcess
			postProcess = postProcess.postProcess
		end
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
	
	-- set the level of verbosity (do not print much for only one time step)
	local verbose = true
	if startTSNo ~= nil and endTSNo ~= nil and endTSNo == startTSNo + 1 then
		verbose = false
	end
	
	-- print newtonSolver setup	
	if verbose then
		print("SolveLinearTimeProblem, Linear Solver setup:")
		print(linSolver:config_string())
	end
	
	-- start
	local time = startTime
	local step = startTSNo or 0
	
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
	if out ~= nil then
		print(">> Writing start values")
		out:print(filename, u, step, time)
	end
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(u:clone(), time)

	if true then -- c++ version
		RequiredPlugins({"Luacpp"})
		RequiredPlugins({"Limex"})
		print("Linear c++rework... time=", time, " endTSNo=", endTSNo, " step=", step, " endTime=", endTime,
		   " minStepSize=", minStepSize, " maxStepSize=", maxStepSize)

		-- matrix and vectors
		local gl = u:grid_level()
		local A = AssembledLinearOperator(timeDisc, gl)
		local b = u:clone()

		-- set order for bdf to 1 (initially)
		if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end

		local assembled_dt = nil

		if true then --- OuterLoop thing (old)
			loop = SolveLinearTimeProblemOuterLoop() -- ad hoc
			-- this is not actually used anywhere. incomplete/untested.
			-- what type is retVal?
			if retValAtOK == true then
				loop:setRetValAtOK(1);
			elseif retValAtOK == false then
				loop:setRetValAtOK(0);
			elseif retValAtOK ~= nil then
				loop:setRetValAtOK(retValAtOK);
			end
			if retValAtMinStepSize ~= nil then
				loop:setRetValAtMinStepSize(retValAtMinStepSize);
			end
			if retValAtSolver == false then
				loop:setRetValAtSolver(0);
			elseif retValAtSolver == true then
				loop:setRetValAtSolver(1);
			elseif retValAtSolver ~= nil then
				loop:setRetValAtSolver(retValAtSolver);
			end

			loop:set_time_disc(timeDisc) -- move "createTimeDisc(domainDisc, timeScheme, orderOrTheta)" to constructor?
			                             -- also, call NewtonSolver:init from there?
			if endTSNo ~= nil then
				print("incomplete, endTSNo", endTSNo)
				loop:setEndTSNo(endTSNo)
			end
			if out ~= nil then
				print("incomplete, setVTKOutput ", out)
				loop:setVTKOutput(out)
			end
			if filename ~= nil then
				print("incomplete, setFilename ", filename)
				loop:setFilename(filename)
			end

			loop:setStep(step)

			loop:set_linear_solver(linSolver)
			loop:setOrderOrTheta(orderOrTheta) -- already in timeDisc?
		elseif true then
			print("incomplete: use limex? -- which one?")

			loop = SimpleTimeIntegrator(timeDisc) -- does not work with linear solver?
			local inv = SimpleInvert()
			inv:set_linear_solver(linSolver)
			inv:set_time_disc(timeDisc)
			loop:set_solver(inv)

			-- loop = LinearTimeIntegrator(timeDisc) -- does not know min/max dt
			-- loop = TimeIntegratorLinearAdaptive(timeDisc) -- ???

		end

		if preProcess ~= nil then
			loop:setPreProcess(preProcess);
		end
		if postProcess ~= nil then
			loop:setPostProcess(postProcess);
		end

		if(reassemble == true) then
			loop:setReassemble(reassemble)
		end

		-- loop:setGl(gl)
		-- loop:set_time_step(step)
		--

		loop:set_time_step(minStepSize) -- really?
		loop:set_dt_min(minStepSize)
		loop:set_dt_max(maxStepSize)
		loop:set_reduction_factor(reductionFactor)

		if endTime == nil then
			print("incomplete, endTime")
			endTime = 1;
		end

		loop:init(u)

		loop:apply(u, endTime, u, startTime)
		step = 1 -- BUG loop:getStep()
		time = endTime -- BUG loop:getTime()
	else -- legacy script version

	print("legacy SolveLinearTimeProblem")
	-- matrix and vectors
	local gl = u:grid_level()
	local A = AssembledLinearOperator(timeDisc, gl)
	local b = u:clone()

	-- set order for bdf to 1 (initially)
	if timeScheme:lower() == "bdf" then timeDisc:set_order(1) end

	local assembled_dt = nil

	while ((endTime == nil) or (time < endTime)) and ((endTSNo == nil) or (step < endTSNo)) do
		step = step + 1
		if verbose then print("++++++ TIMESTEP "..step.." BEGIN (current time: " .. time .. ") ++++++") end
	
		-- initial time step size
		-- assure, that not reaching beyond end of interval and care for round-off
		local currdt = maxStepSize
		if endTime ~= nil then
			if time+currdt > endTime then currdt = endTime - time end
			if ((endTime - (time+currdt))/currdt) < 1e-8 then currdt = endTime - time end
		end
		
		-- try time step
		local bSuccess = false;	
		while bSuccess == false do
			TerminateAbortedRun()
			if verbose then print("++++++ Time step size: "..currdt) end

			if preProcess ~= nil then
				local pp_res = preProcess(u, step, time, currdt)
				if type(pp_res) == "boolean" and pp_res == false then -- i.e. not nil, not something else, but "false"!
					print("\n++++++ preProcess of the time step failed.")
					exit()
				end
			end
			
			-- reassemble matrix if necessary
			if reassemble or not(currdt == assembled_dt) then 
				print("++++++ Assembling Matrix/Rhs for step size "..currdt); 
				timeDisc:prepare_step(solTimeSeries, currdt)
				-- Remark: Do not use assemble_linear here: it cannot keep the old solution
				-- at the Dirichlet boundaries. Thus, it does not work correctly at the
				-- Dirichlet boundary where the boundary condition is specified not explicitely
				-- by a UserData object or LUA function but should be kept as in the initial condition.
				timeDisc:assemble_jacobian(A, u, gl)
				timeDisc:assemble_rhs(b, gl)
				linSolver:init(A, u)
				assembled_dt = currdt
			else
				timeDisc:prepare_step(solTimeSeries, currdt)
				timeDisc:assemble_rhs(b, gl)
			end
			
			-- apply linear solver
			if linSolver:apply(u,b) == false then
				if retValAtSolver ~= nil then return retValAtSolver end
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
				if retValAtMinStepSize ~= nil
				then return retValAtMinStepSize
				else test.require(false, "Time Solver failed.")
				end
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
		
		-- Post processing.
		if postProcess ~= nil then
      local pp_res = postProcess(u, step, time)
      if type(pp_res) == "boolean" and pp_res == false then -- i.e. not nil, not something else, but "false"!
          write("\n++++++ postProcess of the time step failed. ")
      end
	end -- main loop

	end -- legacy script
			
	if verbose then print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++") end
		
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
	
	if retValAtOK ~= nil then return retValAtOK end
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
	local tol = adaptiveStepInfo["TOLERANCE"] or 1e-3	-- tolerance
	local red = adaptiveStepInfo["REDUCTION"] or 0.5 	-- reduction
	local inc_fac = adaptiveStepInfo["INCREASE"] or 1.5	-- increase of time step
	local safety_fac = adaptiveStepInfo["SAFETY"] or 0.8   	-- safety factor
	local errorEst = adaptiveStepInfo["ESTIMATOR"]

	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end
	
	-- if red == nil then red = 0.5 end   -- reduction of time step
	-- if inc_fac == nil then inc_fac = 1.5 end   
	
	-- if errorEst == nil then errorEst = Norm2Estimator() end
	-- if tol == nil then 
	-- if safety_fac == nil then safety_fac = 0.8 end   
	
	
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


			local dtold = currdt
			
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
				-- bSuccess = true;  -- KAUST
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
					eps = timex:get_error_estimate(1)
		
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
				print ("EULEX-DT:\t".. dtEst .."\t"..inc_fac*currdt.."\t"..maxStepSize);
				currdt = math.min(dtEst, inc_fac*currdt, maxStepSize)
				

				if (eps <= tol) then 
					-- accept
					print ("EULEX-ACCEPTING:\t".. time .."\t"..dtold.."\t"..currdt.."\tq=\t2");
					bAcceptStep = true;
					bSuccess =true;
		 			
				else
	    			-- discard
	    		
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
			if bSuccess == false then 
			   print ("EULEX-REJECTING:\t".. time .."\t"..dtold.."\t"..currdt.."\tq=\t2");
			   break 
			end
				
			-- update new time
			time = timeDisc:future_time()
			nlsteps = nlsteps + newtonSolver:num_newton_steps() 	 
			
			-- push oldest solutions with new values to front, and 
			-- pop oldest sol pointer from end		
			oldestSol = solTimeSeries:oldest()
			VecScaleAssign(oldestSol, 1.0, u2)
			solTimeSeries:push_discard_oldest(oldestSol, time)
			
			if (doControl) then
				-- do the same for second 
				oldestSol = solTimeSeries2:oldest()
				VecScaleAssign(oldestSol, 1.0, u2)
				solTimeSeries2:push_discard_oldest(oldestSol, time)
			end
				
			if not (bFinishTimeStep == nil) and (bFinishTimeStep) then 
				timeDisc:finish_step_elem(solTimeSeries, u:grid_level()) 
			end
				
		end	-- while bSuccess == false -- aka 'try time step'
		
		-- step was executed successfully

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
	limexDesc,
	--reductionFactor,
	--tol,
	bFinishTimeStep,
	postProcess,
	startTSNo,
	endTSNo)


	-- read adaptive stuff
	local tol = limexDesc["TOLERANCE"] or  1e-3 
	local red = limexDesc["REDUCTION"] or 0.5           -- reduction of time step
	local inc_fac = limexDesc["INCREASE"] or 1.5        -- increase of time step
	local safety_fac = limexDesc["SAFETY"] or 0.8       -- safety factor
	local errorEst = limexDesc["ESTIMATOR"]
  
	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end

	if errorEst == nil then 
		print "WARNING: Error estimator not set. Default is euclidean norm! "
		errorEst = Norm2ErrorEst() 
	end


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
					eps = timex:get_error_estimate(1)
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


function util.SolveNonlinearProblemLimex(
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
  adaptiveDesc,
  postProcess)


  -- read adaptive stuff
  local inc_fac = adaptiveDesc["INCREASE"] or 1.5        -- increase of time step
 
  -- check parameters
  if filename == nil then filename = "sol" end
  if minStepSize == nil then minStepSize = maxStepSize end


  -- Check input parameters.
  if u == nil or domainDisc == nil or newtonSolver == nil
    or startTime == nil or endTime == nil or maxStepSize == nil then
    print("Wrong usage found. Please specify parameters as below:")
    
    if (u == nil) then print ("Did not find u!"); end;
    if (domainDisc == nil) then print ("Did not find domainDisc!"); end;

    if (startTime == nil) then print ("Did not find endTime!"); end;
    if (endTime == nil) then print ("Did not find endTime!"); end;
    if (maxStepSize == nil) then print ("Did not find maxStepSize!"); end;
    --util.PrintUsageOfSolveTimeProblem()
    exit()
  end


  print ("maxStepSize ="..maxStepSize)
  print ("minStepSize ="..minStepSize)

  print ("startTime ="..startTime)
  print ("endTime ="..endTime)
  
  
  -- Create LIMEX descriptor
  local limexDesc = {

        nstages = adaptiveDesc["STAGES"] or 2,
        steps = {1,2,3,4,5,6},
        nthreads = 1, 
        tol = adaptiveDesc["TOLERANCE"] or  1e-3,
        rhoSafetyOPT = adaptiveDesc["SAFETY"] or 0.25,

        dt = dt,
        dtmin = minStepSize,
        dtmax = maxStepSize,
        dtred = adaptiveDesc["REDUCTION"] or 0.5,  -- reduction of time step

        -- set disc & solver
        domainDisc= domainDisc,
        nonlinSolver = newtonSolver,   
        -- makeConsistent = true,

        matrixCache = true, -- or true,
        -- costStrategyOPT = time.limexDesc.costStrategyOPT,
        debugOPT = 5,

       -- dampScheideggerOPT = time.limexDesc.dampScheideggerOPT or 1.0,
       -- partialVeloMaskOPT = time.limexDesc.partialVeloMaskOPT or 0,
     }
  -- Create LIMEX object
  local limex = util.limex.CreateIntegrator(limexDesc)
      
   limex:set_time_step(limexDesc.dt)
   limex:set_dt_min(limexDesc.dtmin)
   limex:set_dt_max(limexDesc.dtmax)
   limex:set_reduction_factor(limexDesc.dtred)
  
    if (adaptiveDesc["DEBUG"]) then 
      --limex:set_debug(adaptiveDesc["DEBUG"])
      limex:set_debug_for_timestepper(adaptiveDesc["DEBUG"])
    end
  -- Register LUA callback.
  if type(postProcess) == "function" then 
    -- a) LUA functions
    local luaobserver = LuaCallbackObserver()
    
     function __util_LimexLuaCallbackPost(step, t, currdt) 
            local sol=luaobserver:get_current_solution()
            print(postProcess)
            postProcess(sol, step, t, currdt)
            return 1
      end
      luaobserver:set_callback("__util_LimexLuaCallbackPost") 
    limex:attach_observer(luaobserver)
   end
   
   -- Register VTK output callback.
   if type(out) == "userdata" then
    -- b) VTK output
    limex:attach_observer(VTKOutputObserver(filename, out))
   end
   
   
   
   local limexErrorEst 
   limexErrorEst = CompositeGridFunctionEstimator() 
   
   if (type(adaptiveDesc["SPACES"])=="table") then
    for i, _spacei in ipairs(adaptiveDesc["SPACES"]) do 
      print(_spacei)
      limexErrorEst:add(_spacei)
      
    end
   end -- table
  
 
  -- limex:set_space(limexErrorSpace)


   print(limexErrorEst:config_string())
   limex:add_error_estimator(limexErrorEst)
   
   -- Solve problem
    print(">> Solve using LIMEX...")
    
    -- Replace convergence check.
    local limexConvCheck = ConvCheck()
    limexConvCheck:set_maximum_steps(1)
    limexConvCheck:set_minimum_defect(1e-12)
    limexConvCheck:set_reduction(1e-9)
    limexConvCheck:set_verbose(true)
    limexConvCheck:set_supress_unsuccessful(true)
    
    newtonSolver:set_convergence_check(limexConvCheck) 
    newtonSolver:disable_line_search()
    
    -- Execute solver
    local sw = CuckooClock()
    sw:tic()
    print(newtonSolver:config_string())
    limex:apply(u, endTime, u,  startTime)
    print ("CDELTA="..sw:toc())
  return 
end

--[[!
\}
]]--
