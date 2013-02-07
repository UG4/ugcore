
function util.SolveNonlinearTimeProblem(
	u,
	timeDisc,
	newtonSolver,
	out,
	filename,
	startTime,
	endTime,
	maxStepSize,
	minStepSize,
	reductionFactor)

	-- check parameters
	if out == nil then out = VTKOutput() end
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end
	
	-- start
	local time = startTime
	local step = 0
	
	-- write start solution
	print(">> Writing start values")
	out:print(filename, u, step, time)
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(u:clone(), time)
	
	newtonSolver:init(AssembledOperator(timeDisc))
	
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
			print("++++++ Time step size: "..currdt);
			
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
				exit();
			end
		end
		
		-- update new time
		time = solTimeSeries:time(0) + currdt
		
		-- plot solution
		out:print(filename, u, step, time)
		
		-- get oldest solution
		local oldestSol = solTimeSeries:oldest()
	
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
	
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
	end
end



function util.SolveLinearTimeProblem(
	u,
	timeDisc,
	linSolver,
	out,
	filename,
	startTime,
	endTime,
	maxStepSize,
	minStepSize,
	reductionFactor)

	-- check parameters
	if out == nil then out = VTKOutput() end
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end
	
	-- start
	local time = startTime
	local step = 0
	
	-- write start solution
	print(">> Writing start values")
	out:print(filename, u, step, time)
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	solTimeSeries:push(u:clone(), time)

	-- matrix and vectors
	local A = AssembledLinearOperator(timeDisc)
	local b = u:clone()

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
			print("++++++ Time step size: "..currdt);

			-- reassemble matrix if necessary
			if not(currdt == assembled_dt) then 
				print("++++++ Assembling Matrix/Rhs for step size "..currdt); 
				timeDisc:prepare_step(solTimeSeries, currdt)
				timeDisc:assemble_linear(A, b)
				linSolver:init(A)
				assembled_dt = currdt
			else
				timeDisc:assemble_rhs(b)
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
				exit();
			end
		end
		
		-- update new time
		time = solTimeSeries:time(0) + currdt
		
		-- plot solution
		out:print(filename, u, step, time)
		
		-- get oldest solution
		local oldestSol = solTimeSeries:oldest()
	
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
	
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
	end
end