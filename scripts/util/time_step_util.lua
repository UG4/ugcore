
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
	
	else
		print("Time scheme '"..timeScheme.."' not found. Supported Schemes:")
		print("Theta, ImplEuler, ExplEuler, Crank-Nicolson, Alexander, FracStep, BDF")
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
	print("minStepSize      -- (optinal) minimal step sized used")	
	print("reductionFactor  -- (optinal) factor by which the step size is ")
	print("                    reduced, if the problem was not solved. ")
	print("                    Iterated until minStepSize is reached.")
end


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
	reductionFactor)

	if u == nil or domainDisc == nil or newtonSolver == nil or timeScheme == nil
		or startTime == nil or endTime == nil or maxStepSize == nil then
		print("Wrong usage found. Please specify parameters as below:")
		util.PrintUsageOfSolveTimeProblem()
		exit()
	end

	-- check parameters
	if filename == nil then filename = "sol" end
	if minStepSize == nil then minStepSize = maxStepSize end
	if reductionFactor == nil then reductionFactor = 0.5 end
	
	-- create time disc
	local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)
	
	-- start
	local time = startTime
	local step = 0
	
	-- write start solution
	print(">> Writing start values")
	if not (out==nil) then out:print(filename, u, step, time) end
	
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
					exit();
				end
				
				-- start over again if failed
				if bSuccess == false then break end
				
				-- update new time
				time = timeDisc:future_time()
				
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

				if timeDisc:num_stages() > 1 then
					print("      +++ STAGE " .. stage .. " END   ++++++")
				end
			end
		end

		-- save this solution if multistage
		if timeDisc:num_stages() > 1 then
			uOld = u
		end			
		
		-- plot solution
		if not (out==nil) then out:print(filename, u, step, time) end
	
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
	end
end



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
	reductionFactor)

	if u == nil or domainDisc == nil or linSolver == nil or timeScheme == nil
		or startTime == nil or endTime == nil or maxStepSize == nil then
		print("Wrong usage found. Please specify parameters as below:")
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
	
	-- start
	local time = startTime
	local step = 0
	
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
				exit();
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
			
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
	end
end