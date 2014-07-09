--------------------------------------------------------------------------------
--[[!
-- \file scripts/util/checkpoint_util.lua
-- \author Martin Rupp 
-- \brief This provides restart/checkpointing functions
-- 
-- example writing a checkpoint (std is every 10s) 
	util.WriteCheckpointIntervallic(u, time, {time=time, step=step})
	
-- example reading a checkpoint
	local restart = util.HasParamOption("-restart") 
	if restart then
		cp = util.ReadCheckpoint(u)
		
		-- re-get our additional data
		time = cp.myData.time
		step = cp.myData.step
	end
	
-- Note: The used SaveToFile/ReadFromFile functions for u
-- assume that u has exactly the same structure as in the checkpoint.
-- That means approxSpace, grid and so on have to be the same
-- This will normally exclude adaptive calculations from restarting,
-- unless you add code which also reads/writes the adapted grid. 
-- See also util.GetStdCheckpointData for variables which are
-- automatically checked to be the same between runs. 
]]--
--------------------------------------------------------------------------------


util = util or {}

util.checkpoint = util.checkpoint or {}
util.checkpoint.stdName = "myCheckpoint"
util.checkpoint.stdIntervalMS = 10000

--! the standard checkpoint data
--! data in this table will be checked between loaded checkpoints
--! and the current run. e.g. numRefs and numCores should
--! be the same
--! note that values which are nil are NOT saved in the table
--! so it is OK to add values here which are not defined in every script 
function util.GetStdCheckpointData()
	return {
		numRefs=numRefs, 
		numPreRefs=numPreRefs,
		dim=dim,
		numCores=NumProcs()
	}
end

function util.nilstr(s)
	if s == nil then return "(nil)" else return s end
end

--! check that the stdData is same in file and
--! current run 
--! @sa util.GetStdCheckpointData
function util.CheckCheckpointData(cp)
	local err =""
	thisCP = util.GetStdCheckpointData()
	for i, v in pairs(cp.stdData) do
		if thisCP[i] ~= v then
			err = err.."   ERROR: "..i .. " is "..util.nilstr(v).. " in saved checkpoint, but "..util.nilstr(thisCP[i]) .." in current Run.\n"			
		end
	end
	
	if string.len(err) > 0 then
		print("\n  RESTART FAILED!")
		--print("--------------------------------------")
		--print("Saved checkpoint:")
		--print(cp)
		print("--------------------------------------")
		print(" current run's stdData:")
		print(thisCP)
		
		print("--------------------------------------")
		print("  RESTART FAILED! Reason:\n")
		
		print(err)
		
		print("--------------------------------------\n")
		ug_assert(false, err)
	end
end

--! util.WriteCheckpoint
--! 
--! @param u current solution to write (a GridFunction)
--! @param id the id for the current solution (like, the current time)
--! @param myData additional Data in Form of a table. 
--!        Will be available as cp.myData when loading the checkpoint
--!        Can contain all LUA types (tables, arrays, strings, numbers etc.), but no userdata!
--!        Can also be nil.
--! @param name name to use when writing the checkpoint. may be nil. 
function util.WriteCheckpoint(u, id, myData, name)
	if name == nil then name = util.checkpoint.stdName end
	ug_assert(id ~= nil)
	
	-- create a filename
	filename = name..id..".ug4vec"
	SaveToFile(u, filename)
	checkpoint =
	{
		ugargc=ugargc, 
		ugargv=ugargv,
		commandline = util.GetCommandLine(),
		
		stdData=util.GetStdCheckpointData(),
		
		lastId=id,
		lastFilename=filename,
		myData=myData
	}
	
	--print(checkpoint)
	if ProcRank() == 0 then
		LuaWrite(name..".lua", "checkpoint")
	end
	
	print("\nWrote Checkpoint "..name..".lua. id = "..id..", filename = "..filename.."\n")		
	
end

--! util.WriteCheckpointIntervallic
--! @sa util.WriteCheckpoint
--! use this function to write checkpoint data in time intervals
--! @param IntervalMS write checkpoint every IntervalMS milliseconds.
--!        can be null, then defaults to util.checkpoint.stdIntervalMS
function util.WriteCheckpointIntervallic(u, id, myData, name, IntervalMS)
	local timeMS = GetClockS()*1000
	
	if IntervalMS == nil then IntervalMS = util.checkpoint.stdIntervalMS end
	if util.checkpoint.lastCheckpointTimeMS == nil 
	 	or timeMS -util.checkpoint.lastCheckpointTimeMS > IntervalMS then
	 	
	 	util.WriteCheckpoint(u, id, myData, name)
	 	util.checkpoint.lastCheckpointTimeMS = timeMS
	end	
end

--! util.ReadCheckpoint
--! 
--! @param u the GridFunction to read into
--! @param name name to use when writing the checkpoint. may be nil.
--! @return the checkpoint data cp. (see util.WriteCheckpoint)
--!  myData will be saved  cp.myData.
function util.ReadCheckpoint(u, name)
	if name == nil then name = util.checkpoint.stdName end
	
	-- load the lua data
	-- note: ug_load_script is also caring about distributing
	-- data to all cores. 
	ug_load_script(name..".lua")
	-- the script defines a function named "LoadTheCheckpoint"
	-- call it
	local cp = LoadTheCheckpoint()
	
	-- print the checkpoint
	print("Loading Checkpoint "..name..".lua :")
	print(cp)
	
	-- check that the stdData is same in file and
	-- current run (see util.GetStdCheckpointData)
	util.CheckCheckpointData(cp)

	ReadFromFile(u, cp.lastFilename)
	
	util.checkpoint.lastCheckpointTimeMS = GetClockS()*1000
	return cp
end





function util.SolveLinearTimeProblemRestart(
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
	
	
	-- print newtonSolver setup	
	print("SolveLinearTimeProblem, Linear Solver setup:")
	print(linSolver:config_string())
	
	-- start
	local time = startTime
	local step = 0
	
	-- write start solution
	print(">> Writing start values")
	if not (out==nil) then out:print(filename, u, step, time) end
	
	-- store grid function in vector of  old solutions
	local solTimeSeries = SolutionTimeSeries()
	
	----------------------------------------------------------
	--- Read Checkpoint if necessary
	local restart = util.HasParamOption("-restart") 
	if restart then
		cp = util.ReadCheckpoint(u)
		time = cp.myData.time
		step = cp.myData.step
	end
	----------------------------------------------------------
			
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
			
		print("++++++ TIMESTEP "..step.." END   (current time: " .. time .. ") ++++++");
		
		SaveVectorForConnectionViewer(u, "vec"..time..".vec")
		
		----------------------------------------------------------
		--- Write Checkpoint if necessary
		util.WriteCheckpointIntervallic(u, time, {time=time, step=step})
		----------------------------------------------------------
	end
	
	if not (out==nil) then out:write_time_pvd(filename, u) end
end
 