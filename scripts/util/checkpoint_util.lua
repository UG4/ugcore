-- Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
-- Author: Martin Rupp
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
	local filename = name..id..".ug4vec"
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
	
	-- todo: use if file exists here, and return nil if not
	ug_load_script(name..".lua")
	--! LoadTheCheckpoint is a function defined in the the lua binding:
        --! @see ugcore/ugbase/bindings/lua/lua_serialization.cpp
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

--! util.WriteStateCheckpoint
--! 
--! @param gf_names a table of names of grid functions to save (e.g. {u=u, lsf=lsf})
--! @param id the id for the current solution (e.g. the current time)
--! @param myData additional Data in Form of a table. 
--!        Will be available as cp.myData when loading the checkpoint
--!        Can contain all LUA types (tables, arrays, strings, numbers etc.), but no userdata!
--!        Can also be nil.
--! @param name name to use when writing the checkpoint (may be nil). 
function util.WriteStateCheckpoint (gf_names, id, myData, name)
	if type(gf_names) ~= "table" then
		print("WriteStateCheckpoint: The first argument should be a table {...} of names")
		exit()
	end
	
	if name == nil then name = util.checkpoint.stdName end
	ug_assert(id ~= nil)
	
--	create the directory
	local dirname = name.."."..id..".ug4cp"
	print("Writing state to directory '"..dirname.."':")
	if ProcRank() == 0 then
		CreateDirectory(dirname)
	end
	SynchronizeProcesses()
	ChangeDirectory(dirname)
	
--	write the grid functions
	for gf_n, gf in pairs (gf_names) do
		if gf == nil then
			print("WriteStateCheckpoint: No '"..gf_n.."' found!")
			exit()
		end
		local filename = gf_n..".ug4vec"
		print(" ... writing '"..gf_n.."'")
		SaveToFile(gf, filename)
	end
	
--	write the checkpoint data and return to the upper directory
	print(" ... writing environment")
	if ProcRank() == 0 then
		cpenv =
		{
			ugargc=ugargc, 
			ugargv=ugargv,
			commandline = util.GetCommandLine(),
		
			stdData=util.GetStdCheckpointData(),
		
			id=id,
			myData=myData
		}
		LuaWrite("env.lua", "cpenv")
	end
	ChangeDirectory("..")
	print(" ... done.\n")
end

--! util.ReadStateCheckpoint
--! 
--! @param gf_names a table of names of grid functions to save (e.g. {u=u, lsf=lsf})
--! @param id the id for the current solution (e.g. the current time)
--! @param name name used when writing the checkpoint (may be nil). 
function util.ReadStateCheckpoint(gf_names, id, name)
	if name == nil then name = util.checkpoint.stdName end
	
--	go to the directory and load the environment
	local dirname = name.."."..id..".ug4cp"
	print("Reading state from directory '"..dirname.."':")
	if not DirectoryExists(dirname) then
		print("---- Directory '"..dirname.."' does not exist! ----")
		exit()
	end
	ChangeDirectory(dirname)
	
--	load the lua data
	-- note: ug_load_script is also caring about distributing data to all cores. 
	ug_load_script("env.lua")
	local cp = LoadTheCheckpoint()
	
--	read the data
	for gf_n, gf in pairs (gf_names) do
		if gf == nil then
			print("ReadStateCheckpoint: No '"..gf_n.."' object created!")
			exit()
		end
		local filename = gf_n..".ug4vec"
		print(" ... reading '"..gf_n.."'")
		ReadFromFile(gf, filename)
	end
	
--	return to the parent directory
	ChangeDirectory("..")
	print(" ... done.")
	
	return cp
end

-- End of File
