-- Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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


function ug_error(msg, otherMsg, backtraceSkipLevel)
	print("\n")
	print("[->                             ERROR:")
	print("==============================================================================\n!")		
	if msg ~= nil then print(msg) end
	print("!\n==============================================================================")
	if otherMsg ~= nil then print(otherMsg) end
	local f, l = test.getSourceAndLine()
	print("     File:      "..f)
	print("     Line:      "..l)
	if msg ~= nil then print("     Message:   "..msg) end
	print("LUA BACKTRACE:")
	if backtraceSkipLevel == nil then
		backtraceSkipLevel = 2
	end
	DebugBacktrace(backtraceSkipLevel)		
	print("==============================================================================")
	print("                                 ERROR                                     <-]\n")
	
	--assert(false)
	exit()
end

--! use it like ug_assert(numPreRefs <= numRefs, "It must be chosen: numPreRefs <= numRefs")
--! @param condition the condition to assert
--! @param msg (optional) message to be printed if condition is not fulfilled
function ug_assert(condition, msg)
	if condition then
		return
	else
		ug_error(msg, "ASSERTION FAILED:")		
	end
end

function ug_warning(t)
	print(t)
	err_log(t)
end

function ug_cond_warning(condition, text)
	if condition then
		ug_warning(text)
	end
end



--! @param backtraceLevel the number of levels to go up
--! for backtraceLevel = 0, it returns file and line of the
--! calling function. for  backtraceLevel = 1 the
--! file and line of the caller of the calling function and so on.
function util.GetLUAFileAndLine (backtraceLevel)
	local level = 2+backtraceLevel
	local info = debug.getinfo(level, "Sl")
	if not info then return "" end
	if info.what == "C" then   -- is a C function?
		return "C function"
	else
		return string.format("[%s]:%d", info.short_src, info.currentline)
	end
end

--! Creates a grid function debug writer for the utilities. The function reads
--! the settings from the table util.debug. If this is a boolean variable then no
--! special settings are applied. If this variable is undefined the no debug writer
--! is created.
--! @param approxSpace the approximation space for the debug writer
function util.CreateGridFuncDebugWriter (approxSpace)
	if util.debug then
		if type (util.debug) ~= "boolean" and type (util.debug) ~= "table" then
			ug_error ("util.debug should be either boolean or a table.")
		end
		if util.debug == false then
			return nil
		end
		if approxSpace == nil then
			ug_error ("No approximation space specified for the debug writer")
		end
		
		util.debug_writer = GridFunctionDebugWriter(approxSpace)
		
		if type (util.debug) ~= "table" then
			local vtk = true
			local conn_viewer = false
			if debug.vtk ~= nil then vtk = debug.vtk end
			if debug.conn_viewer ~= nil then conn_viewer = debug.conn_viewer end
			util.debug_writer:set_vtk_output(vtk)
			util.debug_writer:set_conn_viewer_output(conn_viewer)
		end
		
		if util.debug_dir ~= nil and type (util.debug_dir) == "string" then
			CreateDirectory(util.debug_dir)
			util.debug_writer:set_base_dir(util.debug_dir)
		end
		
	end
	return nil
end
