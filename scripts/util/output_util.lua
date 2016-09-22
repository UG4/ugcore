-- Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

util = util or {}


--------------------------------------------------------------------------------
-- lua script functions
--------------------------------------------------------------------------------

function util.TableToTextLongHelper(indexPar, valuePar)
	local str=""
	if type(valuePar) == "table" then
		str = str..util.PrintTableHelperIntend .. tostring(indexPar)  .. " = {\n"
		util.PrintTableHelperIntend = util.PrintTableHelperIntend .. "  "
		
		for i,v in pairs(valuePar) do str = str..util.TableToTextLongHelper(i, v) end
		
		util.PrintTableHelperIntend = string.sub(util.PrintTableHelperIntend, 3)
		str = str..util.PrintTableHelperIntend .. "}\n"
	else
		if type(valuePar) == "string" or type(valuePar) == "number" then
			str = str..util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. valuePar .."\n" 
		elseif type(valuePar) == "boolean" then
			str = str..util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. tostring(valuePar) .."\n"
		else
			str = str..util.PrintTableHelperIntend .. " " .. tostring(indexPar) .. " = " .. tostring(valuePar) .."\n"
		end
	end
	return str
end

function util.TableToTextLong(tablePar)
	util.PrintTableHelperIntend = ""
	return util.TableToTextLongHelper("", tablePar)
end

--! to print tables
function util.PrintTable(tablePar)
	print(util.TableToTextLong(tablePar))
end

function util.TableToText(var)
	local out= ""
	local i
	local v
	if type(var) == "table" then
		out = out.." {"
		
		local count = 0
		for _ in pairs(var) do count = count + 1 end
		if count == #var then count = 0 else count = 1 end
		
		local bfirst = true		
		for i,v in pairs(var) do
			if bfirst then bfirst = false else out=out .. ", " end
			if count == 1 then out=out .. "["..tostring(i).."] = " end 
			out=out..util.TableToText(v)
			 
		end
		
		out = out.. "}"
	else
		out = out ..tostring(var)		
	end
	return out
end




function util.Balance(DataToBeWrittenTable)
	-- check that table passed
	if type(DataToBeWrittenTable) ~= "table" then
		print("util.Balance: Data to be written must be given as a table")
		exit()
	end
	
	local verbose = true
	
	-- store outputs as local varible of closure 
	local Filenames = {}
	local VTK = {}
	local PosEvals = {}
	local Integrals = {}
	local Fluxes = {}
	local FuncValues = {}
	local PrintFreq = DataToBeWrittenTable.datafreq or 1
	
	----------------------------------
	-- loop "Data-to-be-written"-Table
	----------------------------------
	for _, DataSet in ipairs(DataToBeWrittenTable) do

		-- get file to write to	
		local file = nil
		local defFileCnt = 1
		if type(DataSet.file) == "string" then file = DataSet.file
		else
			if type(DataSet.file) ~= "nil" then print("util.Balance: file must be string"); exit(); end
			file = "file"..defFileCnt; defFileCnt = defFileCnt + 1;
		end 

		if table.contains(Filenames, file) then
			print("util.Balance: Filename '"..file.."' already used. Choose different name.")
			exit()															
		end
		table.insert(Filenames, file)
		
		----------------------------------
		-- check for vtk-data
		----------------------------------
		if DataSet.vtk ~= nil then
			local vtkOut = VTKOutput()
			local vtk = DataSet.vtk
			
			-- binary or ascii
			if DataSet.binary ~= nil then
				vtkOut:set_binary (DataSet.binary)
			end
			
			-- case: single sting passed
			if type(vtk) == "string" then
				vtkOut:select(vtk, vtk)
			end
			
			-- case: table of data passed
			if type(vtk) == "table" then
	
				local defDataCnt = 1
				for _, vtkItem in ipairs(vtk) do
				
					if type(vtkItem) ~= "table" and type(vtkItem) ~= "string" and type(vtkItem) ~= "userdata" then
						print("util.Balance: vtk must be: string or table of [table, string, userdata].")
						exit()													
					end
				
					-- get name of data in vtk-file
					local vtkName = nil
					if type(vtkItem) == "table" then vtkName = vtkItem[2] end
					if type(vtkItem) == "string" then vtkName = vtkItem end
					if type(vtkItem) == "userdata" then  end -- no default name
					if type(vtkName) ~= "string" then vtkName = "data"..defDataCnt; defDataCnt = defDataCnt+1; end
				
					-- get data to write to vtk-file
				
					local vtkData = nil
					if type(vtkItem) == "table" then vtkData = vtkItem[1] end
					if type(vtkItem) == "string" then vtkData = vtkItem end
					if type(vtkItem) == "userdata" then vtkData = vtkItem end
					if type(vtkData) ~= "string" and type(vtkData) ~= "userdata" then
						print("util.Balance: vtk data must be: string or userdata.")
						exit()													
					end	
					
					-- select vtk data-item
					vtkOut:select(vtkData, vtkName)

				end -- end table loop
			end -- end Parsing Vtk-Data-Table
						
			-- append to vtk datas
			table.insert(VTK, {vtkOut, file, DataSet.subsets})
			
		end -- end Parsing Vtk-Data
	
	
		----------------------------------
		-- check for position data
		----------------------------------
		if DataSet.value ~= nil then
			-- create position evaluation
			local PosEval = {}
			
			-- get data to be evaluated
			PosEval.value = DataSet.value
			if type(PosEval.value) ~= "userdata" then
				print("util.Balance: PointData: value must be userdata.")
				exit()													
			end
			
			-- get positions
			if type(DataSet.point) ~= "table" then
				print("util.Balance: PointData: no point(s) passed.")
				exit()													
			end
			
			local function CheckPointSet(set)
				for _, coord in ipairs(set) do
					if type(coord) ~= "number" then
						print("util.Balance: PointData: coordinate ".. coord .." is not a number")
						exit()																		
					end
				end
			end	
			
			PosEval.points = {}
			if type(DataSet.point[1]) == "table" then
				for _, point in ipairs(DataSet.point) do
					CheckPointSet(point)
					table.insert(PosEval.points, point)		
				end	
			else
				CheckPointSet(DataSet.point)
				table.insert(PosEval.points, DataSet.point)			
			end
		
			-- store file name
			PosEval.file = file

			-- separator 
			PosEval.sep = DataSet.sep
			if type(PosEval.sep) ~= "string" then
				PosEval.sep = "\t"
			end

			-- create file
			local thefile = io.open (file, "w+")
			thefile:write("# Position Evaluating file\n")
			thefile:write("time"..PosEval.sep)
			for _, point in ipairs(PosEval.points) do
				thefile:write("{")
				for i, coord in ipairs(point) do
					if i > 1 then thefile:write(",") end
					thefile:write(coord)
				end
				thefile:write("}",PosEval.sep)
			end						
			thefile:write("\n")
			io.close(thefile)
			
			-- append to pos-eval datas
			table.insert(PosEvals, PosEval)
		end

		----------------------------------
		-- check for integral data
		----------------------------------
		if DataSet.integral ~= nil then
			-- create integral data
			local IntegralData = {}

			-- get data to be evaluated
			IntegralData.data = DataSet.integral
			if type(IntegralData.data) ~= "userdata" then
				print("util.Balance: Integral: value must be userdata.")
				exit()													
			end
			
			-- store file name
			IntegralData.file = file
			
			-- separator 
			IntegralData.sep = DataSet.sep
			if type(IntegralData.sep) ~= "string" then
				IntegralData.sep = "\t"
			end
			
			-- subsets
			IntegralData.subsets = DataSet.subsets

			-- create file
			local thefile = io.open (file, "w+")
			thefile:write("# Integral data file\n")
			thefile:write("time"..IntegralData.sep.."value\n")
			io.close(thefile)
			
			-- append to integral datas
			table.insert(Integrals, IntegralData)			
		end

		----------------------------------
		-- check for flux data
		----------------------------------
		if DataSet.flux ~= nil then
			-- create integral data
			local FluxData = {}

			-- get data to be evaluated
			FluxData.data = DataSet.flux
			if type(FluxData.data) ~= "userdata" then
				print("util.Balance: Flux: value must be userdata.")
				exit()													
			end

			-- get boundary
			FluxData.boundary = DataSet.boundary
			if type(FluxData.boundary) ~= "string" then
				print("util.Balance: Flux: boundary must be specified.")
				exit()																
			end

			-- get boundary
			FluxData.inner = DataSet.inner
			if type(FluxData.inner) ~= "string" then
				print("util.Balance: Flux: inner must be specified.")
				exit()																
			end
			
			-- store file name
			FluxData.file = file

			-- separator 
			FluxData.sep = DataSet.sep
			if type(FluxData.sep) ~= "string" then
				FluxData.sep = "\t"
			end

			-- create file
			local thefile = io.open (file, "w+")
			thefile:write("# Flux data file\n")
			thefile:write("time"..FluxData.sep.."value\n")
			io.close(thefile)
			
			-- append to integral datas
			table.insert(Fluxes, FluxData)			
		end

		----------------------------------
		-- check for values of user-defined functions data
		----------------------------------
		if DataSet.func_data ~= nil then
			-- create a container for the object
			local FuncValue = {}
			FuncValue.func = DataSet.func_data
			FuncValue.data = DataSet.data
			FuncValue.file = file

			-- separator 
			FuncValue.sep = DataSet.sep
			if type(FuncValue.sep) ~= "string" then
				FuncValue.sep = "\t"
			end
			
			-- a comment to write at the beginning of the file
			if DataSet.comment ~= nil then
				-- create the file
				local thefile = io.open (file, "w+")
				thefile:write("# "..DataSet.comment.."\n")
				io.close(thefile)
			else
				-- clean the file
				local thefile = io.open (file, "w+")
				io.close(thefile)
			end
			
			-- append to the other objects
			table.insert(FuncValues, FuncValue)
		end
		
	end -- end DataSet loop
	
	-------------------------------------------------------------------
	-- the function being called, when data writer is invoked
	-------------------------------------------------------------------
	return function(u, step, time)
		
		if (step < 0) or (math.fmod(step, PrintFreq) == 0) then
		if verbose then print(" ******** Start Balancing ********") end		
		
		-- write VTK datas
		for _, vtkData in ipairs(VTK) do
			local vtkOut = vtkData[1]
			local file = vtkData[2]
			local subsets = vtkData[3]

			if verbose then write(" * Write VTK-data to '"..file.."' ... ") end
			
			if type(subsets) == "string" then
				vtkOut:print_subsets(file, u, subsets, step, time) 
				vtkOut:write_time_pvd(file, u) 				
			else
				vtkOut:print(file, u, step, time) 
				vtkOut:write_time_pvd(file, u) 
			end
			
			if verbose then print("done.") end
		end

		-- evaluate Pos-Datas
		for _, PosEval in ipairs(PosEvals) do
			local points = PosEval.points
			local value = PosEval.value
			local filename = PosEval.file
			local sep = PosEval.sep

			if verbose then print(" * Write PosData to '"..filename.."' ... ") end
		
			local file = io.open (filename, "a")
			file:write(time)
			file:write(sep)
			for _, point in ipairs(points) do
				local val = value:evaluate(point)
				--print("Point data is: "..val)
				file:write(val, sep)
			end
			file:write("\n")
			io.close(file)
			
			if verbose then print("done.") end
		end

		-- evaluate Integrals
		for _, IntegralData in ipairs(Integrals) do
			local data = IntegralData.data
			local filename = IntegralData.file
			local quadOrder = 2
			local sep = IntegralData.sep
			local subsets = IntegralData.subsets
	
			if verbose then write(" * Write Integral to '"..filename.."' ... ") end
		
			local file = io.open (filename, "a")

			local val = nil
			if type(subsets) == "string" then
				val = Integral(data, u, subsets, time)
			else
				val = Integral(data, u, time)
			end
			print("Integral is: "..val)
			file:write(time)
			file:write(sep)
			file:write(val)
			file:write("\n")
			io.close(file)
			
			if verbose then print("done.") end
		end

		-- evaluate Fluxes
		for _, FluxData in ipairs(Fluxes) do
			local data = FluxData.data
			local filename = FluxData.file
			local boundary = FluxData.boundary
			local inner = FluxData.inner
			local quadOrder = 2
			local sep = FluxData.sep
	
			if verbose then write(" * Write Flux to '"..filename.."' ... ") end
		
			local file = io.open (filename, "a")

			local val = IntegralNormalComponentOnManifold(data, u, boundary, inner)
			print("Flux is: "..val)
			file:write(time)
			file:write(sep)
			file:write(val)
			file:write("\n")
			io.close(file)
			
			if verbose then print("done.") end
		end
		
		-- evaluate the user-defined function data
		for _, FuncValue in ipairs(FuncValues) do
			local filename = FuncValue.file
			local func = FuncValue.func
			local data = FuncValue.data
			local sep = FuncValue.sep
	
			if verbose then write(" * Write Function Values to '"..filename.."' ... ") end
		
			local file = io.open (filename, "a")
			file:write(time)
			file:write(sep)
			file:write(func(data))
			file:write("\n")
			io.close(file)
			
			if verbose then print("done.") end
		end
		
		if verbose then print(" ******** End   Balancing ********")	end		
		end	
	end
end



--- adding tostring for booleans
_tostring = _tostring or tostring
function tostring(Val)
	if type(Val) == "table" then
   		return util.TableToTextLong(Val)
   	elseif type(Val) == "boolean" then
   		if Val then
   			return "true"
   		else
   			return "false"
   		end
   	else
   		return _tostring(Val)
   	end
end

if print_all == nil then
	function print_all(...)
		local la = GetLogAssistant()
		local opp = la:get_output_process()
		la:set_output_process(-1)
		print(unpack(arg))
		la:set_output_process(opp) 
	end
end    


