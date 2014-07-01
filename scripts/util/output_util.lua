util = util or {}
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
		
	end -- end DataSet loop
	
	-------------------------------------------------------------------
	-- the function being called, when data writer is invoked
	-------------------------------------------------------------------
	return function(u, step, time)
		if verbose then print(" ******** Start Balancing ********") end		
		
		-- write VTK datas
		for _, vtkData in ipairs(VTK) do
			local vtkOut = vtkData[1]
			local file = vtkData[2]
			local subsets = vtkData[3]

			if verbose then write(" * Write VTK-data to '"..file.."' ... ") end

			if type(subsets) == "string" then
				vtkOut:print(file, u, step, time) 
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
		end
		
		if verbose then print(" ******** End   Balancing ********")	end			
	end
end
