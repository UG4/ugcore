-- Create util namespace
util = util or {}

-- Creates a domain for the given dimension. 1, 2, 3 are supported.
function util.CreateDomain(dim)
	if dim == 1 then
		return Domain1d()
	elseif dim == 2 then
		return Domain2d()
	elseif dim == 3 then
		return Domain3d()
	end
	return nil
end

-- returns the standard path at which grids are stored
function util.GetGridPath()
	return ug_get_data_path().."/grids/"
end

-- loads a domain.
-- If the file can not be found, the method tries to find it in ugs data path.
function util.LoadDomain(domain, gridName)
	local dim = domain:get_dim()

	local tname = gridName
--	first try the original name.
	if ug_file_exists(tname) == false then
	--	now try relative to the current path.
		tname = ug_get_current_path() .. "/" .. gridName
		if ug_file_exists(tname) == false then
		--	finally check the default grid path.
			tname = ug_get_data_path().."/grids/" .. gridName
			if ug_file_exists(tname) == false then
				return false
			end
		end
	end

	return LoadDomain(domain, tname, 0)
end


function util.GlobalRefineParallelDomain(domain)
	local dim = domain:get_dim()
	if dim == 1 then
		return GlobalRefineParallelDomain1d(domain)
	elseif dim == 2 then
		return GlobalRefineParallelDomain2d(domain)
	elseif dim == 3 then
		return GlobalRefineParallelDomain3d(domain)
	end
	return false
end

	
-- Creates an Approxmiation space of correct dimesion
function util.CreateApproximationSpace(domain)
	local approxSpace = nil
	if dim == 1 then
		approxSpace = ApproximationSpace1d()
	elseif dim == 2 then
		approxSpace = ApproximationSpace2d()
	elseif dim == 3 then
		approxSpace = ApproximationSpace3d()
	end
	approxSpace:assign_domain(domain)
	
	return approxSpace
end

-- creates Neumann Boundary
function util.CreateNeumannBoundary(approxSpace, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local neumannDisc
	if dim == 1 then
		neumannDisc = FV1NeumannBoundary1d()
	elseif dim == 2 then
		neumannDisc = FV1NeumannBoundary2d()
	elseif dim == 3 then
		neumannDisc = FV1NeumannBoundary3d()
	else
		return nil
	end

	neumannDisc:set_approximation_space(approxSpace)
	neumannDisc:set_subsets(subsets)
	return neumannDisc
end

-- creates Dirichlet Boundary
function util.CreateDirichletBoundary(approxSpace)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local dirichlet
	if dim == 1 then
		dirichlet = DirichletBND1d()
	elseif dim == 2 then
		dirichlet = DirichletBND2d()
	elseif dim == 3 then
		dirichlet = DirichletBND3d()
	else
		return nil
	end

	dirichlet:set_approximation_space(approxSpace)
	return dirichlet
end

-- creates Inner Boundary
function util.CreateInnerBoundary(approxSpace, functions, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local innerDisc
	if dim == 1 then
		innerDisc = FV1InnerBoundary1d()
	elseif dim == 2 then
		innerDisc = FV1InnerBoundary2d()
	elseif dim == 3 then
		innerDisc = FV1InnerBoundary3d()
	else
		return nil
	end

	innerDisc:set_approximation_space(approxSpace)
	innerDisc:set_subsets(subsets)
	innerDisc:set_functions(functions)
	return innerDisc
end

-- creates FV1ConvDiff
function util.CreateFV1ConvDiff(approxSpace, functions, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local elemDisc
	if dim == 1 then
		elemDisc = FV1ConvectionDiffusion1d()
	elseif dim == 2 then
		elemDisc = FV1ConvectionDiffusion2d()
	elseif dim == 3 then
		elemDisc = FV1ConvectionDiffusion3d()
	else
	return nil
	end
	
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

-- creates FV1ConvDiff
function util.CreateFV1ThermohalineFlow(approxSpace, functions, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local elemDisc
	if dim == 1 then
		elemDisc = FV1ThermohalineFlow1d()
	elseif dim == 2 then
		elemDisc = FV1ThermohalineFlow2d()
	elseif dim == 3 then
		elemDisc = FV1ThermohalineFlow3d()
	else
	return nil
	end
	
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

-- creates FV1ConstEq
function util.CreateFV1ConstEq(approxSpace, functions, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local elemDisc
	if dim == 1 then
		elemDisc = FV1ConstantEquation1d()
	elseif dim == 2 then
		elemDisc = FV1ConstantEquation2d()
	elseif dim == 3 then
		elemDisc = FV1ConstantEquation3d()
	else
	return nil
	end
	
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

function util.CreateFE1ConvDiff(approxSpace, functions, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local elemDisc
	if dim == 1 then
		elemDisc = FE1ConvectionDiffusion1d()
	elseif dim == 2 then
		elemDisc = FE1ConvectionDiffusion2d()
	elseif dim == 3 then
		elemDisc = FE1ConvectionDiffusion3d()
	else
		return nil
	end
	
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

function util.CreateUpwind(type, dim)
	local upwind
	
	if type == "no" then
		if 		dim == 1 then upwind = NoUpwind1d()
		elseif	dim == 2 then upwind = NoUpwind2d()
		elseif 	dim == 3 then upwind = NoUpwind3d()
		else return nil end
	elseif type == "full" then
		if 		dim == 1 then upwind = FullUpwind1d()
		elseif	dim == 2 then upwind = FullUpwind2d()
		elseif 	dim == 3 then upwind = FullUpwind3d()
		else return nil end
	elseif type == "part" then
		if 		dim == 1 then upwind = PartialUpwind1d()
		elseif	dim == 2 then upwind = PartialUpwind2d()
		elseif 	dim == 3 then upwind = PartialUpwind3d()
		else return nil end
	elseif type == "weighted" then
		if 		dim == 1 then upwind = WeightedUpwind1d()
		elseif	dim == 2 then upwind = WeightedUpwind2d()
		elseif 	dim == 3 then upwind = WeightedUpwind3d()
		else return nil end
	else
		print("Type of Upwind not recognized. Aborting.")
		exit();
	end
		
	return upwind
end

-- create Geometric Multigrid
function util.CreateGeometricMultiGrid(approxSpace)
	local dim = approxSpace:get_domain():get_dim()
	local gmg = nil;
	if dim == 1 then
		gmg = GeometricMultiGridPreconditioner1d()
	elseif dim == 2 then
		gmg = GeometricMultiGridPreconditioner2d()
	elseif dim == 3 then
		gmg = GeometricMultiGridPreconditioner3d()
	else
		gmg = nil
	end
	
	gmg:set_approximation_space(approxSpace)
	return gmg
end

-- create Prolongation / Restriction
function util.CreateP1Prolongation(approxSpace)
	local dim = approxSpace:get_domain():get_dim()
	local transfer = nil;
	if dim == 1 then
		transfer = P1ProlongationOperator1d()
	elseif dim == 2 then
		transfer = P1ProlongationOperator2d()
	elseif dim == 3 then
		transfer = P1ProlongationOperator3d()
	else
	transfer = nil
	end
	
	transfer:set_approximation_space(approxSpace)
	return transfer
end

-- create Projection
function util.CreateP1Projection(approxSpace)
	local dim = approxSpace:get_domain():get_dim()
	local project = nil;
	if dim == 1 then
		project = P1ProjectionOperator1d()
	elseif dim == 2 then
		project = P1ProjectionOperator2d()
	elseif dim == 3 then
		project = P1ProjectionOperator3d()
	else
	project = nil
	end
	
	project:set_approximation_space(approxSpace)
	return project
end


-- creates a Lua User Matrix using a lua function and returns the Provider
function util.CreateLuaUserMatrix(funcName, dim)
	local mat = nil
	if dim == 1 then
		mat = LuaUserMatrix1d()
	elseif dim == 2 then
		mat = LuaUserMatrix2d()
	elseif dim == 3 then
		mat = LuaUserMatrix3d()
	end
		
	mat:set_lua_callback(funcName)
	return mat
end

-- creates a Const User Matrix using a lua function and returns the Provider
function util.CreateConstDiagUserMatrix(diagVal, dim)
	local mat = nil
	if dim == 1 then
		mat = ConstUserMatrix1d()
	elseif dim == 2 then
		mat = ConstUserMatrix2d()
	elseif dim == 3 then
		mat = ConstUserMatrix3d()
	end

	mat:set_diag_tensor(diagVal)
	return mat
end

-- creates a Lua User Vector using a lua function and returns the Provider
function util.CreateLuaUserVector(funcName, dim)
	local vec = nil
	if dim == 1 then
		vec = LuaUserVector1d()
	elseif dim == 2 then
		vec = LuaUserVector2d()
	elseif dim == 3 then
		vec = LuaUserVector3d()
	end

	vec:set_lua_callback(funcName)
	return vec
end

-- creates a Const User Vector using a lua function and returns the Provider
function util.CreateConstUserVector(val, dim)
	local vec = nil
	if dim == 1 then
		vec = ConstUserVector1d()
	elseif dim == 2 then
		vec = ConstUserVector2d()
	elseif dim == 3 then
		vec = ConstUserVector3d()
	end
	
	vec:set_all_entries(val)
	return vec
end

-- creates a Lua User Number using a lua function and returns the Provider
function util.CreateLuaUserNumber(funcName, dim)
	local number = nil
	if dim == 1 then
		number = LuaUserNumber1d()
	elseif dim == 2 then
		number = LuaUserNumber2d()
	elseif dim == 3 then
		number = LuaUserNumber3d()
	end
	
	number:set_lua_callback(funcName)
	return number
end

-- creates a Const User Number using a lua function and returns the Provider
function util.CreateConstUserNumber(val, dim)
	local number = nil
	if dim == 1 then
		number = ConstUserNumber1d()
	elseif dim == 2 then
		number = ConstUserNumber2d()
	elseif dim == 3 then
		number = ConstUserNumber3d()
	end
	
	number:set(val)
	return number
end


-- creates a Lua Boundary Numver using a lua function and returns the Provider
function util.CreateLuaBoundaryNumber(funcName, dim)
	local number = nil
	if dim == 1 then
		number = LuaBoundaryNumber1d()
	elseif dim == 2 then
		number = LuaBoundaryNumber2d()
	elseif dim == 3 then
		number = LuaBoundaryNumber3d()
	end
	
	number:set_lua_callback(funcName)
	return number
end

-- creates a Const Boundary Number using a lua function and returns the Provider
function util.CreateConstBoundaryNumber(val, dim)
	local number = nil
	if dim == 1 then
		number = ConstBoundaryNumber1d()
	elseif dim == 2 then
		number = ConstBoundaryNumber2d()
	elseif dim == 3 then
		number = ConstBoundaryNumber3d()
	end
	
	number:set(val)
	return number
end


-- creates a GridFunctionDebugWriter
function util.CreateGridFunctionDebugWriter(dim)
	local writer = nil
	if dim == 1 then
		writer = GridFunctionDebugWriter1d()
	elseif dim == 2 then
		writer = GridFunctionDebugWriter2d()
	elseif dim == 3 then
		writer = GridFunctionDebugWriter3d()
	end

	return writer
end


-- creates a VTKWriter
function util.CreateVTKWriter(dim)
	local writer = nil
	if dim == 1 then
		writer = VTKOutput1d()
	elseif dim == 2 then
		writer = VTKOutput2d()
	elseif dim == 3 then
		writer = VTKOutput3d()
	end
	
	return writer
end

--------------------------------------------------------------------------------
-- Navier - Stokes utils
--------------------------------------------------------------------------------

-- creates FV1NavierStokes
function util.CreateFV1NavierStokes(approxSpace, functions, subsets)
	local domain = approxSpace:get_domain()
	local dim = domain:get_dim()
	local elemDisc
	if dim == 1 then
		elemDisc = FV1NavierStokes1d()
	elseif dim == 2 then
		elemDisc = FV1NavierStokes2d()
	elseif dim == 3 then
		elemDisc = FV1NavierStokes3d()
	else
	return nil
	end
	
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end


-- creates NavierStokesNoUpwind
function util.CreateNavierStokesNoUpwind(dim)
	local upwind
	if     dim == 1 then upwind = NavierStokesNoUpwind1d()
	elseif dim == 2 then upwind = NavierStokesNoUpwind2d()
	elseif dim == 3 then upwind = NavierStokesNoUpwind3d()
	else return nil end

	return upwind
end

-- creates NavierStokesFullUpwind
function util.CreateNavierStokesFullUpwind(dim)
	local upwind
	if     dim == 1 then upwind = NavierStokesFullUpwind1d()
	elseif dim == 2 then upwind = NavierStokesFullUpwind2d()
	elseif dim == 3 then upwind = NavierStokesFullUpwind3d()
	else return nil end

	return upwind
end

-- creates NavierStokesSkewedUpwind
function util.CreateNavierStokesSkewedUpwind(dim)
	local upwind
	if     dim == 1 then upwind = NavierStokesSkewedUpwind1d()
	elseif dim == 2 then upwind = NavierStokesSkewedUpwind2d()
	elseif dim == 3 then upwind = NavierStokesSkewedUpwind3d()
	else return nil end

	return upwind
end

-- creates NavierStokesLinearProfileSkewedUpwind
function util.CreateNavierStokesLinearProfileSkewedUpwind(dim)
	local upwind
	if     dim == 1 then upwind = NavierStokesLinearProfileSkewedUpwind1d()
	elseif dim == 2 then upwind = NavierStokesLinearProfileSkewedUpwind2d()
	elseif dim == 3 then upwind = NavierStokesLinearProfileSkewedUpwind3d()
	else return nil end

	return upwind
end

-- creates NavierStokesPositiveUpwind
function util.CreateNavierStokesPositiveUpwind(dim)
	local upwind
	if     dim == 1 then upwind = NavierStokesPositiveUpwind1d()
	elseif dim == 2 then upwind = NavierStokesPositiveUpwind2d()
	elseif dim == 3 then upwind = NavierStokesPositiveUpwind3d()
	else return nil end

	return upwind
end

-- creates NavierStokesFIELDSStabilization
function util.CreateNavierStokesFIELDSStabilization(dim)
	local stab
	if     dim == 1 then stab = NavierStokesFIELDSStabilization1d()
	elseif dim == 2 then stab = NavierStokesFIELDSStabilization2d()
	elseif dim == 3 then stab = NavierStokesFIELDSStabilization3d()
	else return nil end

	return stab
end

-- creates DarcyVelocityLinker
function util.CreateDarcyVelocityLinker(dim)
	local vel
	if     dim == 1 then vel = DarcyVelocityLinker1d()
	elseif dim == 2 then vel = DarcyVelocityLinker2d()
	elseif dim == 3 then vel = DarcyVelocityLinker3d()
	else return nil end

	return vel
end


--- util.CheckSubsets
-- checks if all required subsets are contained in the SubsetHandler
-- @param dom Domain
-- @param neededSubsets List of subsets the SubsetHandler must contain
-- @return true if all subsets are contained, false else
function util.CheckSubsets(dom, neededSubsets)
	sh = dom:get_subset_handler()
	for i, tval in ipairs(neededSubsets) do
		if sh:get_subset_index(tval) == -1 then
			print("Domain does not contain subset '"..tval.."'.")
			return false
		end
	end
	
	return true
end


function util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

	-- create Instance of a Domain
	local dom = util.CreateDomain(dim)
	
	-- load domain
	if util.LoadDomain(dom, gridName) == false then
	   print("Loading Domain failed. Aborting.")
	   exit()
	end
	
	-- create Refiner
	if numPreRefs > numRefs then
		print("numPreRefs must be smaller than numRefs. Aborting.");
		exit();
	end
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner = GlobalDomainRefiner(dom)
	
	-- Performing pre-refines
	for i=1,numPreRefs do
		refiner:refine()
	end
	
	-- Distribute the domain to all involved processes
	if DistributeDomain(dom) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	
	-- Perform post-refine
	for i=numPreRefs+1,numRefs do
		refiner:refine()
	end
	
	-- Now we loop all subsets an search for it in the SubsetHandler of the domain
	if neededSubsets ~= nil then
		if util.CheckSubsets(dom, neededSubsets) then 
			print("Wrong subsets detected.") 
		end
	end
	
	-- TODO: Delete the refiner
	
	-- return the created domain
	return dom
end

--------------------------------------------------------------------------------
-- some auxiliary functions
--------------------------------------------------------------------------------
-- function returns true if the number is a power of two
function util.isPowerOfTwo(n)
	local number compare = 1

	while (compare < n) do
		compare = compare*2
	end

	return compare==n
end

-- function returns true if the number is a natural number
function util.isNaturalNumber(n)
	if n-math.floor(n) == 0 then
		return true
	else
		return false
	end
end

--------------------------------------------------------------------------------
-- Command line functions
--------------------------------------------------------------------------------


--- util.GetParam
-- returns parameter in ugargv after ugargv[i] == name
-- @param name parameter in ugargv to search for
-- @param return_if_unavailable when parameter 'name' is not found, this will be returned
-- @return parameter in ugargv after ugargv[i] == name
function util.GetParam(name, return_if_unavailable)
	local i
	for i = 1, ugargc-1 do
		if ugargv[i] == name then
			return ugargv[i+1]
		end
	end
	return return_if_unavailable; 
end


--- util.GetParamNumber
-- use with CommandLine to get option, like -useAMG
-- if parameter is not a number, returns return_if_unavailable
-- @param name parameter in ugargv to search for
-- @return the number after the parameter 'name'
function util.GetParamNumber(name, return_if_unavailable)
	local param = util.GetParam(name, return_if_unavailable)
	local number = tonumber(param)
	if number == nil then
		print("WARNING: Parameter "..name.." is not a number, using "..return_if_unavailable.." instead\n") 
		return return_if_unavailable
	else
		return number
	end
end

--- util.HasParamOption
-- use with CommandLine to get option, like -useAMG
-- @param name option in argv to search for
-- @return true if option found, else false
function util.HasParamOption(name)
	for i = 1, ugargc do
		if ugargv[i] == name then
			return true
		end
	end
	return false 
end


--------------------------------------------------------------------------------
-- lua script functions
--------------------------------------------------------------------------------

function util.PrintTableHelper(indexPar, valuePar)
	if type(valuePar) == "table" then
		print(util.PrintTableHelperIntend .. tostring(indexPar)  .. " = {")
		util.PrintTableHelperIntend = util.PrintTableHelperIntend .. " "
		table.foreachi (valuePar, util.PrintTableHelper)
		util.PrintTableHelperIntend = string.sub(util.PrintTableHelperIntend, 2)
		print(util.PrintTableHelperIntend .. "}")
	else
		print(util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. valuePar)
	end
end

-- to print tables
function util.PrintTable(tablePar)
	util.PrintTableHelperIntend = ""
	util.PrintTableHelper("", tablePar)
end