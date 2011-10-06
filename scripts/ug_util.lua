-- Create util namespace
util = util or {}

ug_load_script("util/partition_maps.lua")

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

	return LoadDomain(domain, tname)
end

-- loads a GridObject.
-- If the file can not be found, the method tries to find it in ugs data path.
function util.LoadGridObject(gridObj, gridName)
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

	return LoadGridObject(gridObj, tname)
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

-- creates Neumann Boundary
function util.CreateNeumannBoundary(approxSpace, subsets)
	local neumannDisc = FV1NeumannBoundary()
	neumannDisc:set_approximation_space(approxSpace)
	neumannDisc:set_subsets(subsets)
	return neumannDisc
end

-- creates Dirichlet Boundary
function util.CreateDirichletBoundary(approxSpace)
	local dirichlet = DirichletBND()
	dirichlet:set_approximation_space(approxSpace)
	return dirichlet
end

-- creates Inner Boundary
function util.CreateInnerBoundary(approxSpace, functions, subsets)
	local innerDisc = FV1InnerBoundary()
	innerDisc:set_approximation_space(approxSpace)
	innerDisc:set_subsets(subsets)
	innerDisc:set_functions(functions)
	return innerDisc
end

-- creates FV1ConvDiff
function util.CreateFV1ConvDiff(approxSpace, functions, subsets)
	local elemDisc = ConvectionDiffusion()
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	elemDisc:set_disc_scheme("fv1")
	return elemDisc
end

-- creates FV1ConvDiff
function util.CreateFV1ThermohalineFlow(approxSpace, functions, subsets)
	local elemDisc = FV1ThermohalineFlow()
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

-- creates FV1ConstEq
function util.CreateFV1ConstEq(approxSpace, functions, subsets)
	local elemDisc = FV1ConstantEquation()
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

function util.CreateFE1ConvDiff(approxSpace, functions, subsets)
	local elemDisc = ConvectionDiffusion()
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	elemDisc:set_disc_scheme("fe")
	return elemDisc
end

-- creates FV1NavierStokes
function util.CreateFV1NavierStokes(approxSpace, functions, subsets)
	local elemDisc = FV1NavierStokes()	
	elemDisc:set_approximation_space(approxSpace)
	elemDisc:set_subsets(subsets)
	elemDisc:set_functions(functions)
	return elemDisc
end

-- create Geometric Multigrid
function util.CreateGeometricMultiGrid(approxSpace)
	local gmg = GeometricMultiGridPreconditioner()
	gmg:set_approximation_space(approxSpace)
	return gmg
end

-- create Prolongation / Restriction
function util.CreateP1Prolongation(approxSpace)
	local transfer = P1ProlongationOperator()
	transfer:set_approximation_space(approxSpace)
	return transfer
end

-- create Projection
function util.CreateP1Projection(approxSpace)
	local project = P1ProjectionOperator()
	project:set_approximation_space(approxSpace)
	return project
end

--------------------------------------------------------------------------------
-- User Data utils
--------------------------------------------------------------------------------

-- creates a Const User Matrix 2d 
function util.CreateConstUserMatrix2d(m00, m01, m10, m11)
	local mat = ConstUserMatrix2d()
	mat:set_entry(0, 0, m00)
	mat:set_entry(0, 1, m01)
	mat:set_entry(1, 0, m10)
	mat:set_entry(1, 1, m11)	
	return mat
end

-- creates a Const User Matrix 3d
function util.CreateConstUserMatrix3d(m00, m01, m02, m10, m11, m12, m20, m21, m22)
	local mat = ConstUserMatrix3d()
	mat:set_entry(0, 0, m00)
	mat:set_entry(0, 1, m01)
	mat:set_entry(0, 2, m02)
	mat:set_entry(1, 0, m10)
	mat:set_entry(1, 1, m11)
	mat:set_entry(1, 2, m12)
	mat:set_entry(2, 0, m20)
	mat:set_entry(2, 1, m21)
	mat:set_entry(2, 2, m22)	
	return mat
end

-- creates a Const User Vector 2d
function util.CreateConstUserVector2d(v0, v1)
	local vec = ConstUserVector2d()
	vec:set_entry(0, v0)
	vec:set_entry(1, v1)
	return vec
end

-- creates a Const User Vector 3d
function util.CreateConstUserVector3d(v0, v1, v2)
	local vec = ConstUserVector3d()
	vec:set_entry(0, v0)
	vec:set_entry(1, v1)
	vec:set_entry(2, v2)
	return vec
end

--------------------------------------------------------------------------------
-- Subset utils
--------------------------------------------------------------------------------

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
	local dom = Domain()
	
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
		if util.CheckSubsets(dom, neededSubsets) == false then 
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
function util.IsPowerOfTwo(n)
	local number compare = 1

	while (compare < n) do
		compare = compare*2
	end

	return compare==n
end

-- function returns true if the number is a natural number
function util.IsNaturalNumber(n)
	if n-math.floor(n) == 0 then
		return true
	else
		return false
	end
end

-- function to factorise number which has to be a power of 2 in two factors
-- which differ at most by a factor of 2 and returns both
-- (first the smaller one, then the larger one).
function util.FactorizeInPowersOfTwo(n)
	if not util.IsPowerOfTwo(n) then
		print("Number to factorise must be a power of 2. Aborting.")
		exit()
	end

	local number firstFactor = n
	local number secFactor = 1

	while (firstFactor > 2*secFactor) do
		firstFactor = firstFactor/2
		secFactor = secFactor*2
	end

	return secFactor, firstFactor
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
		util.PrintTableHelperIntend = util.PrintTableHelperIntend .. "  "
		
		for i,v in pairs(valuePar) do util.PrintTableHelper(i, v) end
		
		util.PrintTableHelperIntend = string.sub(util.PrintTableHelperIntend, 3)
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

-- adds writeln 
function writeln(...)
	write(...)
	write("\n")
end
