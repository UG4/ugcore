-- Creates a domain for the given dimension. 1, 2, 3 are supported.
function utilCreateDomain(dim)
	if dim == 1 then
		return Domain1d()
	elseif dim == 2 then
		return Domain2d()
	elseif dim == 3 then
		return Domain3d()
	end
	return nil
end

-- loads a domain. Automatically chooses right methods depending on the
-- domains dimension.
-- If the file can not be found, the method tries to find it in ugs data path.
function utilLoadDomain(domain, gridName)
	local dim = domain:get_dim()

	local tname = gridName
	if ug_file_exists(tname) == false then
		tname = ug_get_data_path().."/grids/" .. gridName
		if ug_file_exists(tname) == false then
			return false
		end
	end

	if dim == 1 then
		return LoadDomain1d(domain, tname)
	elseif dim == 2 then
		return LoadDomain2d(domain, tname)
	elseif dim == 3 then
		return LoadDomain3d(domain, tname)
	end
	return false
end

function utilSaveDomain(domain, gridName)
	local dim = domain:get_dim()
	if dim == 1 then
		return SaveDomain1d(domain, gridName)
	elseif dim == 2 then
		return SaveDomain2d(domain, gridName)
	elseif dim == 3 then
		return SaveDomain3d(domain, gridName)
	end
	return false
end

-- Creates anApproximationSpace of the dimension of the given domain.
-- The associated P1ConformFunctionPattern is created and assigned on the fly.
-- domain has to be either of type Domain1d, Domain2d or Domain3d
-- funcNames has to be an array containing strings
function utilCreateApproximationSpaceP1(domain, funcNames)
	local dim = domain:get_dim()
	-- create function pattern, add the functions and lock it
	local pattern = P1ConformFunctionPattern()
	pattern:set_subset_handler(domain:get_subset_handler())
--TODO iterate over the funcNames and add them
	AddP1Function(pattern, funcNames, dim)
	pattern:lock()

	-- create Approximation Space and assign the function pattern
	-- and the domain
	approxSpace, pattern = utilCreateApproximationSpace(domain, pattern)
	return approxSpace, pattern
end
	
-- Creates an Approxmiation space of correct dimesion and assigns the pattern
function utilCreateApproximationSpace(domain, pattern)
	local approxSpace = nil
	if dim == 1 then
		approxSpace = ApproximationSpace1d()
	elseif dim == 2 then
		approxSpace = ApproximationSpace2d()
	elseif dim == 3 then
		approxSpace = ApproximationSpace3d()
	end
	approxSpace:assign_domain(domain)
	approxSpace:assign_function_pattern(pattern)
	
	return approxSpace, pattern
end

-- creates a 2d dirichlet boundary and assigns the given parameters.
-- todo: specify types
function utilCreateDirichletBnd2d(domain, bndFct, fctInd)
	local dirichletBND = DirichletBND2d()
	dirichletBND:set_domain(domain)
	dirichletBND:set_dirichlet_function(bndFct)
	dirichletBND:set_function(fctInd)
	return dirichletBND
end

-- creates Neumann Boundary
function utilCreateNeumannBoundary(domain, userFunc, bndSubset)
	local dim = domain:get_dim()
	local neumannDisc
	if dim == 1 then
		neumannDisc = FV1NeumannBoundaryElemDisc1d()
	elseif dim == 2 then
		neumannDisc = FV1NeumannBoundaryElemDisc2d()
	elseif dim == 3 then
		neumannDisc = FV1NeumannBoundaryElemDisc3d()
	else
		return nil
	end

	neumannDisc:set_domain(dom)
	neumannDisc:set_bnd_cond(userFunc, bndSubset)
	return neumannDisc
end

-- creates FV1ConvDiffElemDisc
function utilCreateFV1ConvDiffElemDiscNoInit(domain)
	local dim = domain:get_dim()
	local elemDisc
	if dim == 1 then
		elemDisc = FV1ConvectionDiffusionElemDisc1d()
	elseif dim == 2 then
		elemDisc = FV1ConvectionDiffusionElemDisc2d()
	elseif dim == 3 then
		elemDisc = FV1ConvectionDiffusionElemDisc3d()
	else
	return nil
	end
	
	elemDisc:set_domain(domain)
	return elemDisc
end

function utilCreateFV1ConvDiffElemDisc(domain, diffusion, velocity, reaction, rhs, upwindAmount)
	local elemDisc = utilCreateFV1ConvDiffElemDiscNoInit(domain)

	elemDisc:set_upwind_amount(upwindAmount)
	elemDisc:set_diffusion_tensor(diffusion)
	elemDisc:set_velocity_field(velocity)
	elemDisc:set_reaction(reaction)
	elemDisc:set_rhs(rhs)
	return elemDisc
end

-- saves matrix for connection viewer
function SaveMatrixForConnectionViewer(gridFunc, linOp, fileName)
	local dim = gridFunc:get_dim()
	if dim == 1 then
		SaveMatrixForConnectionViewer1d(gridFunc, linOp, fileName)
	elseif dim == 2 then
		SaveMatrixForConnectionViewer2d(gridFunc, linOp, fileName)
	elseif dim == 3 then
		SaveMatrixForConnectionViewer3d(gridFunc, linOp, fileName)
	else
		return false
	end
	
	return true
end

-- saves vector for connection viewer
function SaveVectorForConnectionViewer(gridFunc, fileName)
	local dim = gridFunc:get_dim()
	if dim == 1 then
		SaveVectorForConnectionViewer1d(gridFunc, fileName)
	elseif dim == 2 then
		SaveVectorForConnectionViewer2d(gridFunc, fileName)
	elseif dim == 3 then
		SaveVectorForConnectionViewer3d(gridFunc, fileName)
	else
	return false
	end

	return true
end

-- writes grid function in vtk format
function WriteGridFunctionToVTK(gridFunc, fileName)
	local dim = gridFunc:get_dim()
	if dim == 1 then
		WriteGridFunctionToVTK1d(gridFunc, fileName)
	elseif dim == 2 then
		WriteGridFunctionToVTK2d(gridFunc, fileName)
	elseif dim == 3 then
		WriteGridFunctionToVTK3d(gridFunc, fileName)
	else
	return false
	end
	
	return true
end

-- applies the linear solver
function ApplyLinearSolver(linOp, u, b, linSolver)
	local dim = u:get_dim()
	if dim == 1 then
		ApplyLinearSolver1d(linOp, u, b, linSolver)
	elseif dim == 2 then
		ApplyLinearSolver2d(linOp, u, b, linSolver)
	elseif dim == 3 then
		ApplyLinearSolver3d(linOp, u, b, linSolver)
	else
	return false
	end
	
	return true
end

-- applies the linear solver
function utilCreateGeometricMultiGridPreconditioner(approxSpace)
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

-- creates a Lua User Matrix using a lua function and returns the Provider
function utilCreateLuaUserMatrix(funcName, dim)
	local mat = nil
	local provider = nil
	if dim == 1 then
		mat = LuaUserMatrix1d()
		provider = LuaUserMatrixProvider1d()
	elseif dim == 2 then
		mat = LuaUserMatrix2d()
		provider = LuaUserMatrixProvider2d()
	elseif dim == 3 then
		mat = LuaUserMatrix3d()
		provider = LuaUserMatrixProvider3d()
	end
		
	mat:set_lua_callback(funcName)
	provider:set_user_matrix(mat)
	return provider
end

-- creates a Const User Matrix using a lua function and returns the Provider
function utilCreateConstDiagUserMatrix(diagVal, dim)
	local mat = nil
	local provider = nil
	if dim == 1 then
		mat = ConstUserMatrix1d()
		provider = ConstUserMatrixProvider1d()
	elseif dim == 2 then
		mat = ConstUserMatrix2d()
		provider = ConstUserMatrixProvider2d()
	elseif dim == 3 then
		mat = ConstUserMatrix3d()
		provider = ConstUserMatrixProvider3d()
	end

	mat:set_diag_tensor(diagVal)
	provider:set_user_matrix(mat)
	return provider
end

-- creates a Lua User Vector using a lua function and returns the Provider
function utilCreateLuaUserVector(funcName, dim)
	local vec = nil
	local provider = nil
	if dim == 1 then
		vec = LuaUserVector1d()
		provider = LuaUserVectorProvider1d()
	elseif dim == 2 then
		vec = LuaUserVector2d()
		provider = LuaUserVectorProvider2d()
	elseif dim == 3 then
		vec = LuaUserVector3d()
		provider = LuaUserVectorProvider3d()
	end

	vec:set_lua_callback(funcName)
	provider:set_user_vector(vec)
	return provider
end

-- creates a Const User Vector using a lua function and returns the Provider
function utilCreateConstUserVector(val, dim)
	local vec = nil
	local provider = nil
	if dim == 1 then
		vec = ConstUserVector1d()
		provider = ConstUserVectorProvider1d()
	elseif dim == 2 then
		vec = ConstUserVector2d()
		provider = ConstUserVectorProvider2d()
	elseif dim == 3 then
		vec = ConstUserVector3d()
		provider = ConstUserVectorProvider3d()
	end
	
	vec:set_all_entries(val)
	provider:set_user_vector(vec)
	return provider
end

-- creates a Lua User Vector using a lua function and returns the Provider
function utilCreateLuaUserNumber(funcName, dim)
	local number = nil
	local provider = nil
	if dim == 1 then
		number = LuaUserNumber1d()
		provider = LuaUserNumberProvider1d()
	elseif dim == 2 then
		number = LuaUserNumber2d()
		provider = LuaUserNumberProvider2d()
	elseif dim == 3 then
		number = LuaUserNumber3d()
		provider = LuaUserNumberProvider3d()
	end
	
	number:set_lua_callback(funcName)
	provider:set_user_number(number)
	return provider
end

-- creates a Const User Number using a lua function and returns the Provider
function utilCreateConstUserNumber(val, dim)
	local number = nil
	local provider = nil
	if dim == 1 then
		number = ConstUserNumber1d()
		provider = ConstUserNumberProvider1d()
	elseif dim == 2 then
		number = ConstUserNumber2d()
		provider = ConstUserNumberProvider2d()
	elseif dim == 3 then
		number = ConstUserNumber3d()
		provider = ConstUserNumberProvider3d()
	end
	
	number:set(val)
	provider:set_user_number(number)
	return provider
end

