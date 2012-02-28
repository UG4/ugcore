--------------------------------------------------------------------------------
--	tut04_2_disc_laplace.lua
--
--	This file is part of tut04_modular_programming.lua.
--
--	For a more detailed documentation on how to assemble the laplace operator
--	please see tutorial 3: tut03_laplace_with_exact_solver.lua
--
--	This file contains the method AssembleLaplace.
--	The method will set up the required approximation space and discretization
--	objects and will assemble the linear operator.
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("../ug_util.lua")


--------------------------------------------------------------------------------
-- AssembleLaplace assembles the laplace operator. It returns a linear operator
-- containing the created matrix and the underlying approximation space
--
-- You have to specify a domain and the names of the inner and boundary subsets.
-- Optionally you may specify callbacks for the diffusion tensor, the rhs and
-- the dirichlet boundary conditions. Note that you have to pass the names of
-- the callback methods as strings.
--
-- The assembled system operates on the surface-view of the specified domain
-- and is not usable in geometric multigrid methods (this could be enhanced
-- by an additional)
--
-- Params:
--   * dom: A domain object.
--	 * innerSubsets: A string that contains the subset names for the inner
--			subsets. Multiple names can be separated by comma.
--   * boundarySubsets: A string that contains the subset names for the boundary
--			subsets. Multiple names can be separated by comma.
--   * cbDiffTensorName (optional): string that contains the name of the callback
--			which defines the diffusion tensor. Can be nil.
--			The callback has to return 4 number values in 2d and 9 in 3d,
--			defining the coefficients of the diffusion matrix. Parameters
--			are the coordinates and the time. 
--			See LaplaceDefaults_DiffTensor2d / 3d below for an example.
--   * cbRhs (optional): string that contains the name of the callback
--			which defines the right hand side. Can be nil.
--			The callback has to return 1 number value both in 2d and in 3d,
--			defining the value of the rhs at the given point. Parameters
--			are the coordinates and the time.
--			See LaplaceDefaults_RHS2d / 3d below for an example.
--   * cbDirichletBndName (optional): string that contains the name of the callback
--			which defines the dirichlet boundary values. Can be nil.
--			The callback has to return a boolean and 1 number value both in 2d and
--			in 3d, defining the dirichlet boundary values at the given point.
--			Parameters are the coordinates and the time.
--			The returned boolean (true or false) defines whether the specified
--			point is indeed regarded as a dirichlet point.
--			See LaplaceDefaults_DirichletBnd2d / 3d below for an example.
--
-- Returns:
--	 * AssembledLinearOperator: The assembled Linear operator A.
--   * ApproximationSpace: The underlying approximation space.
--
function AssembleLaplace(dom, innerSubsets, boundarySubsets, b,
						cbDiffTensorName, cbRhsName, cbDirichletBndName)
--	init some local variables
	local dim = dom:get_dim()

--	create the approximation space
	local approxSpace = ApproximationSpace(dom) -- creates new object
	approxSpace:add_fct("c", "Lagrange", 1)          -- adds one function
	approxSpace:print_statistic()                    -- write some information
	
--	initialize the callbacks. If callbacks were specified, use them. If not,
--	use the default callbacks.
	local tmpDiffCBName = cbDiffTensorName
	if tmpDiffCBName == nil then
		tmpDiffCBName = "LaplaceDefaults_DiffTensor" .. dim .. "d"
	end
	
	local tmpRhsCBName = cbRhsName
	if tmpRhsCBName == nil then
		tmpRhsCBName = "LaplaceDefaults_RHS" .. dim .. "d"
	end
	
	local tmpBndCBName = cbDirichletBndName
	if tmpBndCBName == nil then
		tmpBndCBName = "LaplaceDefaults_DirichletBnd" .. dim .. "d"
	end

	local cbDiff = LuaUserMatrix(tmpDiffCBName)
	local cbRhs = LuaUserNumber(tmpRhsCBName)
	local cbBnd = LuaBoundaryNumber(tmpBndCBName)
	
	
--	set up the discretization
	local elemDisc = ConvectionDiffusion("c", innerSubsets)
	elemDisc:set_disc_scheme("fv1")
	if dim == 1 then
		upwind = NoUpwind1d() -- create an upwind procedure ("No Upwind")
	elseif dim == 2 then
		upwind = NoUpwind2d() -- create an upwind procedure ("No Upwind")
	elseif dim == 3 then
		upwind = NoUpwind3d() -- create an upwind procedure ("No Upwind")
	end
	elemDisc:set_upwind(upwind)  -- set the upwind procedure
	elemDisc:set_diffusion_tensor(cbDiff)	-- set the diffusion matrix
	elemDisc:set_source(cbRhs)					-- set the right hand side
	
	local dirichletBnd = DirichletBoundary()
	dirichletBnd:add(cbBnd, "c", boundarySubsets)
	
	local domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc)
	domainDisc:add(dirichletBnd)
	
	
--	create the linear operator
	local linOp = AssembledLinearOperator()
	linOp:set_discretization(domainDisc)

	linOp:init_op_and_rhs(b)

	return linOp, approxSpace, domainDisc, dirichletBnd
end



--------------------------------------------------------------------------------
--	The default callbacks
--
--	The default callbacks are chosen so that the exact solution in 2d is
--	f(x, y) = x*x + y*y
--	and in 3d
--	f(x, y, z) = x*x + y*y + z*z
--------------------------------------------------------------------------------

function LaplaceDefaults_DiffTensor2d(x, y, t)
	return	1, 0, 
			0, 1
end

function LaplaceDefaults_DiffTensor3d(x, y, z, t)
	return	1, 0, 0,
			0, 1, 0,
			0, 0, 1
end

-- Dirichlet boundary callbacks
function LaplaceDefaults_DirichletBnd2d(x, y, t)
	return true, x*x + y*y
end

function LaplaceDefaults_DirichletBnd3d(x, y, z, t)
	return true, x*x + y*y + z*z
end

-- RHS callbacks
function LaplaceDefaults_RHS2d(x, y, t)
	return	-4
end

function LaplaceDefaults_RHS3d(x, y, z, t)
	return	-6
end
