-- Test script for the downwind ordering

PluginRequired("LevelSet")

--------------------------------------------------------------------------------
--  Initialization
--------------------------------------------------------------------------------

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")

gridName = "unit_square_unstructured_tris_coarse_left_dirichlet.ugx"

numPreRefs = util.GetParamNumber("-numPreRefs", 0, "Number of pre-refinements")
numRefs = util.GetParamNumber("-numRefs", 4, "Number of refinements")

-- initialize ug with the world dimension and the algebra type
InitUG(2, AlgebraType("CPU", 1))

-- Load a domain without initial refinements.
dom = util.CreateDomain(gridName, 0)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")

util.refinement.CreateRegularHierarchy(dom, numRefs, true)

-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("u", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


--------------------------------------------------------------------------------
--  Velocity
--------------------------------------------------------------------------------


function velocity(x, y)
	local x_dir = 1
	local y_dir = 0
	return x_dir, y_dir
end

--------------------------------------------------------------------------------
--  Discretization
--------------------------------------------------------------------------------

upwind = FullUpwind()

elemDisc = FV1_Convection(upwind, "u", "Inner")
elemDisc:set_velocity(velocity)

DirichletBnd = DirichletBoundary()
DirichletBnd:add(1.0, "u", "Dirichlet")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(DirichletBnd)

--------------------------------------------------------------------------------
--  Linear solver (choose!) - using 'util/solver_util.lua'
--------------------------------------------------------------------------------


ilu = ILU()
ilu:set_damp(1)

--algo = BoostCuthillMcKeeOrdering() --old version not using UndirectedMatrix
algo = BoostCuthillMcKeeNewOrdering() --using UndirectedMatrix

ilu:set_ordering_algorithm(algo)


-- iterative
ilu_solverDesc = {
	type = "linear",
	precond = ilu,
	convCheck = {
		type		= "standard",
		iterations	= 100,
		absolute	= 1e-12,
		reduction	= 1e-10,
		verbose	= true
	}
}
solver = util.solver.CreateSolver(ilu_solverDesc)


--------------------------------------------------------------------------------
--  Assemble the discretized problem
--------------------------------------------------------------------------------


print("Assembling...")

A = AssembledLinearOperator(domainDisc)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)


print("--------------------------------------------------------------------------")

--------------------------------------------------------------------------------
--  Solve the discretized linear problem
--------------------------------------------------------------------------------

print("Solving...")

solver:init(A, u)
solver:apply(u, b)

--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------

solFileName = "ordertest"

print("writing solution to '" .. solFileName .. ".vtu'...")


vtk_out = VTKOutput()
vtk_out:clear_selection()
vtk_out:select_nodal(GridFunctionNumberData(u, "u"), "u")
vtk_out:select_nodal(LuaUserVector("velocity"), "vel")
vtk_out:set_binary(false)
vtk_out:set_write_proc_ranks(true)

vtk_out:print(solFileName, u)

print("done")


