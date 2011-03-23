----------------------------------------------------------
--
--   Lua - 	Script to perform the Laplace-Problem on
--	 		smooth multigrid hierarchies.
--
--   Authors: Andreas Vogel, Sebastian Reiter
--
----------------------------------------------------------

-- filenames are interpreted relative to APP_PATH/../scripts,
-- where APP_PATH is the path in which ug_shell is located.
-- If the script is not found, then the method will attempt to
-- load it relative to the current path.
ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 2
gridName = "volgrid_b.ugx"
numRefs = 4
useLoopScheme = true

-- name function
_C_ = 0

-- name subsets
_INNER_ = 0
_BND_ = 1

-- user functions
function ourDiffTensor(x, y, t)
return	1, 0, 
0, 1
end

function ourVelocityField(x, y, t)
return	0, 0
end

function ourReaction(x, y, t)
return	0
end

function ourRhs(x, y, t)
s = 2*math.pi
return	s*s*(math.sin(s*x) + math.sin(s*y))
end

	
-- create Instance of a Domain
print("Create Domain.")
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == true then
	-- create Refiner
	print("Create Hierarchy")
	if useLoopScheme == true then
		CreateSmoothHierarchy(dom:get_grid(), numRefs)
	else
		CreateSemiSmoothHierarchy(dom:get_grid(), numRefs)
	end
	
	print("Saving domain for debug reasons")
	-- write grid to file for test purpose
	SaveDomain(dom, "refined_grid.ugx")
	SaveGridHierarchy(dom:get_grid(), "refined_grid_hierarchy.ugx")

	print("creating approximation space...")
	-- create Approximation Space
	approxSpace = util.CreateApproximationSpace(dom)
	approxSpace:add_fct("c", "Lagrange", 1)
	approxSpace:init()

	-- name subsets
	sh = dom:get_subset_handler()
	if sh:num_subsets() ~= 2 then 
	   print("Domain must have 2 Subsets for this problem.")
	   exit()
	end
	sh:set_subset_name("Inner", 0)
	sh:set_subset_name("Boundary", 1)

	-- create Discretization
	domainDisc = DomainDiscretization()

	print("creating dirichlet boundary...")
	-- create dirichlet boundary and add it to discretization
	dirichletBND = util.CreateDirichletBnd2d(dom, SinusDirichletBoundaryFunction2d(), _C_)
	domainDisc:add_dirichlet_bnd(dirichletBND, _BND_)

	------------------------------------------
	--  Setting up user functions
	-------------------------------------------
	
	-- Diffusion Tensor setup
		diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor", dim)
	-- Velocity Field setup
		velocityField = util.CreateLuaUserVector("ourVelocityField", dim)
	-- Reaction setup
		reaction = util.CreateLuaUserNumber("ourReaction", dim)
	-- rhs setup
		rhs = util.CreateLuaUserNumber("ourRhs", dim)

	print("creating element discretization...")
	-- create Finite-Volume Element Discretization for Convection Diffusion Equation
	elemDisc = util.CreateFV1ConvDiffElemDisc(dom, diffusionMatrix, velocityField, reaction, rhs, 0)
	
	-- add Element Discretization to discretization
	domainDisc:add(elemDisc, pattern, "c", "Inner")

	print("creating linear operator...")	
	-- create operator from discretization
	linOp = AssembledLinearOperator()
	linOp:export_rhs(true)
	linOp:set_discretization(domainDisc)
	linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

	-- get grid function
	u = approxSpace:create_surface_function()
	b = approxSpace:create_surface_function()

	-- set initial value
	u:set(1.0)

	-- init Operator
	linOp:init()

	-- write matrix for test purpose
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")

	-- set dirichlet values in start iterate
	linOp:set_dirichlet_values(u)
	b:assign(linOp:get_rhs())

	-- create algebraic Preconditioner
	jac = JacobiPreconditioner()
	jac:set_damp(0.8)
	gs = GSPreconditioner()
	sgs = SGSPreconditioner()
	bgs = BGSPreconditioner()
	ilu = ILUPreconditioner()
	ilut = ILUTPreconditioner()

	-- create GMG
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(500)
	baseConvCheck:set_minimum_defect(1e-8)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(false)

	-- base = LapackLUSolver()
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)

	gmg = util.CreateGeometricMultiGridPreconditioner(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_surface_level(numRefs)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)

	-- create Convergence Check
	convCheck = StandardConvergenceCheck()
	convCheck:set_maximum_steps(1000)
	convCheck:set_minimum_defect(1e-12)
	convCheck:set_reduction(1e-12)

	-- create Linear Solver
	linSolver = LinearSolver()
	linSolver:set_preconditioner(gmg)
	linSolver:set_convergence_check(convCheck)

	-- create CG Solver
	cgSolver = CGSolver()
	cgSolver:set_preconditioner(ilu)
	cgSolver:set_convergence_check(convCheck)

	-- Apply Solver
	ApplyLinearSolver(linOp, u, b, linSolver)

	-- Output
	WriteGridFunctionToVTK(u, "Solution")

else	-- LoadDomain
	print("Load failed.")
end		--LoadDomain