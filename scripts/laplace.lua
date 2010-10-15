----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

-- make sure that ug_util is in the right path.
-- currently only the path in which you start your application is valid.
dofile("../scripts/ug_util.lua")

-- constants
dim = 2
gridName = "../scripts/unit_square_tri_neumann.ugx"
numRefs = 4

-- name function (will be removed in the future, do not use them)
_C_ = 0
_INNER_ = 0
_BND_ = 1
_NEUMANN_BND_ = 2

--------------------------------
-- User Data Functions (begin)
--------------------------------

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
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

function ourNeumannBnd(x, y, t)
	local s = 2*math.pi
	return -s*math.cos(s*x)
end

--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom = utilCreateDomain(dim)

-- load domain
print("Load Domain from File.")
if utilLoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end

-- get subset handler
sh = dom:get_subset_handler()
if sh:num_subsets() ~= 3 then 
	print("Domain must have 3 Subsets for this problem.")
	exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("Boundary", 1)

-- create Refiner
print("Create Refiner")
refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom:get_grid())
for i=1,numRefs do
	refiner:refine()
end

-- write grid to file for test purpose
utilSaveDomain(dom, "refined_grid.ugx")

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1Function(pattern, "c", 2)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-------------------------------------------
--  Setup User Functions
-------------------------------------------

-- Diffusion Tensor setup
	--diffusionMatrix = utilCreateLuaUserMatrix("ourDiffTensor", dim)
	diffusionMatrix = utilCreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
	--velocityField = utilCreateLuaUserVector("ourVelocityField", dim)
	velocityField = utilCreateConstUserVector(0.0, dim)

-- Reaction setup
	--reaction = utilCreateLuaUserNumber("ourReaction", dim)
	reaction = utilCreateConstUserNumber(0.0, dim)

-- rhs setup
	rhs = utilCreateLuaUserNumber("ourRhs", dim)
	--rhs = utilCreateConstUserNumber(0.0, dim)

-- neumann setup
	neumann = utilCreateLuaUserNumber("ourNeumannBnd", dim)
	--neumann = utilCreateConstUserNumber(0.0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = utilCreateFV1ConvDiffElemDiscNoInit(dom)
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_rhs(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDisc = utilCreateNeumannBoundary(dom, neumann, _NEUMANN_BND_)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

-- create dirichlet boundary
bndFct = SinusDirichletBoundaryFunction2d()
dirichletBND = DirichletBND2d()
dirichletBND:set_domain(dom)
dirichletBND:set_dirichlet_function(bndFct)
dirichletBND:set_function(_C_)

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add(elemDisc, pattern, "c", "Inner")
domainDisc:add(neumannDisc, pattern, "c", "Inner")
domainDisc:add_dirichlet_bnd(dirichletBND, _BND_)

-------------------------------------------
--  Algebra
-------------------------------------------

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:export_rhs(true)
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function("u", true)
b = approxSpace:create_surface_function("b", true)

-- set initial value
u:set(1.0)

-- init Operator
linOp:init()

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)
b:assign(linOp:get_rhs())

-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.mat")

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

gmg = utilCreateGeometricMultiGridPreconditioner(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_surface_level(numRefs)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)

amg = AMGPreconditioner()
amg:set_nu1(2)
amg:set_nu2(2)
amg:set_gamma(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_debug(u)

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