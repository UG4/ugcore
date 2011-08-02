----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Andreas Vogel / Martin Rupp
--
----------------------------------------------------------

-- make sure that ug_util is in the right path.
-- currently only the path in which you start your application is valid.
dofile("../scripts/ug_util.lua")

-- constants
dim = 2
gridName = "unit_square_tri_neumann.ugx"
numRefs = 4

-- name function (will be removed in the future, do not use them)
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
--local s = 2*math.pi
--return	s*s*(math.sin(s*x) + math.sin(s*y))
return -2*y
end

function ourNeumannBnd(x, y, t)
--local s = 2*math.pi
--return -s*math.cos(s*x)
return true, -2*x*y
end

function ourDirichletBnd(x, y, t)
--local s = 2*math.pi
--return true, math.sin(s*x) + math.sin(s*y)
return true, x*x*y
end

--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
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
sh:set_subset_name("DirichletBoundary", 1)
sh:set_subset_name("NeumannBoundary", 2)

-- create Refiner
print("Create Refiner")
refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom:get_grid())
for i=1,numRefs do
refiner:refine()
end

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()

-------------------------------------------
--  Setup User Functions
-------------------------------------------

-- Diffusion Tensor setup
diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor", dim)
--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
velocityField = util.CreateLuaUserVector("ourVelocityField", dim)
--velocityField = util.CreateConstUserVector(0.0, dim)

-- Reaction setup
reaction = util.CreateLuaUserNumber("ourReaction", dim)
--reaction = util.CreateConstUserNumber(0.0, dim)

-- rhs setup
rhs = util.CreateLuaUserNumber("ourRhs", dim)
--rhs = util.CreateConstUserNumber(0.0, dim)

-- neumann setup
neumann = util.CreateLuaBoundaryNumber("ourNeumannBnd", dim)
--neumann = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd", dim)
--dirichlet = util.CreateConstBoundaryNumber(0.0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiffElemDisc(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDisc = util.CreateNeumannBoundary(approxSpace, "Inner")
neumannDisc:add(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add(dirichlet, "c", "DirichletBoundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:set_approximation_space(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(neumannDisc)
domainDisc:add(dirichletBND)

-------------------------------------------
--  Algebra
-------------------------------------------

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

-- set initial value
u:set(1.0)

-- init Operator
linOp:init_op_and_rhs(b)

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)

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

-- base = LU()
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

transfer = P1ProlongationOperator2d()
transfer:set_approximation_space(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)

projection = P1ProjectionOperator2d()
projection:set_approximation_space(approxSpace)

gmg = util.CreateGeometricMultiGridPreconditioner(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_surface_level(numRefs)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_prolongation(transfer)
gmg:set_projection(projection)

amg = AMGPreconditioner()
amg:set_nu1(2)
amg:set_nu2(2)
amg:set_gamma(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
amg:set_max_levels(5)
-- amg:set_debug(u)

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(3)
convCheck:set_minimum_defect(1e-9)
convCheck:set_reduction(1e-9)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(amg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CGSolver()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- Apply Solver
ApplyLinearSolver(linOp, u, b, linSolver)

-- Output
WriteGridFunctionToVTK(u, "Solution")