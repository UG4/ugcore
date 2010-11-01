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
gridName = "unit_square_tri.ugx"
numPreRefs = 0
numRefs = 4

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
	--return -2*y
end

function ourNeumannBnd(x, y, t)
	--local s = 2*math.pi
	--return -s*math.cos(s*x)
	return true, -2*x*y
end

function ourDirichletBnd(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
	--return true, x*x*y
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
--if sh:num_subsets() ~= 2 then 
--	print("Domain must have 2 Subsets for this problem.")
--	exit()
--end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- create Refiner
print("Create Refiner")
if numPreRefs > numRefs then
	print("numPreRefs must be smaller/equal than numRefs");
	exit();
end

refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom:get_grid())
for i=1,numPreRefs do
	refiner:refine()
end

if DistributeDomain2d(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	GlobalRefineParallelDomain2d(dom)
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
print ("Setting up Assembling")

-- Diffusion Tensor setup
	diffusionMatrix = utilCreateLuaUserMatrix("ourDiffTensor", dim)
	--diffusionMatrix = utilCreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
	velocityField = utilCreateLuaUserVector("ourVelocityField", dim)
	--velocityField = utilCreateConstUserVector(0.0, dim)

-- Reaction setup
	reaction = utilCreateLuaUserNumber("ourReaction", dim)
	--reaction = utilCreateConstUserNumber(0.0, dim)

-- rhs setup
	rhs = utilCreateLuaUserNumber("ourRhs", dim)
	--rhs = utilCreateConstUserNumber(0.0, dim)

-- neumann setup
	neumann = utilCreateLuaBoundaryNumber("ourNeumannBnd", dim)
	--neumann = utilCreateConstUserNumber(0.0, dim)

-- dirichlet setup
	dirichlet = utilCreateLuaBoundaryNumber("ourDirichletBnd", dim)
	--dirichlet = utilCreateConstBoundaryNumber(0.0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = utilCreateFV1ConvDiffElemDisc(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_rhs(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = utilCreateNeumannBoundary(approxSpace, "Inner")
--neumannDisc:add_boundary_value(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = utilCreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(dirichlet, "c", "DirichletBoundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
--domainDisc:add_elem_disc(neumannDisc)
domainDisc:add_post_process(dirichletBND)

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

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
print ("Assemble Operator ... ")
linOp:init()
print ("done")

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
-- ilut = ILUTPreconditioner()

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

transfer = P1ProlongationOperator2d()
transfer:set_approximation_space(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)

projection = P1ProjectionOperator2d()
projection:set_approximation_space(approxSpace)

gmg = utilCreateGeometricMultiGridPreconditioner(approxSpace)
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

if false then
amg = AMGPreconditioner()
amg:set_nu1(2)
amg:set_nu2(2)
amg:set_gamma(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
--amg:set_debug(u)
end

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(jac)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CGSolver()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStabSolver()
bicgstabSolver:set_preconditioner(jac)
bicgstabSolver:set_convergence_check(convCheck)

-- Apply Solver
ApplyLinearSolver(linOp, u, b, linSolver)

-- Output
WriteGridFunctionToVTK(u, "Solution")