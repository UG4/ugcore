----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Martin Rupp / Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = util.GetParam("-dim", 2)


if dim == 2 then
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
end
if dim == 3 then
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs", 2)

--------------------------------
-- User Data Functions (begin)
--------------------------------
function ourDiffTensor2d(x, y, t)
return	1, 0, 
0, 1
end

function ourVelocityField2d(x, y, t)
return	0, 0
end

function ourReaction2d(x, y, t)
return	0
end

function ourRhs2d(x, y, t)
local s = 2*math.pi
return	s*s*(math.sin(s*x) + math.sin(s*y))
--return -2*y
--return 0;
end

function ourNeumannBnd2d(x, y, t)
--local s = 2*math.pi
--return -s*math.cos(s*x)
return true, -2*x*y
end

function ourDirichletBnd2d(x, y, t)
local s = 2*math.pi
return true, math.sin(s*x) + math.sin(s*y)
--return true, x*x*y
--return true, 2.5
end

function ourDiffTensor3d(x, y, z, t)
return	1, 0, 0,
0, 1, 0,
0, 0, 1
end

function ourVelocityField3d(x, y, z, t)
return	0, 0, 0
end

function ourReaction3d(x, y, z, t)
return	0
end

function ourRhs3d(x, y, z, t)
--local s = 2*math.pi
--return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
return 0;
end

function ourNeumannBnd3d(x, y, t)
--local s = 2*math.pi
--return -s*math.cos(s*x)
return true, -2*x*y*z
end

function ourDirichletBnd3d(x, y, z, t)
--local s = 2*math.pi
--return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
return true, x
end
--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
print("Loading Domain failed.")
exit()
end

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

if DistributeDomain(dom) == false then
print("Error while Distributing Grid.")
exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
util.GlobalRefineParallelDomain(dom)
end

-- get subset handler
sh = dom:get_subset_handler()
if sh:num_subsets() ~= 2 then 
print("Domain must have 2 Subsets for this problem.")
exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

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
print ("Setting up Assembling")

-- Diffusion Tensor setup
if dim == 2 then
diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor2d", dim)
elseif dim == 3 then
diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor3d", dim)
end
--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
if dim == 2 then
velocityField = util.CreateLuaUserVector("ourVelocityField2d", dim)
elseif dim == 3 then
velocityField = util.CreateLuaUserVector("ourVelocityField3d", dim)
end 
--velocityField = util.CreateConstUserVector(0.0, dim)

-- Reaction setup
if dim == 2 then
reaction = util.CreateLuaUserNumber("ourReaction2d", dim)
elseif dim == 3 then
reaction = util.CreateLuaUserNumber("ourReaction3d", dim)
end
--reaction = util.CreateConstUserNumber(0.0, dim)

-- rhs setup
if dim == 2 then
rhs = util.CreateLuaUserNumber("ourRhs2d", dim)
elseif dim == 3 then
rhs = util.CreateLuaUserNumber("ourRhs3d", dim)
end
--rhs = util.CreateConstUserNumber(0.0, dim)

-- neumann setup
if dim == 2 then
neumann = util.CreateLuaBoundaryNumber("ourNeumannBnd2d", dim)
elseif dim == 3 then
neumann = util.CreateLuaBoundaryNumber("ourNeumannBnd3d", dim)
end
--neumann = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
if dim == 2 then
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd2d", dim)
elseif dim == 3 then
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd3d", dim)
end
--dirichlet = util.CreateConstBoundaryNumber(3.2, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = util.CreateNeumannBoundary(approxSpace, "Inner")
--neumannDisc:add_boundary_value(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
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
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

print ("Reset initial value")
-- set initial value
u:set(0.0)

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
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()

-- create GMG ---
-----------------

-- Base Solver
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-8)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose_level(false)
-- base = LapackLUSolver()
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

-- Transfer and Projection
transfer = util.CreateP1Prolongation(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)
projection = util.CreateP1Projection(approxSpace)

-- Gemoetric Multi Grid
gmg = util.CreateGeometricMultiGrid(approxSpace)
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

-- create AMG ---
-----------------

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(jac)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(jac)
bicgstabSolver:set_convergence_check(convCheck)

-- Apply Solver
tBefore = os.clock()

ApplyLinearSolver(linOp, u, b, linSolver)

tAfter = os.clock()

if(GetProfilerAvailable()) then
	print("Profiler is available!\n")
else
	print("Profiler not available.\n")
	-- however os.clock will work
end


-- Get profiler node information
assemblePN = GetProfileNode("assembleLinearMatrix")
-- functions are assemblePN:get_avg_total_time_ms(), assemblePN:get_avg_self_time_ms(), assemblePN:get_avg_entry_count() und assemblePN:is_valid()

 
if(assemblePN:is_valid()) then
	print("assembleLinearMatrix found!")
else
	print("profile node for assembleLinearMatrix not found!")
end

s = string.format("%d\t%d\t%d\t%.2f\t%.2f\n", dim, numPreRefs, numRefs, tAfter-tBefore, assemblePN:get_avg_total_time_ms())
output = io.open("output.txt", "a")
output:write(s)
print("dim\tnumPreRefs\tnumRefs\tt\tassembleLinearMatix");
print(s)



-- Output
