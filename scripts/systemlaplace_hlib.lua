----------------------------------------------------------
--
--   Lua - Script to perform the several Laplace-Problems
--
--	 Description:	For a testing purpose we discretize
--					#<blocks> completely decoupled laplace equations
--
--   Author: Andreas Vogel - testwise solving with HLib arithmetic added (14032011ih)
--
----------------------------------------------------------
-- ugshell -ex ../scripts/systemlaplace_hlib.lua -block 3 -numPreRefs 0 -numRefs 0

ug_load_script("ug_util.lua")


nSystems = util.GetParamNumber("-block", 3)

-- choose algebra
algebra = CPUAlgebraChooser()
algebra:set_fixed_blocksize(nSystems)
InitAlgebra(algebra)
-- InitAlgebra also loads all discretization functions and classes


-- constants
if util.HasParamOption("-3d") == true then
	dim = 3
else
	dim = 2
end

if dim == 2 then
gridName = "unit_square_tri.ugx" -- "unit_square_quads_8x8.ugx"
end
if dim == 3 then
gridName = "unit_cube_hex.ugx"
--gridName = "unit_cube_tets_regular.ugx"
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs = util.GetParamNumber("-numRefs", 0)

print("\nSYSTEMLAPLACE, nSystems = "..nSystems..", Dim = "..dim..", numRefs = "..numRefs.."\n")

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

-- Performing pre-refines
for i=1,numPreRefs do
	refiner:refine()
end

-- Distribute the domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
util.GlobalRefineParallelDomain(dom)
end

-- write grid to file for test purpose
-- SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
for i=1, nSystems do
	approxSpace:add_fct("c"..i, "Lagrange", 1)
end
approxSpace:init()
approxSpace:print_statistic()

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
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
if dim == 2 then
	velocityField = util.CreateLuaUserVector("ourVelocityField2d", dim)
elseif dim == 3 then
	velocityField = util.CreateLuaUserVector("ourVelocityField3d", dim)
end 
velocityField = util.CreateConstUserVector(0.0, dim)

-- Reaction setup
if dim == 2 then
	reaction = util.CreateLuaUserNumber("ourReaction2d", dim)
elseif dim == 3 then
	reaction = util.CreateLuaUserNumber("ourReaction3d", dim)
end
reaction = util.CreateConstUserNumber(0.0, dim)

-- rhs setup
if dim == 2 then
	rhs = util.CreateLuaUserNumber("ourRhs2d", dim)
elseif dim == 3 then
	rhs = util.CreateLuaUserNumber("ourRhs3d", dim)
end
rhs = util.CreateConstUserNumber(0.0, dim)

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
--dirichlet = util.CreateConstBoundaryNumber(0.0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------
elemDisc = {}
for i=1, nSystems do
elemDisc[i] = util.CreateFV1ConvDiff(approxSpace, "c"..i, "Inner")
elemDisc[i]:set_upwind_amount(0.0)
elemDisc[i]:set_diffusion_tensor(diffusionMatrix)
elemDisc[i]:set_velocity_field(velocityField)
elemDisc[i]:set_reaction(reaction)
elemDisc[i]:set_source(rhs)
end

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = util.CreateNeumannBoundary(approxSpace, "Inner")
--neumannDisc:add_boundary_value(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
for i=1, nSystems do
dirichletBND:add_boundary_value(dirichlet, "c"..i, "DirichletBoundary")
end

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
for i=1, nSystems do
domainDisc:add_elem_disc(elemDisc[i])
end
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
u:set(0.0) -- 'systemlaplace.lua: '1.0'

-- init Operator
print ("Assemble Operator ... ")
linOp:init()
print ("done")

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)
b:assign(linOp:get_rhs())

-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, linOp, "SysLaplace-Stiffness.mat")
-- SaveVectorForConnectionViewer(b, "SysLaplace-Rhs.mat")

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
	base = LU()
	--base = LinearSolver()
	--base:set_convergence_check(baseConvCheck)
	--base:set_preconditioner(jac)
	
	-- Transfer and Projection
	transfer = util.CreateP1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	projection = util.CreateP1Projection(approxSpace)
	
	-- Geometric Multi Grid
	gmg = util.CreateGeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
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
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(jac)
bicgstabSolver:set_convergence_check(convCheck)

-- create HLIBSolver - Test (14032011ih)
linSolver2 = HLIBSolver()
linSolver2:set_hlib_accuracy_H(1.e-4)  -- default: 1.e-4
linSolver2:set_hlib_accuracy_LU(1.e-1) -- default: 1.e-4
-- define construction of cluster tree
--   first  arg: "clustering type" \in [algebraic | geometric (not yet implemented)]; algebraic is default
--   second arg: "clustering mode" \in [nested dissection | empty/everything else]; nested dissection is default 
linSolver2:set_clustering_method("algebraic", "nested dissection")
                                                                   
linSolver2:set_hlib_verbosity(4) -- '>= 2': create HLIB related postscripts; '>=3' also create plots of matrix entries
--linSolver2:set_ps_basename("SysLP")

-- Apply Solver
--ApplyLinearSolver(linOp, u, b, linSolver)
ApplyLinearSolver(linOp, u, b, linSolver2)

-- Output
WriteGridFunctionToVTK(u, "SystemLaplace-Solution")

--linSolver2:check_crs_matrix()