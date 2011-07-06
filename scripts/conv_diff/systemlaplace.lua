----------------------------------------------------------
--
--   Lua - Script to perform the several Laplace-Problems
--
--	 Description:	For a testing purpose we discretize
--					#<blocks> completely decoupled laplace equations
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")


nSystems = util.GetParamNumber("-block", 1)

-- choose algebra
algebra = CPUAlgebraSelector()
algebra:set_fixed_blocksize(nSystems)
InitAlgebra(algebra)
-- InitAlgebra also loads all discretization functions and classes

-- constants
if 		util.HasParamOption("-1d") == true then dim = 1
elseif 	util.HasParamOption("-3d") == true then dim = 3
else dim = 2 end

if 		dim == 2 then gridName = "unit_square/unit_square_quads_8x8.ugx"
elseif 	dim == 3 then gridName = "unit_square/unit_cube_hex.ugx"
					--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

numPreRefs = util.GetParamNumber("-numPreRefs", 2)
numRefs = util.GetParamNumber("-numRefs", 3)

print("\nSYSTEMLAPLACE, nSystems = "..nSystems..", Dim = "..dim..", numRefs = "..numRefs.."\n")

--------------------------------
-- User Data Functions (begin)
--------------------------------
function ourDiffTensor2d(x, y, t)
	return	1, 0, 
			0, 1
end

function ourRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

function ourDirichletBnd2d(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
end

function ourDiffTensor3d(x, y, z, t)
	return	1, 0, 0,
			0, 1, 0,
			0, 0, 1
end


function ourRhs3d(x, y, z, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
end

function ourNeumannBnd3d(x, y, t)
	local s = 2*math.pi
	return -s*math.cos(s*x)
end

function ourDirichletBnd3d(x, y, z, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
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

-- Now we loop all subsets an search for it in the SubsetHandler of the domain
neededSubsets = {"Inner", "Boundary"}
if util.CheckSubsets(dom, neededSubsets) then print("Wrong subsets detected.") end

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
--diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- rhs setup
rhs = util.CreateLuaUserNumber("ourRhs"..dim.."d", dim)
--rhs = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
--dirichlet = util.CreateConstBoundaryNumber(0.0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------
elemDisc = {}
upwind = {}
for i=1, nSystems do
	upwind[i] = FullUpwind2d()
	elemDisc[i] = util.CreateFV1ConvDiff(approxSpace, "c"..i, "Inner")
	elemDisc[i]:set_upwind(upwind[i])
	elemDisc[i]:set_diffusion_tensor(diffusionMatrix)
	elemDisc[i]:set_source(rhs)
end

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
for i=1, nSystems do
	dirichletBND:add_boundary_value(dirichlet, "c"..i, "Boundary")
end

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
for i=1, nSystems do
	domainDisc:add_elem_disc(elemDisc[i])
end
domainDisc:add_post_process(dirichletBND)

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

print ("Reset initial value")
-- set initial value
u:set(1.0)

-- init Operator
print ("Assemble Operator ... ")
linOp:init_op_and_rhs(b)
print ("done")

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)

-- write matrix for test purpose
-- SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
-- SaveVectorForConnectionViewer(b, "Rhs.mat")

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
amg = RSAMGPreconditioner()
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

-- Apply Solver
ApplyLinearSolver(linOp, u, b, linSolver)

-- Output
-- WriteGridFunctionToVTK(u, "Solution")