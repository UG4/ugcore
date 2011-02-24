----------------------------------------------------------
--
--   Lua - 	Script to perform the Laplace-Problem on
--	 		adaptive multigrid hierarchies.
--
--   Authors: Andreas Vogel, Sebastian Reiter (only a small part!)
--
--	This script intends to demonstrate how adaptive-multigrid
--	hierarchies with hanging nodes can be used to solve
--	the laplace problem.
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraChooser());

------------------------
-- CONSTANTS
dim = 2
gridName = "open_circle.ugx"
--gridName = "unit_square_tri.ugx"

--refinement (default is 5)
numRefs    = GetParamNumber("-numRefs",    1)

-- all elements connected to vertices in this sphere will be refined
refCenterX = 0.0
refCenterY = 0.0
refCenterZ = 0
initialRadius = 0.1
-- in every refinement iteration the radius shrinks with this factor
radiusFalloff = 0.75

--------------------------------
-- User Data Functions (begin)
--------------------------------

-- 2D example, cf. Braess, Finite Elemente, p. 33
-- 
--	- \Delta u = 0                 in \Omega
--           u = sin(2/3 * \phi)   on circle bnd
--           u = 0                 on line bnd

function ourDiffTensor2d(x, y, t)
	return	1, 0, 
			0, 1
end
		
function ourDirichletBnd2d(x, y, t)
	-- line bnd
	local small = 0.0001
	if (-small < x and x < small) and y <= 0 then return true, 0.0 end
	if (-small < y and y < small) and x >= 0 then return true, 0.0 end
	
	-- circle bnd
	local phi;
	local r = math.sqrt(x*x + y*y)
	if y >= 0 then phi = math.acos(x/r) end
	if y < 0 then phi =2*math.pi - math.acos(x/r) end
	
	return true, math.sin((2*phi)/3)
end

-- 3D

function ourDiffTensor3d(x, y, z, t)
	return	1, 0, 0,
			0, 1, 0,
			0, 0, 1
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
dom = utilCreateDomain(dim)

-- load domain
print("Load Domain from File.")
if utilLoadDomain(dom, gridName) == false then print("Loading Domain failed."); exit() end

-- Distribute the domain to all involved processes
if DistributeDomain(dom) == false then print("Error while Distributing Grid."); exit() end

-- get subset handler and make sure that the right number of subsets is contained.
sh = dom:get_subset_handler()
if sh:num_subsets() ~= 2 then print("Domain must have 2 Subsets for this problem."); exit() end

-- create refiner
print("Create Hierarchy")
refiner = HangingNodeDomainRefiner(dom);

-- refine
local radius = initialRadius
for i = 1, numRefs do
	MarkForRefinement_VerticesInSphere(refiner, refCenterX, refCenterY, refCenterZ, radius)
	
	refiner:refine()
	radius = radius * radiusFalloff
end

-- write grid to file for test purpose
print("Saving domain grid and hierarchy.")
SaveDomain(dom, "refined_grid.ugx")
SaveGridHierarchy(dom:get_grid(), "refined_grid_hierarchy.ugx")

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1Function(pattern, "c", dim)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Diffusion Tensor setup
diffusionMatrix = utilCreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
--diffusionMatrix = utilCreateConstDiagUserMatrix(1.0, dim)

-- dirichlet setup
-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDirichletBnd2d")
dirichlet = utilCreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

-- create a Hanging Node Finite Volume Assembling
-- Note, that only the diffusion tensor is set. All other possible user data
-- (such as Convection Velocity, Reaction Term, Source Term) are not set and
-- the discretization uses the default (i.e. zero) value
elemDisc = utilCreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)

-----------------------------------------------------------------
--  Setup Constraints
-----------------------------------------------------------------

constraints = OneSideP1Constraints()

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = utilCreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(dirichlet, "c", "Boundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
domainDisc:add_post_process(dirichletBND)
domainDisc:add_post_process(constraints)

-------------------------------------------
--  Setup Operator
-------------------------------------------

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:export_rhs(true)
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- debug writer
dbgWriter = utilCreateGridFunctionDebugWriter(dim)
dbgWriter:set_vtk_output(false)

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
	-- base = LUSolver()
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
	
	-- Transfer and Projection
	transfer = utilCreateP1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	projection = utilCreateP1Projection(approxSpace)
	
	-- Gemoetric Multi Grid
	gmg = utilCreateGeometricMultiGrid(approxSpace)
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
	gmg:set_debug(dbgWriter)

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- exact Solver
exSolver = LU()

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

-------------------------------------------
--  Apply Solver
-------------------------------------------

-- choose some solver
solver = linSolver

-- create grid function
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()
dbgWriter:set_reference_grid_function(u)

-- set initial value
u:set(0.0)

-- 1. init operator
print("Init operator (i.e. assemble matrix).")
if linOp:init() ~= true then print("Cannot init operator"); exit() end

-- 2. set dirichlet values in solution
linOp:set_dirichlet_values(u)

-- 3. set right hand side (assembled together with Operator for performance reasons)
b:assign(linOp:get_rhs())

-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.mat")

-- 4. init solver for linear Operator
print("Init solver for operator.")
if solver:init(linOp) ~= true then print("Cannot init solver"); exit() end

-- 5. apply solver
print("Apply solver.")
if solver:apply_return_defect(u,b) ~= true then print("Cannot apply solver"); exit() end

-------------------------------------------
--  Output
-------------------------------------------
WriteGridFunctionToVTK(u, "Solution")