----------------------------------------------------------
--
--   Lua - 	Script to perform the Laplace-Problem on
--	 		adaptive multigrid hierarchies.
--
--   Authors: Andreas Vogel, Sebastian Reiter
--
--	This script intends to demonstrate how adaptive-multigrid
--	hierarchies with hanging nodes can be used to solve
--	the laplace problem.
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

------------------------
-- CONSTANTS
dim = util.GetParamNumber("-dim",    2)
if dim == 2 then
	gridName = "grids/open_circle.ugx"
--	gridName = "unit_square/unit_square_quads_8x8.ugx"
--	gridName = "unit_square/unit_square_tri_8x8.ugx"
--	gridName = "unit_square/unit_square_tri_2x2.ugx"
elseif dim == 3 then
	gridName = "grids/open_cube_hex.ugx"
else
	print("Dimension " .. dim .. " not supported. Aborting.")
	exit()
end

--refinement (default is 5)
numRefs    = util.GetParamNumber("-numRefs",    1)

-- all elements connected to vertices in this sphere will be refined
refCenterX = 0.0
refCenterY = 0.0
refCenterZ = 0
initialRadius = 0.5
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
	
--	return true, 10;
end

-- 3D

function ourDiffTensor3d(x, y, z, t)
	return	1, 0, 0,
			0, 1, 0,
			0, 0, 1
end

function ourDirichletBnd3d(x, y, z, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
--	return true, x
end
--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then print("Loading Domain failed."); exit() end

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
local refFrom = 0.3
local refStep = 0.2
local refCurr = refFrom

local TriagCenterX = 1/6
local TriagCenterY = 1/6
local TriagScale = 0.5;
local TriagRadius = 0.4;
for i = 1, numRefs do
	MarkForRefinement_VerticesInSphere(refiner, refCenterX, refCenterY, refCenterZ, radius)
--	MarkForRefinement_VerticesInSquare(refiner, refCurr, 1, -1, 1, 0, 0)
--	MarkForRefinement_VerticesInSquare(refiner, 0.2*i, 1, -1, 1, 0, 0)
--	MarkForRefinement_FacesInSphere(refiner, TriagCenterX, TriagCenterY, 0.0, TriagRadius)
	
	refiner:refine()
	
	refCurr = refCurr + refStep/(i+2)
	radius = radius * radiusFalloff
	TriagRadius = TriagRadius * TriagScale
end

-- write grid to file for test purpose
print("Saving domain grid and hierarchy.")
SaveDomain(dom, "refined_grid.ugx")
SaveGridHierarchy(dom:get_grid(), "refined_grid_hierarchy.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()

SaveGridHierarchyTransformed(dom:get_grid(), approxSpace:get_surface_view(), "surface_view_extrude.ugx", 0.25)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Diffusion Tensor setup
diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- dirichlet setup
-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDirichletBnd2d")
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

-- create a Hanging Node Finite Volume Assembling
-- Note, that only the diffusion tensor is set. All other possible user data
-- (such as Convection Velocity, Reaction Term, Source Term) are not set and
-- the discretization uses the default (i.e. zero) value
elemDisc = util.CreateFE1ConvDiff(approxSpace, "c", "Inner")
--elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)

-----------------------------------------------------------------
--  Setup Constraints
-----------------------------------------------------------------

constraints = SymP1Constraints()

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
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
dbgWriter = util.CreateGridFunctionDebugWriter(dim)
dbgWriter:set_vtk_output(false)

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(1.0)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()

-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(100)
	baseConvCheck:set_minimum_defect(1e-10)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(true)
	base = LU()
	--base = LinearSolver()
	--base:set_convergence_check(baseConvCheck)
	--base:set_preconditioner(jac)
	
	-- Transfer and Projection
	transfer = util.CreateP1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	transfer:set_restriction_damping(1.0)
	projection = util.CreateP1Projection(approxSpace)
	
	-- Gemoetric Multi Grid
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
	--gmg:set_debug(dbgWriter)


amg = AMGPreconditioner()
--amg:enable_aggressive_coarsening_A(2)

amg:set_num_presmooth(3)
amg:set_num_postsmooth(3)
amg:set_cycle_type(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
amg:set_max_levels(5)

amg:set_max_nodes_for_base(300)
amg:set_max_fill_before_base(0.7)
amg:set_fsmoothing(true)

--amg:tostring()


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(200)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-12)

-- exact Solver
exSolver = LU()

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)
--linSolver:set_debug(dbgWriter)

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
function StartValue(x, y, t)
	if y > 0 then return y end
	return -y
end

LuaStartValue = util.CreateLuaUserNumber("StartValue", dim)

--InterpolateFunction(LuaStartValue, u, "c", 0.0)
u:set(1.0)

for i = 1,1 do
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
SaveVectorForConnectionViewer(u, "Sol.mat")

-- 4. init solver for linear Operator
print("Init solver for operator.")
if solver:init(linOp) ~= true then print("Cannot init solver"); exit() end

-- 5. apply solver
print("Apply solver.")
if solver:apply_return_defect(u,b) ~= true then print("Cannot apply solver"); end -- exit() end
end
-------------------------------------------
--  Output
-------------------------------------------
WriteGridFunctionToVTK(u, "Solution")