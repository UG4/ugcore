----------------------------------------------------------
--
--   Lua - 	Script to perform the Laplace-Problem on
--	 		adaptive multigrid hierarchies using error estimator
--
--   Authors: Andreas Vogel, Sebastian Reiter
--
--	This script intends to demonstrate how adaptive-multigrid
--	hierarchies with hanging nodes can be used to solve
--	the laplace problem. An error indicator is used.
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

------------------------
-- CONSTANTS
dim = util.GetParamNumber("-dim",    2)
if dim == 2 then
	gridName = "unit_square_quads_2x2.ugx"
--	gridName = "unit_square_01/unit_square_01_tri.ugx"
	gridName = "unit_square_01/unit_square_01_quads_8x8.ugx"
elseif dim == 3 then
	gridName = "open_cube_hex.ugx"
else
	print("Dimension " .. dim .. " not supported. Aborting.")
	exit()
end

--refinement (default is 5)
numRefs    = util.GetParamNumber("-numRefs",    2)

--------------------------------
-- User Data Functions (begin)
--------------------------------

-- 2D example, cf. Braess, Finite Elemente, p. 33
-- 
--	- \Delta u = 0                 in \Omega
--           u = sin(2/3 * \phi)   on circle bnd
--           u = 0                 on line bnd

function ourSource2d(x, y, t)
	local r = math.sqrt(x*x + y*y)
	return	-200*(1.0 + math.exp(-200*r+160) - 200*r + 200*math.exp(-200*r+160)*r)
				* math.exp(-200*r+160)/math.pow(1.0 + math.exp(-200*r+160), 3)/r;
end
		
function ourDirichletBnd2d(x, y, t)
	local r = math.sqrt(x*x + y*y)
	
	return true, 1.0/(1 + math.exp(-200*(r-0.8)))
end

-- 3D

function ourSource3d(x, y, z, t)
	local r = math.sqrt(x*x + y*y + z*z)
	return	-400 * (1 + math.exp(-200 * r + 160) - 100 * r + 100 * math.exp(-200*r+160)*r)
				 * math.exp(-200*r+160) * math.pow(1 + exp(-200 * r + 160), -3) * math.pow(r, -1);

end

function ourDirichletBnd3d(x, y, z, t)
	local r = math.sqrt(x*x + y*y + z*z)
	return true, 1.0/(1 + math.exp(-200*(r-0.8)))
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
--refiner = GlobalDomainRefiner(dom);

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
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Source setup
sourceValue = util.CreateLuaUserNumber("ourSource"..dim.."d", dim)

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
elemDisc:set_source(sourceValue)

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

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- debug writer
dbgWriter = util.CreateGridFunctionDebugWriter(dim)
dbgWriter:set_vtk_output(true)

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


if false then
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

amg:tostring()
end

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-12)

-- exact Solver
exSolver = LU()

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(ilu)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(jac)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ilu)
bicgstabSolver:set_convergence_check(convCheck)

-------------------------------------------
--  Apply Solver
-------------------------------------------

-- choose some solver
solver = bicgstabSolver

approxSpace:init()
approxSpace:print_statistic()
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- create surface view
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

for i = 1,30 do

print(" #######  START Adaption " .. i .."  #######");
u:set(0.0)

-- 1. init operator
print("Init operator (i.e. assemble matrix).")
if linOp:init() ~= true then print("Cannot init operator"); exit() end

-- 2. set dirichlet values in solution
print("Setting Dirichlet values.")
linOp:set_dirichlet_values(u)

-- 3. set right hand side (assembled together with Operator for performance reasons)
print("Setting rhs values.")
b:assign(linOp:get_rhs())

-- debug output
print("Saving for debug.")
SaveMatrixForConnectionViewer(u, linOp, "Stiffness".. i ..".mat")
SaveVectorForConnectionViewer(b, "Rhs".. i ..".mat")

-- 4. init solver for linear Operator
print("Init solver for operator.")
if solver:init(linOp) ~= true then print("Cannot init solver"); exit() end

-- 5. apply solver
print("Apply solver.")
if solver:apply_return_defect(u,b) ~= true then print("Cannot apply solver"); exit() end

WriteGridFunctionToVTK(u, "Solution"..i)

-- 6. estimate error and mark
MarkForRefinement_GradientIndicator(refiner, u, 1e-8, 0.5);
print("Elements marked.")

-- 7. refine
refiner:refine()
print("Grid refined.")
refiner:clear_marks()
SaveDomain(dom, "refined_grid_" .. i .. "_p" .. GetProcessRank() .. ".ugx")
SaveGridHierarchy(dom:get_grid(), "refined_grid_"..i.."hierarchy.ugx")
print("Refined Grid written to file.")

-- 8. defragment approximation space
approxSpace:defragment()
print("Defragmentate");
approxSpace:print_statistic()
print("Statistic done");

SaveGrid(dom:get_grid(), approxSpace:get_surface_view(), "surface_view_" .. i .. ".ugx")
SaveGridHierarchyTransformed(dom:get_grid(), approxSpace:get_surface_view(), "surface_view_extrude_" .. i .. "_p".. GetProcessRank() .. ".ugx", 0.25)

--if CheckSurfaceView(approxSpace:get_surface_view()) ~= true then 
--print("Surface View is not correct. Aborting."); exit(); end

-- 9. write statistik on new grid
--PrintGridElementNumbers(dom:get_grid())

print(" #######  END   Adaption " .. i .."  #######");

end
