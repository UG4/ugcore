-------------------------------------------------------------------------------
--
--	Jump-Diffusion in 3d
--
--	This script is used to set up a 3d testing problem as it has been 
-- 	implemented in ug3.9 with the following setting:
-- 		- The problem domain is the unit cube
--		- Diffusion equation with
--		  a) zero right-hand side (no source)
--		  b) diffusion coeff diagonal, 
--				D == 10^-6	if 1/3 <= |x| <= 3/4
--			    D == 1	   	else
--		- Neumann zero boundary at front, back, upper and bottom side
--		- Dirichlet 1.0 value at left side
--		- Dirichlet 0.0 value at right side
--
--   Author: Andreas Vogel
--
-------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 3

if dim == 3 then
	gridName = util.GetParam("-grid", "grids/unit_cube_2x2x2_jump_diff.ugx")
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    4)
numLoop	   = util.GetParamNumber("-numLoop",    1)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

--------------------------------------------------------------------------------
--  Domain setup
--------------------------------------------------------------------------------

-- Lets define a list of all subsets that we need
neededSubsets = {"Inner", "Top", "Bottom", "Front", "Back", "Left", "Right"}

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- write grid to file for test purpose
print("saving refined domain to refined_grid.ugx")
SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------

function jumpDiffTensor3d(x, y, z, t)

	local r = math.sqrt(x*x+y*y+z*z)
	local d = 1.0
	--if r >= 1/3 and r <= 3/4 then d = 1e-6 end

	return	d, 0, 0,
			0, d, 0,
			0, 0, d
end

-- Diffusion Tensor setup
luaDiffusionMatrix = util.CreateLuaUserMatrix("jumpDiffTensor"..dim.."d", dim)
--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- dirichlet setup
constLeftDirichlet = util.CreateConstBoundaryNumber(1.0, dim)
constRightDirichlet = util.CreateConstBoundaryNumber(0.0, dim)
	
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(luaDiffusionMatrix)

--------------------------------------------------------------------------------
--  Setup Dirichlet Boundary
--------------------------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add(constLeftDirichlet, "c", "Left")
dirichletBND:add(constRightDirichlet, "c", "Right")

--------------------------------------------------------------------------------
--  Setup Domain Discretization
--------------------------------------------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

-- debug writer
dbgWriter = util.CreateGridFunctionDebugWriter(dim)
dbgWriter:set_reference_grid_function(u)
dbgWriter:set_vtk_output(false)

-- create algebraic Preconditioner
jac = Jacobi()
--jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
--ilu:set_debug(dbgWriter)
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
	
	-- Gemoetric Multi Grid
	gmg = util.CreateGeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(ilu)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(1)
	gmg:set_num_postsmooth(1)
	gmg:set_prolongation(transfer)
	gmg:set_projection(projection)
	--gmg:set_debug(dbgWriter)

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create ILU Solver
iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(jac)
bicgstabSolver:set_convergence_check(convCheck)

-- choose some solver
solver = linSolver

--------------------------------------------------------------------------------
--  Apply Solver
--------------------------------------------------------------------------------

for i=1,numLoop do
	
	print(" ###### LOOP " .. i .. " #######")
	
	-- 0. set initial value
	u:set(0.0)
	
	-- 1. init operator
	print("Init operator (i.e. assemble matrix).")
	tAssembleStart = os.clock() 
	if linOp:init_op_and_rhs(b) == false then print("Could assemble operator"); exit(); end
	tAssembleEnd = os.clock()
	print("Assembling took " .. tAssembleEnd - tAssembleStart .. " seconds.");
	
	-- set dirichlet values in start iterate
	linOp:set_dirichlet_values(u)
	
	-- write matrix for test purpose
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.mat")
	
	-- 2. init solver for linear Operator
	print("Init solver for operator.")
	tInitStart = os.clock() 
	if solver:init(linOp) == false then print("Could not init solver"); exit(); end
	tInitEnd = os.clock()
	print("Solver Init took " .. tInitEnd - tInitStart .. " seconds.");
	
	-- 3. apply solver
	print("Apply solver.")
	
	-- Start some energy consumption measurement
	--ExecuteSystemCommand("echo \"Start script\"")
	
	tSolveStart = os.clock() 
	if solver:apply_return_defect(u,b) == false then print("Solver failed."); exit(); end
	tSolveEnd = os.clock()

	-- Stop some energy consumption measurement
	--ExecuteSystemCommand("echo \"Stop script\"")
	
	print("Solving took " .. tSolveEnd - tSolveStart .. " seconds.");
	
end

--------------------------------------------------------------------------------
--  Output
--------------------------------------------------------------------------------

WriteGridFunctionToVTK(u, "Solution")