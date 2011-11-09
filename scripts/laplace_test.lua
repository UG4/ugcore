----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- constants
dim = 2

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = "unit_square_tri_four_dirichlet_nodes.ugx"
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    3)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

-- Init UG for dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

--------------------------------------------------------------------------------
-- Domain Setup
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)
test.require(dom ~= nil, "Domain loaded.")

-- create Approximation Space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee
test.require(OrderCuthillMcKee(approxSpace, true), 
	"ordering Cuthill-McKee");

--------------------------------------------------------------------------------
-- User Data Setup
--------------------------------------------------------------------------------
	
-- Diffusion Tensor setup
diffusionMatrix = ConstUserMatrix(1.0)

function userRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

function userRhs3d(x, y, z, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
end

-- rhs setup
rhs = LuaUserNumber("userRhs"..dim.."d")

function userDirichletBnd2d(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
end
		
function userDirichletBnd3d(x, y, z, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
end

-- dirichlet setup
dirichlet = LuaBoundaryNumber("userDirichletBnd"..dim.."d")

--------------------------------------------------------------------------------
-- FV Disc setup
--------------------------------------------------------------------------------

elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_upwind(NoUpwind())
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------

-- matrix and vectors
A = MatrixOperator()
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)

-- Base Solver
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-8)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose_level(false)

base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

-- Transfer and Projection
transfer = P1Prolongation(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)
projection = P1Projection(approxSpace)

-- Gemoetric Multi Grid
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_prolongation(transfer)
gmg:set_projection(projection)


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- choose some solver
solver = linSolver

-------------------------------------------
--  Apply Solver
-------------------------------------------

-- 1. init operator
print("Init operator (i.e. assemble matrix).")
domainDisc:assemble_linear(A, b)

-- 2. set initial value
print ("Reset initial value")
u:set(0.0)
domainDisc:adjust_solution(u)

-- 3. init solver for linear Operator
print("Init solver for operator.")
test.require(solver:init(A), "Initializing Solver.")

-- 4. apply solver
print("Apply solver.")
test.require(solver:apply_return_defect(u,b), "Applying Solver.")

-------------------------------------------
--  Output
-------------------------------------------
WriteGridFunctionToVTK(u, "Solution")