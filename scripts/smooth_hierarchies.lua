----------------------------------------------------------
--
--   Lua - 	Script to perform the Laplace-Problem on
--	 		smooth multigrid hierarchies.
--
--   Authors: Andreas Vogel, Sebastian Reiter
--
----------------------------------------------------------

-- filenames are interpreted relative to APP_PATH/../scripts,
-- where APP_PATH is the path in which ug_shell is located.
-- If the script is not found, then the method will attempt to
-- load it relative to the current path.
ug_load_script("ug_util.lua")

-- constants
dim = 2

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

gridName = "volgrid_b.ugx"
requiredSubsets = {"Inner", "Boundary"}
numRefs = 4
useLoopScheme = true


-- user functions
function DiffTensor(x, y, t)
	return	1, 0, 
			0, 1
end

function VelocityField(x, y, t)
	return	0, 0
end

function Reaction(x, y, t)
	return	0
end

function Rhs(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

function DirichletBnd(x, y, t)
	local s = 2*math.pi
	return	true, math.sin(s*x) + math.sin(s*y)
end


-- Diffusion Tensor setup
diffusionMatrix = LuaUserMatrix("DiffTensor")
-- Velocity Field setup
velocityField = LuaUserVector("VelocityField")
-- Reaction setup
reaction = LuaUserNumber("Reaction")
-- rhs setup
rhs = LuaUserNumber("Rhs")
-- dirichlet
dirichlet = LuaBoundaryNumber("DirichletBnd")

--------------------------------------------------------------------------------
-- Domain preparation
--------------------------------------------------------------------------------
-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
if LoadDomain(dom, gridName) == false then
	print("LoadDomain failed. Make sure the specified file exists!")
	exit()
end

-- create Refiner
print("Create Hierarchy")
if useLoopScheme == true then
	CreateSmoothHierarchy(dom:grid(), numRefs)
else
	CreateSemiSmoothHierarchy(dom:grid(), numRefs)
end

print("Saving domain for debug reasons")
-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")
SaveGridHierarchy(dom:grid(), "refined_grid_hierarchy.ugx")


--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------
print("creating approximation space...")
-- create Approximation Space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-- make sure that the required subsets exist	
if util.CheckSubsets(dom, requiredSubsets) == false then 
   print("Subsets missing. Aborting")
   exit()
end

-- create Discretization
print("creating element discretization...")
upwind = WeightedUpwind(); upwind:set_weight(0.0)
elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_upwind(upwind)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction_rate(reaction)
elemDisc:set_source(rhs)

print("creating dirichlet boundary...")
-- create dirichlet boundary and add it to discretization
dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary")


print("creating domain discretization...")
-- The domain discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

print("filling stiffness matrix")	
-- create the stiffness matrix
A = MatrixOperator()
-- create the grid functions (vectors for the solution and rhs)
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- assemble matrix and rhs
domainDisc:assemble_linear(A, b)

-- set dirichlet values
u:set(0)
domainDisc:adjust_solution(u)

-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")


--------------------------------------------------------------------------------
-- Solvers
--------------------------------------------------------------------------------
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
baseConvCheck:set_verbose(false)

-- base = LU()
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_surface_level(numRefs)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CGSolver()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- Apply Solver
ApplyLinearSolver(A, u, b, linSolver)

-- Output
WriteGridFunctionToVTK(u, "Solution")
