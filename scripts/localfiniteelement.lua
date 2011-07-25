-------------------------------------------------------------------------------
--
--   Lua - Script to test local finite element spaces
--
--   Author: Andreas Vogel
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2)

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- refinements:
numRefs    = util.GetParamNumber("-numRefs",    2)

-- Display parameters (or defaults):
print(" General parameters chosen:\n")
print("    dim        = " .. dim .. "\n")
print("    grid       = " .. gridName .. "\n")
print("    numRefs    = " .. numRefs .. "\n")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());
	
-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, 0, 0, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 2)
approxSpace:init()
approxSpace:print_local_dof_statistic(2)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------

function Rhs2d(x, y, t)
	return	-12.0 *x*x
end

function ExactSolution2d(x, y, t)
	return x*x*x*x
end

--diffusionMatrix = util.CreateLuaUserMatrix("DiffTensor"..dim.."d", dim)
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

rhs = util.CreateLuaUserNumber("Rhs"..dim.."d", dim)
--rhs = util.CreateConstUserNumber(0.0, dim)

exactSolution = util.CreateLuaUserNumber("ExactSolution"..dim.."d", dim)
--exactSolution = util.CreateConstBoundaryNumber(3.2, dim)
	
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

-- Select upwind
if dim == 2 then 
--upwind = NoUpwind2d()
--upwind = FullUpwind2d()
upwind = WeightedUpwind2d(); upwind:set_weight(0.0)
--upwind = PartialUpwind2d()
elseif dim == 3 then 
--upwind = NoUpwind3d()
--upwind = FullUpwind3d()
upwind = WeightedUpwind3d(); upwind:set_weight(0.0)
--upwind = PartialUpwind3d()
else print("Dim not supported for upwind"); exit() end


elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
if elemDisc:set_upwind(upwind) == false then exit() end
elemDisc:set_disc_scheme("fe")
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(exactSolution, "c", "Boundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
--domainDisc:add_elem_disc(neumannDisc)
domainDisc:add_post_process(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------

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

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_verbose_level(true)

-- create CG Solver
cgSolver = CG()
ilu = ILU()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- choose some solver
solver = cgSolver

--------------------------------------------------------------------------------
--  Apply Solver - using method defined in 'operator_util.h',
--  to get separate profiling for assemble and solve
--------------------------------------------------------------------------------

l2error = {}

for i=1,numRefs do

-- 0. refine space
local refiner = GlobalDomainRefiner(dom)
refiner:refine()
approxSpace:defragment()
approxSpace:print_statistic()

-- 1. init operator
print("Init operator (i.e. assemble matrix).")
if linOp:init_op_and_rhs(b) == false then print("Could assemble operator"); exit(); end

-- set dirichlet values in start iterate
u:set(0.0)
linOp:set_dirichlet_values(u)

InterpolateFunction(exactSolution, u, "c", 0.0)

-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.mat")
SaveVectorForConnectionViewer(u, "StartSol.mat")

-- 2. init solver for linear Operator
print("Init solver for operator.")
solver:init(linOp)

-- 3. apply solver
print("Apply solver.")
solver:apply_return_defect(u,b)

-- 4. compute error
l2error[i] = L2Error(exactSolution, u, "c", 0.0)
print("L2Error on Level "..i.." is "..l2error[i] .."\n");
end

print("L2 Error result:\n")
for i=1,numRefs do
print(i..": "..l2error[i].. "  factor: ") 
if i == 1 then print(" --- \n") 
else print(l2error[i-1]/l2error[i].."\n") end
end
--------------------------------------------------------------------------------
--  Output of computed solution
--------------------------------------------------------------------------------
WriteGridFunctionToVTK(u, "Solution")
SaveVectorForConnectionViewer(u, "Solution.mat")

