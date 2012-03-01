-------------------------------------------------------------------------------
--
--   Lua - Script to test local finite element spaces
--
--   Author: Andreas Vogel
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = util.GetParamNumber("-dim", 2)
numRefs    = util.GetParamNumber("-numRefs",    0)
numLoop    = util.GetParamNumber("-numLoop",    0)
trialOrder = util.GetParamNumber("-order",    1)
discType   = util.GetParam("-type", "fe")
elemType   = util.GetParam("-elem", "simplex")

if dim == 2 then
	if elemType == "simplex" then 
	  gridName = "unit_square_01/unit_square_01_tri_2x2.ugx"
	elseif elemType == "cube" then
	  gridName = "unit_square_01/unit_square_01_quads_2x2.ugx"
	else print("elemType not found"); exit(); end
end
if dim == 3 then
	if elemType == "simplex" then 
	  	gridName = "unit_square_01/unit_cube_01_tets.ugx"
	elseif elemType == "cube" then
	  	gridName = "unit_square_01/unit_cube_01_hex_1x1x1.ugx"
	elseif elemType == "prism" then
	 	gridName = "unit_square_01/unit_cube_01_prism_2x2x2.ugx"
	else print("elemType not found"); exit(); end
end

-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numLoop    = " .. numLoop)
print("    order      = " .. trialOrder)
print("    type       = " .. discType)
print("    elem       = " .. elemType)

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));
	
-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", trialOrder)
approxSpace:init_level()
approxSpace:init_surface()
approxSpace:print_local_dof_statistic(2)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------

function Rhs2d(x, y, t)
--	return -12.0 *x*x
	return 10*10*(math.sin(10*x) + math.sin(10*y))
--	return 0
end

function ExactSolution2d(x, y, t)
--	return x*x*x*x
	return math.sin(10*x) + math.sin(10*y)
--	return x*y
end

function Rhs3d(x, y, z, t)
--	return	-12.0 *x*x
	return 0
end

function ExactSolution3d(x, y, z, t)
--	return x*x*x*x
	return x
end

--diffusionMatrix = LuaUserMatrix("DiffTensor"..dim.."d")
diffusionMatrix = ConstUserMatrix(1.0)

rhs = LuaUserNumber("Rhs"..dim.."d")
--rhs = ConstUserNumber(0.0)

exactSolution = LuaUserNumber("ExactSolution"..dim.."d")
--exactSolution = ConstBoundaryNumber(3.2)
	
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

-- Select upwind
--upwind = NoUpwind()
--upwind = FullUpwind()
upwind = WeightedUpwind(); upwind:set_weight(0.0)
--upwind = PartialUpwind()

elemDisc = ConvectionDiffusion("c", "Inner")
if discType ~= "fe" and discType ~= "fv" then print("type of discretization not found"); exit(); end

if elemDisc:set_upwind(upwind) == false then exit() end
elemDisc:set_disc_scheme(discType)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(exactSolution, "c", "Boundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------

-- create operator from discretization
linOp = AssembledLinearOperator(domainDisc)

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_verbose_level(true)

-- create CG Solver
solver = BiCGStab()
ilu = ILU()
solver:set_preconditioner(ilu)
solver:set_convergence_check(convCheck)

--------------------------------------------------------------------------------
--  Apply Solver - using method defined in 'operator_util.h',
--  to get separate profiling for assemble and solve
--------------------------------------------------------------------------------

if numLoop == 0 then

	write("Assemble Operator...")
	if linOp:init_op_and_rhs(b) == false then print("Could assemble operator"); exit(); end
	write("done.\n")

	write("Set Start Solution...")
	u:set(0.0)
	linOp:set_dirichlet_values(u)
	write("done.\n")
		
	--InterpolateFunction(exactSolution, u, "c", 0.0)
	
	-- write matrix for test purpose
	write("Write Matrix/Vector for debug ...")
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
--	SaveVectorForConnectionViewer(b, "Rhs.vec")
	SaveVectorForConnectionViewer(u, "StartSol.vec")
	write("done.\n")

	write("Init solver...")
	solver:init(linOp)
	write("done.\n")
	
	-- 3. apply solver
	write("Apply solver...")
	solver:apply_return_defect(u,b)
	write("done.\n")
	
	-- 4. compute error
	error = L2Error(exactSolution, u, "c", 0.0, "Inner")
	write("L2Error is "..error  .."\n");
	
	WriteGridFunctionToVTK(u, "Solution")
	SaveVectorForConnectionViewer(u, "Solution.vec")
	
	exit()
end




l2error = {}

for i=1, numLoop do
	
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
	
	--InterpolateFunction(exactSolution, u, "c", 0.0)
	
	-- write matrix for test purpose
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.vec")
	SaveVectorForConnectionViewer(u, "StartSol.vec")
	
	-- 2. init solver for linear Operator
	print("Init solver for operator.")
	solver:init(linOp)
	
	-- 3. apply solver
	print("Apply solver.")
	solver:apply_return_defect(u,b)
	
	-- 4. compute error
	l2error[i] = L2Error(exactSolution, u, "c", 0.0)
	write("L2Error on Level "..i.." is "..l2error[i] .."\n");
end

print("L2 Error result:\n")
for i=1,numLoop do
	write(i..": "..l2error[i].. "  factor: ") 
	if i == 1 then write(" --- \n") 
	else write(l2error[i-1]/l2error[i].."\n") end
end
--------------------------------------------------------------------------------
--  Output of computed solution
--------------------------------------------------------------------------------
WriteGridFunctionToVTK(u, "Solution")
SaveVectorForConnectionViewer(u, "Solution.vec")

