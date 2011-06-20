--------------------------------------------------------------------------------
--
--   Lua - Script for scalability studies solving the
--   Laplace-Problem.
--
--   For some flexibility several definitions concerning
--   the solver can also be given via command line parameters
--
--   Authors: Ingo Heppner, Sebastian Reiter, Andreas Vogel
--
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- Execute on cekon via e.g.
-- salloc -n 32 mpirun ./ugshell  -ex ../scripts/tests/scalability_test.lua -dim 2 -grid unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs 8

-- Execute on JuGene via e.g. ('mpirun' call has to be specified in a LoadLeveler script!):
-- mpirun -np 32 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs  8"
-- Or the same job interactively:
-- llrun -v -np 32 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid  unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs 8"

--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: this sets
						    -- 'fetiSolver:set_debug(dbgWriter)'

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
dim = util.GetParamNumber("-dim", 2)

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end

if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numRefs    = util.GetParamNumber("-numRefs",    3)

-- parameters concerning the linear solver:
lsIterator = util.GetParam("-lsIterator",     "gmg")
lsMaxIter  = util.GetParamNumber("-lsMaxIter", 100)

-- MG solver related parameters
if lsIterator == "gmg" then
	baseSolverType = util.GetParam("-bs", "exact") -- choose one in ["exact" | "cg"  |Ê"jac"]
	baseLevel      = util.GetParamNumber("-bl",    0)
else
	baseSolverType = "exact" -- no meaning here, just to leave gmg stuff untouched
	baseLevel      = numRefs -- no meaning here, just to leave gmg stuff untouched
end

-- parallelisation related stuff
-- way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect") -- [grid2d | bisect | metis]

-- number of processes per node (only used if distType == grid2d)
-- should be a square number
numProcsPerNode = util.GetParamNumber("-numPPN", 1)


-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Linear solver related parameters chosen:")
print("    lsIterator = " .. lsIterator)
print("    lsMaxIter  = " .. lsMaxIter)

if lsIterator == "gmg" then
	print("    bs (baseSolverType) = " .. baseSolverType)
	print("    bl (baseLevel)      = " .. baseLevel)
end

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)
print("    numPPN (numProcsPerNode) = " .. numProcsPerNode)

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------


-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- create Instance of a Domain
print("Create Domain.")
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end

-- create Refiner
print("Create Refiner")
if numPreRefs > numRefs then
	print("It must be choosen: numPreRefs <= numRefs");
	exit();
end

-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
refiner = GlobalDomainRefiner(dom)

-- Performing pre-refines
print("Perform (non parallel) pre-refinements of grid")
for i=1,numPreRefs do
	print( "PreRefinement step " .. i .. " ...")
	refiner:refine()
	print( "... done!")
end

-- get number of processes
numProcs = GetNumProcesses()

-- Distribute the domain to all involved processes
-- Since only process 0 loaded the grid, it is the only one which has to
-- fill a partitionMap.
partitionMap = PartitionMap()

if GetProcessRank() == 0 then
	if distributionType == "bisect" then
		util.PartitionMapBisection(partitionMap, dom, numProcs)
		
	elseif distributionType == "grid2d" then
		local numNodesX, numNodesY = util.FactorizeInPowersOfTwo(numProcs / numProcsPerNode)
		util.PartitionMapLexicographic2D(partitionMap, dom, numNodesX,
										 numNodesY, numProcsPerNode)

	elseif distributionType == "metis" then
		util.PartitionMapMetis(partitionMap, dom, numProcs)
										 
	else
	    print( "distributionType not known, aborting!")
	    exit()
	end
	
-- save the partition map for debug purposes
	if verbosity >= 1 then
		SavePartitionMap(partitionMap, dom, "partitionMap_p" .. GetProcessRank() .. ".ugx")
	end
end

if RedistributeDomain(dom, partitionMap, true) == false then
	print("Redistribution failed. Please check your partitionMap.")
	exit()
end

--------------------------------------------------------------------------------
-- end of partitioning
--------------------------------------------------------------------------------

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	print( "Refinement step " .. i .. " ...")
	refiner:refine()
	print( "... done!")
end

-- get subset handler
sh = dom:get_subset_handler()
if sh:num_subsets() ~= 2 then 
print("Domain must have 2 Subsets for this problem.")
exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- write grid to file for test purpose
if verbosity >= 1 then
	refinedGridOutName = "refined_grid_p" .. GetProcessRank() .. ".ugx"
	print("saving domain to " .. refinedGridOutName)
	if SaveDomain(dom, refinedGridOutName) == false then
		print("Saving of domain to " .. refinedGridOutName .. " failed. Aborting.")
		    exit()
	end
	
	hierarchyOutName = "hierachy_p" .. GetProcessRank() .. ".ugx"
	print("saving hierachy to " .. hierarchyOutName)
	if SaveGridHierarchy(dom:get_grid(), hierarchyOutName) == false then
		print("Saving of domain to " .. hierarchyOutName .. " failed. Aborting.")
		    exit()
	end
end

print("NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee
if OrderCuthillMcKee(approxSpace, true) == false then
	print("ERROR when ordering Cuthill-McKee"); exit();
end

--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------

function ourRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
	--return -2*y
	--return 0;
end

function ourDirichletBnd2d(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
	--return true, x*x*y
	--return true, 2.5
end

function ourRhs3d(x, y, z, t)
	--local s = 2*math.pi
	--return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
	return 0;
end

function ourDirichletBnd3d(x, y, z, t)
	--local s = 2*math.pi
	--return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
	return true, x
end

-- We need a constant Identity matrix for the Diffusion.
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- The right hand side is given by a lua function
rhs = util.CreateLuaUserNumber("ourRhs"..dim.."d", dim)

-- The Dirichlet values are also given by lua values
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
	
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

-- Select upwind 
if dim == 2 then 
--upwind = NoUpwind2d()
--upwind = FullUpwind2d()
upwind = WeightedUpwind2d(); upwind:set_weight(0.0)
--upwind = PartialUpwind2d()
else print("Dim not supported for upwind"); exit() end


elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
if elemDisc:set_upwind(upwind) == false then exit() end
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(dirichlet, "c", "DirichletBoundary")

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
domainDisc:add_post_process(dirichletBND)

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
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilu:set_debug(dbgWriter)
ilut = ILUT()

-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(500)
	baseConvCheck:set_minimum_defect(1e-8)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(false)

-- choose base solver
if baseSolverType == "exact" then
	base = LU()
elseif baseSolverType == "jac" then
	base = LinearSolver()
	baseJac = Jacobi()
	base:set_preconditioner(baseJac)
elseif baseSolverType == "cg" then
	base = CG() 
	base:set_preconditioner(jac)
else
	print ("base solver not specified ==> exit")
	exit()
end
	base:set_convergence_check(baseConvCheck)
	
	-- Transfer and Projection
	transfer = util.CreateP1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	projection = util.CreateP1Projection(approxSpace)
	
	-- Gemoetric Multi Grid
	gmg = util.CreateGeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(baseLevel) -- variable defining baselevel, set from script
	gmg:set_base_solver(base)
	gmg:set_parallel_base_solver(false)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_prolongation(transfer)
	gmg:set_projection(projection)
	if activateDbgWriter >= 1 then
		gmg:set_debug(dbgWriter)
	end

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(lsMaxIter)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_verbose_level(true)

-- create Linear Solver
linSolver = LinearSolver()
if lsIterator == "gmg" then 	linSolver:set_preconditioner(gmg) 
elseif lsIterator == "jac" then	linSolver:set_preconditioner(jac) 
else print ("linear solver iterator not specified ==> exit"); exit(); end

linSolver:set_convergence_check(convCheck)

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
--solver = cgSolver

--------------------------------------------------------------------------------
--  Apply Solver - using method defined in 'operator_util.h',
--  to get separate profiling for assemble and solve
--------------------------------------------------------------------------------

if lsIterator == "gmg" then
	print( "   baseSolverType is " .. baseSolverType .. ",  baseLevel is " .. baseLevel )
end

-- 0. Reset start solution
u:set(0.0)

-- 1. init operator
print("Init operator (i.e. assemble matrix).")
if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
	print("Could not assemble operator"); exit(); 
end

-- 1.b write matrix for test purpose
if verbosity >= 1 then
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.mat")
end

-- 2. apply solver
print("Apply solver.")
if ApplyLinearSolver(linOp, u, b, solver) == false then
	print("Could not apply linear solver.");
end

--------------------------------------------------------------------------------
--  Output of computed solution
--------------------------------------------------------------------------------
if verbosity >= 1 then
	WriteGridFunctionToVTK(u, "Solution")
end

--------------------------------------------------------------------------------
--  Print Profiling
--------------------------------------------------------------------------------

-- check if profiler is available
if GetProfilerAvailable() == true then
    -- get node
    pn = GetProfileNode("main")
--    pn2 = GetProfileNode("GMG_lmgc")
    -- check if node is valid
    if pn:is_valid() then
	    print(pn:call_tree())
--        print(pn2:total_time_sorted())
    else
        print("main is not known to the profiler.")
    end
else
    print("Profiler not available.")
end 
