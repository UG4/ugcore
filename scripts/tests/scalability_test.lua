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

-- some parameters for renaming of logfiles (in the moment only used for FETI)
str_problem = "laplace"
str_startgrid = ""

-- output command line arguments
line = "#ANALYZER INFO: cmd args:"
for _, arg in ipairs(ugargv) do
	line = line .. " " .. arg
end
print(line)

if PclDebugBarrierEnabled() then
	print("#ANALYZER INFO: PCLDebugBarrier is enabled. Expect slowdowns.")
end

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
dim = util.GetParamNumber("-dim", 2)

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_8x8.ugx")
	str_startgrid = "8x8-quads"
	--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_2x2x2.ugx")
	str_startgrid = "2x2x2-quads"
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numRefs    = util.GetParamNumber("-numRefs",    3)

-- parameters concerning the linear solver:
lsType     = util.GetParam("-lsType",         "gmg") -- choose one in ["gmg" | "feti" | "hlib"]
lsIterator = util.GetParam("-lsIterator",     "gmg")
lsMaxIter  = util.GetParamNumber("-lsMaxIter", 100)

-- MG solver related parameters
if lsIterator == "gmg" then
	baseSolverType = util.GetParam("-bs", "exact") -- choose one in ["exact" | "cg"  | "jac"]
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

-- amount of output
verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: this sets
						    -- 'fetiSolver:set_debug(dbgWriter)'



renameLogfileAfterRun = 0	  
renameLogfileAfterRun = util.GetParamNumber("-rlf", 0)
logfileName = "" -- empty at start; will be built by concatenating some parameters

--------------------------------------------------------------------------------
-- HIERARCHICAL REDISTRIBUTION
-- Note that hierarchical redistribution is not compatible with
-- distributionType == "grid2d"
-- If this is set to a level > numPreRefs and <= numRefs then hierarchical
-- redistribution is enabled. If it is set to -1, it is disabled (default).
-- Defines the level in which hierarchical redistribution starts.
-- Note that hierarchical redistribution is only performed, if a refinement follows.
hRedistFirstLevel = util.GetParamNumber("-hRedistFirstLevel", -1)

-- The number of processes to which a grid is redistributed in a
-- hierarchical redistribution step.
hRedistNewProcsPerStep = util.GetParamNumber("-hRedistNewProcsPerStep", math.pow(2, dim))

-- defines the number of refinements between hierarchical redistributions
-- if < 1 the hierarchical redistribution is disabled.
-- After s steps, the grid will be redistributed to
-- hRedistNewProcsPerStep processes.
-- Only has effect, if hierarchicalRedistFirstLevel ~= -1.
hRedistStepSize = util.GetParamNumber("-hRedistStepSize", 1)

-- the number of processes to which we will distribute the grid during
-- initial distribution has to be calculated now.
-- It depends on whether we perform hierarchical redistribution or not.
numInitialDistProcs = GetNumProcesses()
hRedistEnabled = false
if hRedistFirstLevel ~= -1 then
--	make sure that all parameters are valid
	if hRedistFirstLevel <= numPreRefs or hRedistFirstLevel >= numRefs then
		print("HORIZONTAL REDISTRIBUTION ERROR:")
		print("  Make sure that numPreRefs < hRedistFirstLevel < numRefs")
		exit()
	end

--	make sure that distribution type grid2d is not used (the util-script is
--	not suited for hierarchical redistribution)
	if distributionType == "grid2d" then
		print("HORIZONTAL REDISTRIBUTION ERROR:")
		print("  distributionType 'grid2d' is not supported for hierarchical"
			  .. " redistribution.")
		exit()
	end
	
	hRedistEnabled = true
	
	local procs = GetNumProcesses()
	local refinements = numRefs - hRedistFirstLevel
	while refinements > 0 do
		refinements = refinements - hRedistStepSize
		if procs / hRedistNewProcsPerStep < 1 then
			break
		else
			procs = math.floor(procs / hRedistNewProcsPerStep)
		end
	end
	
	numInitialDistProcs = procs
	print("#ANALYZER INFO: hierarchical redistribution: initial distribution at level "
			.. numPreRefs .. " to " .. procs .. " processes.")
elseif hierarchicalRedistFirstLevel ~= -1 then
	print("WARNING: hierarchical redistribution disabled. To enable it, "
		  .. "set hierarchicalRedistFirstLevel to a value > numPreRefs and <= numRefs.")
end


--------------------------------------------------------------------------------
-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Linear solver related parameters chosen:")
print("    lsType     = " .. lsType .. " (TMP - should be integrated better ...)")
print("    lsIterator = " .. lsIterator)
print("    lsMaxIter  = " .. lsMaxIter)

if lsIterator == "gmg" then
	print("    bs (baseSolverType) = " .. baseSolverType)
	print("    bl (baseLevel)      = " .. baseLevel)
end

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)
print("    numPPN (numProcsPerNode) = " .. numProcsPerNode)
print("    hRedistFirstLevel = " .. hRedistFirstLevel)
print("    hRedistNewProcsPerStep = " .. hRedistNewProcsPerStep)
print("    hRedistStepSize = " .. hRedistStepSize)

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------


-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
if LoadDomain(dom, gridName) == false then
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
	write( "PreRefinement step " .. i .. " ...")
	refiner:refine()
	print( "  done.")
end


-- redistribution... the optional hierarchical approach makes it somewhat complicated.
numProcs = GetNumProcesses()
numDistProcs = numInitialDistProcs
numProcsWithGrid = 1
numCurRefs = numPreRefs

-- Distribute the domain to all involved processes
-- Only processes which already have a grid will fill their partition maps.
partitionMap = PartitionMap()

while numDistProcs > 0 do
	partitionMap:clear()
	if GetProcessRank() < numProcsWithGrid then
	--	add target procs. Make sure to keep a portion on the local process
		partitionMap:add_target_proc(GetProcessRank())
		if numDistProcs > 1 then
			local firstNewProc = numProcsWithGrid + GetProcessRank() * (numDistProcs - 1)
			partitionMap:add_target_procs(firstNewProc, numDistProcs - 1)
		end
	
		if distributionType == "bisect" then
			PartitionDomain_Bisection(dom, partitionMap, 0)
	
		elseif distributionType == "grid2d" then
		--	Note that grid2d can not be used if hierarchical redistribution is active.
			local numNodesX, numNodesY = util.FactorizeInPowersOfTwo(numDistProcs / numProcsPerNode)
			util.PartitionMapLexicographic2D(dom, partitionMap, numNodesX,
											 numNodesY, numProcsPerNode)
	
		elseif distributionType == "metis" then
			PartitionDomain_MetisKWay(dom, partitionMap, numProcs, baseLevel, 1, 1)
	
		else
		    print( "distributionType not known, aborting!")
		    exit()
		end

	-- save the partition map for debug purposes
		if verbosity >= 1 then
			print("saving partition map to 'partitionMap_p" .. GetProcessRank() .. ".ugx'")
			SavePartitionMap(partitionMap, dom, "partitionMap_p" .. GetProcessRank() .. ".ugx")
		end
	end
	
	print("Distribute domain with 'distributionType' = '" .. distributionType .. "' ...")
	if RedistributeDomain(dom, partitionMap, true) == false then
		print("Redistribution failed. Please check your partitionMap.")
		exit()
	end
	print("  ... domain distributed.")

	numProcsWithGrid = numProcsWithGrid * numDistProcs
	numDistProcs = 0
	
	-- check whether we have to perform another hierarchical distribution.
	-- calculate the number of required refinements in this step.
	maxRefsInThisStep = numRefs
	if hRedistEnabled == true then
		if numCurRefs < hRedistFirstLevel then
			maxRefsInThisStep = hRedistFirstLevel
		else
			maxRefsInThisStep = numCurRefs + hRedistStepSize
		end
		if maxRefsInThisStep >= numRefs then
			numDistProcs = 0
			maxRefsInThisStep = numRefs
		else
			numDistProcs = hRedistNewProcsPerStep
			if numProcsWithGrid * numDistProcs > numProcs then
				numDistProcs = 0
				maxRefsInThisStep = numRefs
			end
		end
	end
	
	-- Perform post-refine
	print("Refine Parallel Grid")
	for i = numCurRefs + 1, maxRefsInThisStep do
		print( "Refinement step " .. i .. " ...")
		refiner:refine()
		print( "  ... done.")
	end
	numCurRefs = maxRefsInThisStep
end

-- clean up
delete(partitionMap)

--------------------------------------------------------------------------------
-- end of partitioning
--------------------------------------------------------------------------------

-- get subset handler
sh = dom:get_subset_handler()
subsetsFine = (sh:num_subsets() == 2)
if subsetsFine == true then subsetsFine = util.CheckSubsets(dom, {"Inner", "Boundary"}) end
if AllProcsTrue(subsetsFine) == false then 
	print("Domain must have 2 Subsets for this problem ('Inner' and 'Boundary').")
	print("Make sure that every process received a part of the grid during distribution!")
	exit()
end

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

print("#ANALYZER INFO: NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
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
diffusionMatrix = ConstUserMatrix(1.0)

-- The right hand side is given by a lua function
rhs = LuaUserNumber("ourRhs"..dim.."d")

-- The Dirichlet values are also given by lua values
dirichlet = LuaBoundaryNumber("ourDirichletBnd"..dim.."d")
	
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

-- Select upwind
if dim == 2 then 
--upwind = NoUpwind2d()
--upwind = FullUpwind2d()
upwind = WeightedUpwind2d(); upwind:set_weight(0.0)
--upwind = PartialUpwind2d()
elseif dim == 3 then upwind = NoUpwind3d()
else print("Dim not supported for upwind"); exit() end


elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
if elemDisc:set_upwind(upwind) == false then exit() end
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
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
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- debug writer
dbgWriter = GridFunctionDebugWriter()
dbgWriter:set_reference_grid_function(u)
dbgWriter:set_vtk_output(true) -- TMP 06102011, war false

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
	baseConvCheck:set_minimum_defect(1e-12)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(false)

-- choose base solver
isBaseSolverParallel = false
if baseSolverType == "exact" then
	base = LU()
elseif baseSolverType == "jac" then
	base = LinearSolver()
	baseJac = Jacobi()
	base:set_preconditioner(baseJac)
elseif baseSolverType == "cg" then
	base = CG() 
	base:set_preconditioner(jac)
	isBaseSolverParallel = true
else
	print ("base solver not specified ==> exit")
	exit()
end
	base:set_convergence_check(baseConvCheck)
	
	-- Transfer and Projection
	transfer = P1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	projection = P1Projection(approxSpace)
	
	-- Geometric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(baseLevel) -- variable defining baselevel, set from script
	gmg:set_base_solver(base)
	gmg:set_parallel_base_solver(isBaseSolverParallel)
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

if lsType == "feti" then
	print("Loading FETI solver setup ...")
	ug_load_script("setup_fetisolver.lua")

	solver, logfileName = SetupFETISolver(str_problem,str_startgrid,
					      dim,
					      lsMaxIter,
					      numProcs,
					      activateDbgWriter,
					      dbgWriter,
					      verbosity, logfileName)

elseif lsType == "hlib" then
	print("Loading HLIB solver setup ...")
	ug_load_script("setup_hlibsolver.lua")
	solver = SetupHLIBSolver(lsMaxIter,
				 activateDbgWriter,
				 dbgWriter,
				 verbosity)
end

--------------------------------------------------------------------------------
--  Apply Solver - using method defined in 'operator_util.h',
--  to get separate profiling for assemble and solve
--------------------------------------------------------------------------------

if lsType == "gmg" and lsIterator == "gmg" then
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
if verbosity >= 2 then
	print("Save matrix and rhs for connection viewer ...")
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.vec")
end

-- 2. apply solver
print("Apply solver.")
if ApplyLinearSolver(linOp, u, b, solver) == false then
	print("Could not apply linear solver.");
end

if lsType == "feti" then
	solver:print_statistic_of_inner_solver()
	if verbosity >= 2 then
		--print("Testing standard interfaces ...")
		--TestDomainInterfaces(dom)

		print("# Testing FETI layouts ...")
		solver:test_layouts(false) -- 'true' prints indices
	end
	if renameLogfileAfterRun then
		print("Renaming logfile (name built in 'setup_fetisolver.lua') to '" .. logfileName .. "' ...")
		if GetLogAssistant():rename_log_file(logfileName) == true then
			print(".. done!")
		else
			print(".. could not rename - no logfile open! Try again with '-logtofile <name>'!")
		end
	end
end

--------------------------------------------------------------------------------
--  Output of computed solution
--------------------------------------------------------------------------------
if verbosity >= 1 then
	print("Write vtk file for solution ...")
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
