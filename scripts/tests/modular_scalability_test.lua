--------------------------------------------------------------------------------
--
--   Lua - Script for scalability studies solving the Laplace-Problem.
--
--   For some flexibility several definitions concerning the solver
--   can also be given via command line parameters
--
--   Authors: Ingo Heppner, Sebastian Reiter, Andreas Vogel
--
--------------------------------------------------------------------------------
-- Some usage infos:
-- * For an automatic (re)naming of logfile after simulation run add the
--   following parameters / options:
--
--      '-logtofile bla -rlf'
--
--   'bla' is a temporary dummy name, '-rlf' ("Rename Log File") eventually
--   does the job. The new name is (gradually through the script) assembled
--   by concatenating some significant parameter values, e.g. problem,
--   dimension, startgrid, refinements, solver ...
--   Note: Not tested on JuGene (there a logfile is created by default and
--   named by the load leveler)!

--------------------------------------------------------------------------------
-- Execute on cekon via e.g.
-- salloc -n 32 mpirun ./ugshell  -ex ../scripts/tests/scalability_test.lua -dim 2 -grid unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs 8

-- Execute on JuGene via e.g. ('mpirun' call has to be specified in a LoadLeveler script!):
-- mpirun -np 32 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs  8"
-- Or the same job interactively:
-- llrun -v -np 32 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid  unit_square/unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs 8"

--------------------------------------------------------------------------------
print("Execution of script 'modular_scalability_test.lua' started at " .. os.date())
PrintBuildConfiguration()



ug_load_script("ug_util.lua")
ug_load_script("domain_distribution_util.lua")

-- get number of processes of this run
numProcs = GetNumProcesses()

-- some parameters for renaming of logfiles (in the moment only used for FETI)
str_problem = "scaltest-laplace"
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
-- choose dimension
dim = util.GetParamNumber("-dim", 2)

-- choose grid
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
-- choose number of pre-Refinements (before sending grid onto different processes)	
numPreRefs = util.GetParamNumber("-numPreRefs", 1)

-- choose number of total Refinements (incl. pre-Refinements)
numRefs    = util.GetParamNumber("-numRefs",    3)

if numPreRefs > numRefs then
	print("It must be choosen: numPreRefs <= numRefs");
	exit();
end


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

-- amount of output
verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.
printSol  = util.GetParamNumber("-ps",   0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for VTK output of solution

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: call 'set_debug(dbgWriter)'
						    -- for the main solver ('gmg', 'fetiSolver')

if util.HasParamOption("-rlf") then
	renameLogfileAfterRun = true
else
	renameLogfileAfterRun = false
end
logfileName = "" -- empty at start; will be built by concatenating some relevant parameters


-- Here all parameters related to refinement and distribution are parsed.
-- -numPreRefs, -numRefs, -distType, -numPPN, -hRedistFirstLevel,
-- -numRedistNewProcsPerStep, -hRedistStepSize.
-- Check the documentation of domain_distribution_util.lua for more information.
-- distribution related parameters
if(ddu.ParseAndInitializeParameters(dim) == false) then
	print("An error occured during ddu.ParseAndInitializeParameters. Aborting")
	exit()
end

--------------------------------------------------------------------------------
-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numPreRefs = " .. numPreRefs)
print("    numRefs    = " .. numRefs)

print("    verbosity (verb)         = " .. verbosity)
print("    ps (printSol)            = " .. printSol)
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
print("    numProcs   = " .. numProcs)
ddu.PrintParameters("    ") -- displays also 'numRefs' and 'numPreRefs', already displayed above!

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1))

-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
if LoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end
--------------------------------------------------------------------------------
-- refine, partition and distribute the domain.
-- Parameters were initialized above.
ddu.RefineAndDistributeDomain(dom, verbosity)
--------------------------------------------------------------------------------
-- End of partitioning
--------------------------------------------------------------------------------

-- get subset handler -- TODO: koennte - und sollte dann m.E. auch - Check fuer Subsets frueher erfolgen? Oder muss das nach dem Verteilen erfolgen?
sh = dom:subset_handler()
subsetsFine = (sh:num_subsets() == 2)
if subsetsFine == true then subsetsFine = util.CheckSubsets(dom, {"Inner", "Boundary"}) end

if AllProcsTrue(subsetsFine) == false then 
	print("Domain must have 2 Subsets for this problem ('Inner' and 'Boundary').")
	print("Make sure that every process received a part of the grid during distribution!")
	exit()
end
--------------------------------------------------------------------------------

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
	if SaveGridHierarchy(dom:grid(), hierarchyOutName) == false then
		print("Saving of domain to " .. hierarchyOutName .. " failed. Aborting.")
		    exit()
	end
end
--------------------------------------------------------------------------------

ddu.PrintAnalyzerInfo()
print("#ANALYZER INFO: grid = " .. gridName)



-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
--approxSpace:init_levels() 
approxSpace:init_top_surface() -- init top surface for DoF statistic before execution of solver
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee (last arg: flag if "reverse Cuthill-McKee" is used)
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------
-------------------------------------------
print ("Setting up Assembling")

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


elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_upwind(upwind)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
print("Adding bnd conds to global problem.")
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

-- create algebraic Preconditioner
print("Creating Preconditioner.")
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
if activateDbgWriter >= 1 then
	ilu:set_debug(dbgWriter)
end
ilut = ILUT()

-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(500)
	baseConvCheck:set_minimum_defect(1e-12)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose(false)

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
	
	-- Geometric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(baseLevel) -- variable defining baselevel, set from parameters
	gmg:set_base_solver(base)
	gmg:set_parallel_base_solver(isBaseSolverParallel)
	gmg:set_smoother(jac)
	gmg_gamma = 1
	gmg_nu1   = 3
	gmg_nu2   = 3
	gmg:set_cycle_type(gmg_gamma)
	gmg:set_num_presmooth(gmg_nu1)
	gmg:set_num_postsmooth(gmg_nu2)
	if activateDbgWriter >= 1 then
		gmg:set_debug(dbgWriter)
	end

-- create AMG ---
-----------------

	if false then
	amg = RSAMGPreconditioner()
	amg:set_nu1(2)
	amg:set_nu2(2)
	amg:set_gamma(1)
	amg:set_presmoother(jac)
	amg:set_postsmoother(jac)
	amg:set_base_solver(base)
	--amg:set_debug(u)
	end

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(lsMaxIter)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_verbose(true)

-- create Linear Solver
linSolver = LinearSolver()
if lsIterator == "gmg" then 	linSolver:set_preconditioner(gmg) 
elseif lsIterator == "jac" then	linSolver:set_preconditioner(jac) 
else print ("linear solver iterator not specified ==> exit"); exit(); end

linSolver:set_convergence_check(convCheck) -- will be overwritten if 'lsType \in ["feti" | "hlib"]

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

-- Create new name of logfile (will be used in 'scalability_test.lua', 'GetLogAssistant:rename_log_file()')
logpostfix    = ".txt" -- ".log"

-- concatenate general part of logfile name
str_problem = str_problem .. "-" .. dim .. "d"
logfileName_tmp = str_problem .. "_" .. str_startgrid

str_refs = "refs-" .. ddu.numPreRefs .. "-" .. ddu.numRefs
logfileName_tmp = logfileName_tmp .. "_" .. str_refs
logfileName_tmp = logfileName_tmp .. "_" .. lsType

if lsType == "feti" then
	print("Loading FETI solver setup ('setup_fetisolver.lua') ...")
	ug_load_script("setup_fetisolver.lua")

	solver, lsConvCheck, logfileName = SetupFETISolver(str_problem,
					      dim,
					      lsMaxIter,
					      numProcs,
					      u,                         -- for testvector writer for FAMG (created by 'CreateAMGTestvector()')
					      dirichletBND, approxSpace, -- for testvector writer for FAMG (created by 'CreateAMGTestvectorDirichlet0()')
					      activateDbgWriter,
					      verbosity, logfileName_tmp)

elseif lsType == "hlib" then
	print("Loading HLIB solver setup ('setup_hlibsolver.lua') ...")
	ug_load_script("setup_hlibsolver.lua")
	solver, lsConvCheck, logfileName = SetupHLIBSolver(str_problem,
					      dim,
					      lsMaxIter,
					      activateDbgWriter,
					      dbgWriter,
					      verbosity, logfileName_tmp)

	-- maybe one would like to improve naming in this case ...
	logfileName = logfileName_tmp
else -- (lsType == "gmg")
	lsConvCheck = convCheck
	-- maybe one would like to improve naming in this case ...
	cycle_type  = "-gamma-" .. gmg_gamma .. "-nu1-" .. gmg_nu1 .. "-nu2-" .. gmg_nu2
	logfileName = logfileName_tmp .. "_" .. lsType .. cycle_type
end
print ("Setup of Algebra Solver finished!")

str_pe = "pe-" .. numProcs
logfileName = logfileName .. "_" .. str_pe
logfileName = logfileName .. logpostfix
-- Creation of logfile name END

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
print("Init operator (i.e. assemble matrix) ...")
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
--approxSpace:print_layout_statistic() -- TMP
--approxSpace:print_statistic() -- TMP

if lsType == "feti" then
	solver:print_statistic_of_inner_solver(false) -- true: print only averages
	print("Averages for subproblem solvers:")
	solver:print_statistic_of_inner_solver(true) -- true: print only averages
	print("")
	if verbosity >= 2 then
		--print("Testing standard interfaces ...")
		--TestDomainInterfaces(dom)

		print("# Testing FETI layouts ...")
		solver:test_layouts(false) -- 'true' prints indices
	end
end

--lsConvCheck = solver:convergence_check() -- this doesn't work; returns only an object of type of the base class ('IConvergenceCheck') ...
--SetDebugShell(?) -- no registered method til now (21032012)
--breakpoint() -- available only in serial case
--error("break before printSol")
print("#ANALYZER INFO: linear solver converged in " .. lsConvCheck:step() .. " steps (defect reached: " .. lsConvCheck:defect() ..")")

--------------------------------------------------------------------------------
--  Output of computed solution
--------------------------------------------------------------------------------
if printSol >= 1 then
	print("Write vtk file for solution ...")
	WriteGridFunctionToVTK(u, "Solution")
end

--------------------------------------------------------------------------------
--  Print Profiling
--------------------------------------------------------------------------------

-- check if profiler is available
if GetProfilerAvailable() == true then
    print("")
    -- get node
    pn = GetProfileNode("main")
--    pn2 = GetProfileNode("GMG_lmgc")
    -- check if node is valid
    if pn:is_valid() then
	    print(pn:call_tree(0.0))
--        print(pn2:total_time_sorted())
    else
        print("main is not known to the profiler.")
    end
else
    print("Profiler not available.")
end 

--------------------------------------------------------------------------------
--  Rename logfile
--------------------------------------------------------------------------------
if renameLogfileAfterRun == true then
	print("Renaming logfile to '" .. logfileName .. "' ...")
	if GetLogAssistant():rename_log_file(logfileName) == true then
		print(".. done!")
	else
		print(".. could not rename - no logfile open! Try again with '-logtofile <name>'!")
	end
else
	print("To automatically rename logfile to '" .. logfileName .. "' add option '-rlf'!")
end

if util.HasParamOption("-stats") then
	if GetProcessRank() == 0 then
		stats = {			
			{ "procs", GetNumProcesses() },
			{ "numPreRefs", numPreRefs},
			{ "numRefs", numRefs },
			{ "lastReduction", convCheck:defect()/convCheck:previous_defect()},
			{ "steps", convCheck:step()},
			{ "main [ms]", GetProfileNode("main"):get_avg_total_time_ms()},
			{ "perform_refinement [ms]", GetProfileNode("perform_refinement"):get_avg_total_time_ms()},
			{ "AssembleLinearOperatorRhsAndSolution [ms]", GetProfileNode("ASS_AssembleLinearOperatorRhsAndSolution"):get_avg_total_time_ms()},
			{ "InitLinearSolver [ms]", GetProfileNode("ALS_InitLinearSolver"):get_avg_total_time_ms()},
			{ "ApplyLinearSolver [ms]", GetProfileNode("ALS_ApplyLinearSolver"):get_avg_total_time_ms()},
			{ "lsType", lsType},
			{ "lsIterator", lsIterator},
			{ "lsMaxIter", lsMaxIter},
			{ "dim", dim},
			{ "gridName", gridName},
			{ "baseSolverType", baseSolverType},
			{ "baseLevel", baseLevel},
			{ "date", os.date("y%Ym%md%d") },
			{ "SVN Revision", GetSVNRevision()},
			{"host",GetBuildHostname()},
			{"logfileName", logfileName},
			{"commandline", util.GetCommandLine() }
		} 
			
		util.printStats(stats)
		statsdir = util.GetParam("-stats", ".")
		util.writeFileStats(stats, statsdir.."/modular_scalability_test.txt")

		fromfile=NULL
		if renameLogfileAfterRun == true then
			fromfile = logfileName				
		elseif util.HasParamOption("-logtofile") then
			fromfile = util.GetParam("-logtofile")
		end
		if fromfile then
			t=io.open(fromfile, "r"):read("*all")
			io.open(statsdir.."/"..logfileName, "w"):write(t)				
		end						
	end 
end

print("Execution finished at " .. os.date())
