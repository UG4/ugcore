----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------
PrintBuildConfiguration()

ug_load_script("ug_util.lua")

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
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
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    3)

if numPreRefs > numRefs then
	print("It must be choosen: numPreRefs <= numRefs");
	exit();
end

-- parallelisation related stuff
-- way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect") -- [grid2d | bisect | metis]

-- amount of output
verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: call 'set_debug(dbgWriter)'
						    -- for the main solver ('gmg')


-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));


--------------------------------
-- User Data Functions (begin)
--------------------------------
	function ourDiffTensor2d(x, y, t)
		return	1, 0, 
				0, 1
	end
	
	function ourVelocityField2d(x, y, t)
		return	0, 0
	end
	
	function ourReaction2d(x, y, t)
		return	0
	end
	
	function ourRhs2d(x, y, t)
		local s = 2*math.pi
		return	s*s*(math.sin(s*x) + math.sin(s*y))
		--return -2*y
--		return 0;
	end
	
	function ourNeumannBnd2d(x, y, t)
		--local s = 2*math.pi
		--return -s*math.cos(s*x)
		return true, -2*x*y
	end
	
	function ourDirichletBnd2d(x, y, t)
		local s = 2*math.pi
		return true, math.sin(s*x) + math.sin(s*y)
		--return true, x*x*y
--		return true, 2.5
	end

	function ourDiffTensor3d(x, y, z, t)
		return	1, 0, 0,
				0, 1, 0,
				0, 0, 1
	end
	
	function ourVelocityField3d(x, y, z, t)
		return	0, 0, 0
	end
	
	function ourReaction3d(x, y, z, t)
		return	0
	end
	
	function ourRhs3d(x, y, z, t)
		--local s = 2*math.pi
		--return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
		return 0;
	end
	
	function ourNeumannBnd3d(x, y, t)
		--local s = 2*math.pi
		--return -s*math.cos(s*x)
		return true, -2*x*y*z
	end
	
	function ourDirichletBnd3d(x, y, z, t)
		--local s = 2*math.pi
		--return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
		return true, x
	end
--------------------------------
-- User Data Functions (end)
--------------------------------

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
-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
refiner = GlobalDomainRefiner(dom)

-- Performing pre-refines
print("Perform (non parallel) pre-refinements of grid")
for i=1,numPreRefs do
	write( "PreRefinement step " .. i .. " ...")
	refiner:refine()
	print( " done.")
end

-- get number of processes
numProcs = GetNumProcesses()

-- Distribute the domain to all involved processes
-- Since only process 0 loaded the grid, it is the only one which has to
-- fill a partitionMap (but every process needs one and has to return his map
-- by calling 'RedistributeDomain()', even if in this case the map is empty
-- for all processes but process 0).
if numProcs > 1 then
	print("Distribute domain with 'distributionType' = " .. distributionType .. "...")
	partitionMap = PartitionMap()
	
	if GetProcessRank() == 0 then
		if distributionType == "bisect" then
			util.PartitionMapBisection(dom, partitionMap, numProcs)
			
		elseif distributionType == "grid2d" then
			local numNodesX, numNodesY = util.FactorizeInPowersOfTwo(numProcs / numProcsPerNode)
			util.PartitionMapLexicographic2D(dom, partitionMap, numNodesX,
											 numNodesY, numProcsPerNode)
	
		elseif distributionType == "metis" then
			util.PartitionMapMetis(dom, partitionMap, numProcs)
											 
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
	
	print("Redistribute domain with 'distributionType' = '" .. distributionType .. "' ...")
	if RedistributeDomain(dom, partitionMap, true) == false then
		print("Redistribution failed. Please check your partitionMap.")
		exit()
	end
	print("... domain distributed!")
	delete(partitionMap)
end


--------------------------------------------------------------------------------
-- end of partitioning
--------------------------------------------------------------------------------

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	write( "Refinement step " .. i .. " ...")
	refiner:refine()
	print( " done!")
end


-- Make sure, that the required subsets are present
requiredSubsets = {"Inner", "Boundary"}
if util.CheckSubsets(dom, requiredSubsets) == false then 
   print("Subsets missing. Aborting")
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
	if SaveGridHierarchy(dom:grid(), hierarchyOutName) == false then
		print("Saving of domain to " .. hierarchyOutName .. " failed. Aborting.")
		    exit()
	end
end

print("NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_level()
approxSpace:init_surface()
approxSpace:print_local_dof_statistic(2)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------
-------------------------------------------
print ("Setting up Assembling")

-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDiffTensor2d")
-- Diffusion Tensor setup
diffusionMatrix = LuaUserMatrix("ourDiffTensor"..dim.."d")
--diffusionMatrix = ConstUserMatrix(1.0)

-- Velocity Field setup
velocityField = LuaUserVector("ourVelocityField"..dim.."d")
--velocityField = ConstUserVector(0.0)

-- Reaction setup
reaction = LuaUserNumber("ourReaction"..dim.."d")
--reaction = ConstUserNumber(0.0)

-- rhs setup
rhs = LuaUserNumber("ourRhs"..dim.."d")
--rhs = ConstUserNumber(0.0)

-- neumann setup
neumann = LuaBoundaryNumber("ourNeumannBnd"..dim.."d")
--neumann = ConstUserNumber(0.0)

-- dirichlet setup
dirichlet = LuaBoundaryNumber("ourDirichletBnd"..dim.."d")
--dirichlet = ConstBoundaryNumber(3.2)
	
--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

-- Select upwind
if dim == 2 then 
--upwind = NoUpwind2d()
--upwind = FullUpwind2d()
upwind = WeightedUpwind(); upwind:set_weight(0.0)
--upwind = PartialUpwind2d()
elseif dim == 3 then 
--upwind = NoUpwind3d()
--upwind = FullUpwind3d()
upwind = WeightedUpwind3d(); upwind:set_weight(0.0)
--upwind = PartialUpwind3d()
else print("Dim not supported for upwind"); exit() end


elemDisc = ConvectionDiffusion("c", "Inner")
if elemDisc:set_upwind(upwind) == false then exit() end
elemDisc:set_disc_scheme("fv1")
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = FV1NeumannBoundary("Inner")
--neumannDisc:add(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
--domainDisc:add(neumannDisc)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator(domainDisc)

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(true)

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
	baseConvCheck:set_maximum_steps(1000)
	baseConvCheck:set_minimum_defect(1e-11)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(false)
	--base = LU()
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)

	-- Geometric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_parallel_base_solver(true)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
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
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_verbose_level(true)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create ILU Solver
iluSolver = LinearSolver()
iluSolver:set_preconditioner(jac)
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
--  Apply Solver - using method defined in 'operator_util.h',
--  to get separate profiling for assemble and solve
--------------------------------------------------------------------------------

-- 0. Reset start solution
u:set(0.0)

-- 1. init operator
print("Init operator (i.e. assemble matrix).")
if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
	print("Could not assemble operator"); exit(); 
end

-- 1.b write matrix for test purpose
if verbosity >= 1 then
	print("Save matrix and rhs for connection viewer ...")
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.vec")
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
