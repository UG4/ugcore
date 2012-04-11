----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------
GetLogAssistant():set_output_process(-1)
ug_load_script("ug_util.lua")

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
dim = util.GetParamNumber("-dim", 2)

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = "unit_square_tri_four_dirichlet_nodes.ugx"
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numRefs    = util.GetParamNumber("-numRefs",    3)

-- parallelisation related stuff
-- way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect") -- [grid2d | bisect | metis]

-- amount of output
verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: this sets
						    -- 'fetiSolver:set_debug(dbgWriter)'


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

-- choose dimension and algebra
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
		return 0
	end
	
	function ourDirichletBnd2d(x, y, t)
		return true, 0
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
-- fill a partitionMap (but every process needs one and has to return his map
-- by calling 'RedistributeDomain()', even if in this case the map is empty
-- for all processes but process 0).
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
		SavePartitionMap(partitionMap, dom, "partitionMap_p" .. GetProcessRank() .. ".ugx")
	end
end

if RedistributeDomain(dom, partitionMap, true) == false then
	print("Redistribution failed. Please check your partitionMap.")
	exit()
end
print("... domain distributed!")

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
sh = dom:subset_handler()
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
--upwind = NoUpwind()
--upwind = FullUpwind()
upwind = WeightedUpwind(); upwind:set_weight(0.0)
--upwind = PartialUpwind()

elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_upwind(upwind)
elemDisc:set_diffusion(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction_rate(reaction)
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
dirichletBND:add(dirichlet, "c", "DirichletBoundary")

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
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
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
	baseConvCheck:set_maximum_steps(1000)
	baseConvCheck:set_minimum_defect(1e-11)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose(false)
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
--convCheck:set_verbose(true)

linOp:init_op_and_rhs(b)

B = MatrixOperator()
-- domainDisc:assemble_stiffness_matrix(A, v)
v = GridFunction(approxSpace)
domainDisc:assemble_mass_matrix(B, v)
SaveMatrixForConnectionViewer(v, B, "B.mat") 
SaveMatrixForConnectionViewer(v, linOp, "A.mat")
v:set(1.0)
linOp:set_dirichlet_values(v)
SaveVectorForConnectionViewer(v, "v.mat") 

srand(1)
nev = 3

eig = EigenSolver()
eig:set_linear_operator_A(linOp)
eig:set_linear_operator_B(B)
eig:set_max_iterations(100)
eig:set_precision(1e-5)
eig:set_preconditioner(gmg)
eig:set_pinvit(3)
ev = {}
for i=1,nev do
	print("adding ev "..i)
	ev[i] = GridFunction(approxSpace)
	ev[i]:set_random(-1.0, 1.0)
	linOp:set_dirichlet_values(ev[i])
	eig:add_vector(ev[i])
end

eig:apply()

for i=1,nev do
	WriteGridFunctionToVTK(ev[i], "ev_"..i)
	SaveVectorForConnectionViewer(ev[i], "ev_"..i..".mat") 
end


