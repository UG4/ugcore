----------------------------------------------------------
--
--   Lua - Script for scalability studies solving the
--   Laplace-Problem.
--
--   For some flexibility several definitions concerning
--   the solver can also be given via command line parameters
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: this sets
						    -- 'fetiSolver:set_debug(dbgWriter)'

----------------------------------
-- Checking for parameters (begin)
----------------------------------
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

-- way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "grid2d") -- [grid2d | bisect]

-- parameters concerning the linear solver:
lsIterator = util.GetParam("-lsIterator",     "gmg")
lsMaxIter  = util.GetParamNumber("-lsMaxIter", 100)

-- MG solver related parameters
if lsIterator == "gmg" then
	baseSolverType = util.GetParam("-bs", "exact") -- choose one in ["exact" | "cg" ]
	baseLevel      = util.GetParamNumber("-bl",    0)
else
	baseSolverType = "exact" -- no meaning here, just to leave gmg stuff untouched
	baseLevel      = numRefs -- no meaning here, just to leave gmg stuff untouched
end

print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

print("    distType   = " .. distributionType)

print("    verbosity  =        " .. verbosity)
print("    activateDbgWriter = " .. activateDbgWriter)

print(" Linear solver related parameters chosen:")
print("    lsIterator = " .. lsIterator)
print("    lsMaxIter  = " .. lsMaxIter)

print("    baseSolverType = " .. baseSolverType)
print("    baseLevel      = " .. baseLevel)

----------------------------------
-- Checking for parameters (end)
----------------------------------


-- choose algebra
InitAlgebra(CPUAlgebraSelector());


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
		--return 0;
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
		--return true, 2.5
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
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end

-- create Refiner
print("Create Refiner")
if numPreRefs >= numRefs then
	print("numPreRefs must be smaller than numRefs");
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

--------------------------------------------------------------------------------
-- New lexicographic partitioning (cf. 'tutorials/tut06_flexible_domain_distribution.lua',
-- first stage; 04052011ih)
-- Execute on cekon via e.g.
-- salloc -n  8 mpirun ./ugshell  -ex ../scripts/mylaplace2.lua -dim 2 -grid unit_square_quads_8x8.ugx -numPreRefs 0 -numRefs 8
-- With "second stage":
-- salloc -n 32 mpirun ./ugshell  -ex ../scripts/mylaplace2.lua -dim 2 -grid unit_square_quads_8x8.ugx -numPreRefs 0 -numRefs 8
-- Execute on JuGene:
-- mpirun -np   32 -exe ./ugshell -mode VN -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args "-ex ../scripts/laplace_jacobi-fixed-it.lua -outproc 0 -grid unit_square_quads_8x8.ugx -numPreRefs 4 -numRefs  8"

--------------------------------------------------------------------------------
-- Distribute the domain to all involved processes
if distributionType == "bisect" then
	if DistributeDomain(dom) == false then
		print("Error while Distributing Grid.")
		exit()
	end
-- it would be nice to also have the partition map as output in this case for comparison ...
--	-- We'll save the partition map. This should only be done for debugging.
--	if verbosity >= 1 then
--		SaveGrid(dom:get_grid(), partitionMap:get_partition_handler(), -- how to get these parameters?
--			"partitionMap_bisect.ugx")
--	end
elseif distributionType == "grid2d" then

--	ug_load_script("lexmap_processes.lua")
	ug_load_script("lexmap_processes2.lua")
	print(" lexicographic right up partitioning done!")

else

        print( "distributionType not known, aborting!")
        exit()

end

--------------------------------------------------------------------------------
-- Ende Partitionierung (04052011ih)
--------------------------------------------------------------------------------
--exit() -- TMP break

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

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDiffTensor2d")
-- Diffusion Tensor setup
diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
velocityField = util.CreateLuaUserVector("ourVelocityField"..dim.."d", dim)
--velocityField = util.CreateConstUserVector(0.0, dim)

-- Reaction setup
reaction = util.CreateLuaUserNumber("ourReaction"..dim.."d", dim)
--reaction = util.CreateConstUserNumber(0.0, dim)

-- rhs setup
rhs = util.CreateLuaUserNumber("ourRhs"..dim.."d", dim)
--rhs = util.CreateConstUserNumber(0.0, dim)

-- neumann setup
neumann = util.CreateLuaBoundaryNumber("ourNeumannBnd"..dim.."d", dim)
--neumann = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
--dirichlet = util.CreateConstBoundaryNumber(3.2, dim)
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

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
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = util.CreateNeumannBoundary(approxSpace, "Inner")
--neumannDisc:add_boundary_value(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(dirichlet, "c", "DirichletBoundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
--domainDisc:add_elem_disc(neumannDisc)
domainDisc:add_post_process(dirichletBND)

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:export_rhs(true)
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

print ("Reset initial value")
-- set initial value
u:set(0.0)

-- init Operator
print ("Assemble Operator ... ")
tAssembleStart = os.clock() 
if linOp:init() == false then print("Could assemble operator"); exit(); end
tAssembleEnd = os.clock()
print("Assembling took " .. tAssembleEnd - tAssembleStart .. " seconds.");

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)
b:assign(linOp:get_rhs())

-- write matrix for test purpose
if verbosity >= 1 then
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.mat")
end

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

	-- base = LapackLUSolver()
	base = LinearSolver()
	base:set_preconditioner(jac)

elseif baseSolverType == "cg" then

	base = CG() -- parallel solver on baselevel for testing two grid method (29032011ih) - or as parallel solver on baselevel 'numPreRefs+1' (13042011ih)
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
	gmg:set_base_level(baseLevel) -- new variable defining baselevel, set from script (default 0; 14042011ih)
--        gmg:set_base_level(numRefs-1) -- two grid method for testing (29032011ih)
	gmg:set_base_solver(base)
ilusmoother = ILU()
	gmg:set_smoother(jac)
--	gmg:set_smoother(ilusmoother) -- jac (29932911ih)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_prolongation(transfer)
	gmg:set_projection(projection)
	if activateDbgWriter >= 1 then
		gmg:set_debug(dbgWriter)
	end

-- create AMG ---
-----------------

	if false then
	amg = AMGPreconditioner()
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
convCheck:set_maximum_steps(lsMaxIter) -- until 18042011: literal '100'
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_verbose_level(true) -- added 28032011ih

-- create Linear Solver
linSolver = LinearSolver()
if lsIterator == "gmg" then
	linSolver:set_preconditioner(gmg) -- until 09062011: literal 'gmg'
elseif lsIterator == "jac" then
	linSolver:set_preconditioner(jac) -- until 09062011: literal 'gmg'
else
	print ("linear solver iterator not specified ==> exit")
	exit()
end
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
--solver = cgSolver

-------------------------------------------
--  Apply Solver - using method defined in 'operator_util.h', to get separate profiling for assemble and solve (05042011ih)
-------------------------------------------
-- 1. init operator
--print("Init operator (i.e. assemble matrix).")
--linOp:init() -- really necessary? Already done above!

-- 2. init solver for linear Operator
--print("Init solver for operator.")
--solver:init(linOp)

-- 3. apply solver
--print("Apply solver.")
--solver:apply_return_defect(u,b)

print( "   baseSolverType is " .. baseSolverType .. ",  baseLevel is " .. baseLevel )

print("Init operator and solver, then apply solver.")
tBefore = os.clock()
ApplyLinearSolver(linOp, u, b, solver)
tAfter = os.clock()

-------------------------------------------
--  Output
-------------------------------------------
if verbosity >= 1 then
	WriteGridFunctionToVTK(u, "Solution")
end
