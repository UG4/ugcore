----------------------------------------------------------
--
--   Lua - Script to test FETI on the Laplace-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

-- make sure that ug_util is in the right path.
-- currently only the path in which you start your application is valid.
ug_load_script("ug_util.lua")

verbosity = 0	-- set to 0 i.e. for time measurements,
		-- >= 1 for writing matrix files etc.

verbosity = util.GetParamNumber("-verb", 0)

activateDbgWriter = 0	-- set to 0 i.e. for time measurements,
		        -- >= 1 for writing matrix files etc. by setting
		        -- 'fetiSolver:set_debug(dbgWriter)'
activateDbgWriter = util.GetParamNumber("-dbgw", 0)

-- constants
dim = 2

-- choose dimension and algebra
InitUG(dim, CPUAlgebraSelector());

if dim == 2 then
--	gridName = "rect_tri.ugx"
--	gridName = "unit_square_duplicated_tri.ugx"
	gridName = "unit_square_01/unit_square_01_tri_2x2.ugx"
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = "unit_square/unit_cube_hex.ugx"
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

numPreRefs = 2
numRefs = 4

numPreRefs = util.GetParamNumber("-numPreRefs", 2)
numRefs    = util.GetParamNumber("-numRefs",    4)

--------------------------------
-- User Data Functions (begin)
--------------------------------
	function ourDiffTensor2d(x, y, t)
		return	1, 0, 
				0, 1
	end
	
	function ourRhs2d(x, y, t)
		local s = 2*math.pi
		return	s*s*(math.sin(s*x) + math.sin(s*y))
		--return 0
		--return 0;
	end
	
	function ourDirichletBnd2d(x, y, t)
		local s = 2*math.pi
		return true, math.sin(s*x) + math.sin(s*y)
		--return true, x + 1
		--return true, 2.5
	end

--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end

-- create Refiner
print("Create Refiner")
if numPreRefs > numRefs then
	print("numPreRefs must be smaller/equal than numRefs");
	exit();
end

refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom:get_grid())
for i=1,numPreRefs do
	refiner:refine()
end

if DistributeDomain(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	util.GlobalRefineParallelDomain(dom)
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
	SaveDomain(dom, "refined_grid.ugx")
end

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()

--------------------------------------------------------------------------------
-- Gather info for domain decomposition
--------------------------------------------------------------------------------
-- get number of processes
numProcs = GetNumProcesses()
if numProcs < 2 then
	print("number of processes is smaller than 2 - huh??")
end

-- get number of processes per subdomain
numProcsPerSubdomain = util.GetParamNumber("-nPPSD", 1)

if not util.isPowerOfTwo(numProcsPerSubdomain) then
	print( "WARNING: nPPSD = '" .. numProcsPerSubdomain .. "' is not a power of 2!" )
--	return
end

print( "Check if nPPSD = '" .. numProcsPerSubdomain .. "' process(es) per subdomain makes sense ..." )

-- compute number of subdomains
numSubdomains = numProcs / numProcsPerSubdomain

-- check that numSubdomains is greater 1 && \in \N && a power of 2.
if numSubdomains < 2 then
	print( "ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is smaller than 2!? Aborting!" )
	return
end

if not util.isNaturalNumber(numSubdomains) then
	print( "ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is NOT a natural number!? Aborting!" )
	return
end

if not util.isPowerOfTwo(numSubdomains) then
	print( "WARNING: numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is not a power of 2! Continuing ..." )
-- TODO: Maybe switch to a default value then
--	return -- in this case the partition can be quite erratic (at least on small (triangular) grids)..
end

print( "NumProcs is " .. numProcs .. ", NumSubDomains is " .. numSubdomains )
--------------------------------------------------------------------------------

-- create subdomain info
print("Create domainDecompInfo")
domainDecompInfo = StandardDomainDecompositionInfo()
domainDecompInfo:set_num_subdomains(numSubdomains)

approxSpace:init()
approxSpace:print_statistic()
approxSpace:print_layout_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Diffusion Tensor setup
	diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor2d", dim)
	--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- rhs setup
	rhs = util.CreateLuaUserNumber("ourRhs2d", dim)
	--rhs = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
	dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd2d", dim)
	--dirichlet = util.CreateConstBoundaryNumber(3.2, dim)
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add(dirichlet, "c", "DirichletBoundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add(elemDisc)
--domainDisc:add(neumannDisc)
domainDisc:add(dirichletBND)

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function("u", true)
b = approxSpace:create_surface_function("b", true)

-- New creation of subdomains and layouts (since 30012011):
-- test one to many interface creation
if verbosity >= 1 then
	for i=0,GetNumProcesses()-1 do
		print("subdom of proc " .. i .. ": " .. domainDecompInfo:map_proc_id_to_subdomain_id(i))
	end
end

-- BuildDomainDecompositionLayoutsTest2d(u, domainDecompInfo);
-- OneToManyTests2d(u)


-- set initial value
u:set(0.0)

-- init Operator
print ("Assemble Operator ... ")
linOp:init_op_and_rhs(b)
print ("done")

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)

-- write matrix for test purpose
if verbosity >= 1 then
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.mat")
end

-- debug writer
dbgWriter = GridFunctionDebugWriter()
dbgWriter:set_reference_grid_function(u)
dbgWriter:set_vtk_output(false)

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)
jac2 = Jacobi()
jac2:set_damp(0.8)
jac3 = Jacobi()
jac3:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilu2 = ILU()
ilu3 = ILU()
ilut = ILUT()

-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(500)
	baseConvCheck:set_minimum_defect(1e-8)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(false)
	-- base = LapackLUSolver()
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
	
	-- Transfer and Projection
	transfer = util.CreateP1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	projection = util.CreateP1Projection(approxSpace)
	
	-- Geometric Multi Grid
	gmg = util.CreateGeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
--	gmg:set_surface_level(numRefs) -- deprecated!
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_prolongation(transfer)
	gmg:set_projection(projection)

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

-- exact Solver
exactConvCheck = StandardConvergenceCheck() -- only for debugging
exactConvCheck:set_maximum_steps(10)
exactConvCheck:set_minimum_defect(1e-10)
exactConvCheck:set_reduction(1e-16)
exactConvCheck:set_verbose_level(false)
exactSolver = LU()
exactSolver:set_convergence_check(exactConvCheck)

-- create coarse problem Solver
cpConvCheck = StandardConvergenceCheck()
cpConvCheck:set_maximum_steps(2000)
cpConvCheck:set_minimum_defect(1e-10)
cpConvCheck:set_reduction(1e-16)
cpConvCheck:set_verbose_level(false)
coarseproblemSolver = CG()
coarseproblemSolver:set_preconditioner(ilu) -- Absturz: "Signal: Bus error (10)", "Signal code:  (5583)", etc ...
--coarseproblemSolver:set_preconditioner(jac) -- "ERROR in 'JacobiPreconditioner::apply': Cannot change parallel status of correction to consistent."
--coarseproblemSolver:set_preconditioner(gs) -- "ERROR in 'Gauss-Seidel::apply': Cannot change parallel storage type of correction to consistent."
-- Ganz ohne preconditioner: "Cannot convert z to consistent vector."
coarseproblemSolver:set_convergence_check(cpConvCheck)

-- create Neumann Solver
neumannConvCheck = StandardConvergenceCheck()
neumannConvCheck:set_maximum_steps(2000)
neumannConvCheck:set_minimum_defect(1e-10)
neumannConvCheck:set_reduction(1e-16)
neumannConvCheck:set_verbose_level(false)
neumannSolver = CG()
--neumannSolver = BiCGStab()
neumannSolver:set_preconditioner(ilu)
neumannSolver:set_convergence_check(neumannConvCheck)

-- create Dirichlet Solver
dirichletConvCheck = StandardConvergenceCheck()
dirichletConvCheck:set_maximum_steps(2000)
dirichletConvCheck:set_minimum_defect(1e-10)
dirichletConvCheck:set_reduction(1e-16)
dirichletConvCheck:set_verbose_level(false)
dirichletSolver = CG()
dirichletSolver:set_preconditioner(ilu2)
dirichletSolver:set_convergence_check(dirichletConvCheck)

-- create FETI Solver
fetiSolver = FETI()

fetiConvCheck = StandardConvergenceCheck()
fetiConvCheck:set_maximum_steps(100)
fetiConvCheck:set_minimum_defect(1e-7)
fetiConvCheck:set_reduction(1e-16)

fetiSolver:set_convergence_check(fetiConvCheck)
fetiSolver:set_domain_decomp_info(domainDecompInfo)

fetiSolver:set_neumann_solver(neumannSolver)
fetiSolver:set_dirichlet_solver(dirichletSolver)
fetiSolver:set_coarse_problem_solver(exactSolver)
--fetiSolver:set_coarse_problem_solver(coarseproblemSolver)
if activateDbgWriter >= 1 then
	fetiSolver:set_debug(dbgWriter)
end

-- Apply Solver
print( "   numPreRefs is " .. numPreRefs .. ",  numRefs is " .. numRefs)
print( "   numProcs   is " .. numProcs   .. ",  NumSubDomains is " .. numSubdomains )

tBefore = os.clock()
ApplyLinearSolver(linOp, u, b, fetiSolver)
tAfter = os.clock()

-- Output
output = io.open("feti-profile.txt", "a")

assemblePNinit  = GetProfileNode("initLinearSolver")
assemblePNapply = GetProfileNode("applyLinearSolver")

print("\n")
print("dim\tnumPreRefs\tnumRefs\tnumProcs\tnumSD\tnumProcsPerSD\t\tCPU time (s)\tinitLinearSolver (ms)\tapplyLinearSolver (ms)");
s = string.format("%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t\t%.2f\t\t%.2f\t\t\t%.2f\n",
		  dim, numPreRefs, numRefs,
		  numProcs, numSubdomains, numProcsPerSubdomain,
		  tAfter-tBefore,
		  assemblePNinit:get_avg_total_time_ms(),
		  assemblePNapply:get_avg_total_time_ms())
output:write(s)
print(s)

fetiSolver:print_statistic_of_inner_solver()

if verbosity >= 1 then
	WriteGridFunctionToVTK(u, "Solution")
end

-- check
print( "domainDecompInfo:get_num_subdomains(): " .. domainDecompInfo:get_num_subdomains())
