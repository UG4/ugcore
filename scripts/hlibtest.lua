----------------------------------------------------------
--
--   Lua - Script to apply HLIBSolver to the Laplace-Problem
--
--   Author: Ingo Heppner
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

verbosity = 0	-- set to 0 i.e. for time measurements,
		-- >= 1 for writing matrix files etc.

verbosity = util.GetParamNumber("-verb", 0)

activateDbgWriter = 0	-- set to 0 i.e. for time measurements,
		        -- >= 1 for writing matrix files etc. by setting
		        -- 'fetiSolver:set_debug(dbgWriter)'
activateDbgWriter = util.GetParamNumber("-dbgw", 0)

-- choose algebra
InitAlgebra(CPUAlgebraChooser());

-- constants
dim = 2

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_tri.ugx")
	--gridName = "unit_square_tri_four_dirichlet_nodes.ugx"
	--gridName = "unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_cube_hex.ugx")
	--gridName = "unit_cube_tets_regular.ugx"
end

numPreRefs = util.GetParamNumber("-numPreRefs", 2)
numRefs    = util.GetParamNumber("-numRefs",    4)

--------------------------------
-- HLIB related
--------------------------------
if util.HasParamOption("-geom") == true then
	clustering = "geom"
else
	clustering = "alg"
end

print(" Parameters chosen:")
print("    numPreRefs =        " .. numPreRefs)
print("    numRefs    =        " .. numRefs)
print("    grid       =        " .. gridName)
print("    verbosity  =        " .. verbosity)
print("    activateDbgWriter = " .. activateDbgWriter)

print(" HLIB parameters chosen:")
print("    clustering =        " .. clustering)
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
if numPreRefs > numRefs then
	print("numPreRefs must be smaller/equal than numRefs");
	exit();
end

-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
refiner = GlobalDomainRefiner(dom)

-- Performing pre-refines
for i=1,numPreRefs do
	refiner:refine()
end

-- Distribute the domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	refiner:refine()
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
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

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

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
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
linOp:init()
print ("done")

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
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_prolongation(transfer)
	gmg:set_projection(projection)
	--gmg:set_debug(dbgWriter)

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

-- create HLIB Solver
hlibSolver = HLIBSolver()
hlibSolver:set_hlib_accuracy_H(1.e-4)  -- default: 1.e-4
hlibSolver:set_hlib_accuracy_LU(1.e-4) -- default: 1.e-4

-- define construction of cluster tree
--   first  arg: "clustering type" \in [algebraic | geometric (not yet implemented)]; algebraic is default
--   second arg: "clustering mode" \in [nested dissection | empty/everything else]; nested dissection is default 
hlibSolver:set_clustering_method("algebraic", "nested dissection")

hlibSolver:set_hlib_verbosity(4) -- '>= 2': create HLIB related postscripts; '>=3' also create plots of matrix entries
hlibSolver:set_ps_basename("hlib")

hlibConvCheck = StandardConvergenceCheck()
hlibConvCheck:set_maximum_steps(600)
hlibConvCheck:set_minimum_defect(1e-10)
hlibConvCheck:set_reduction(1e-16)
hlibConvCheck:set_verbose_level(false)
hlibSolver:set_convergence_check(hlibConvCheck)

if activateDbgWriter >= 1 then
--	print( "activating debug writer for hlibSolver")
	hlibSolver:set_debug(dbgWriter)
end

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator
--print("Init operator (i.e. assemble matrix).")
--linOp:init() - Das ist doch doppelt gemoppelt (wird oben schon ausgefuehrt)!?

-- 2. init solver for linear Operator
--print("Init solver for operator.")
--solver:init(linOp)

-- 3. apply solver
print("Apply solver.")
tBefore = os.clock()
ApplyLinearSolver(linOp, u, b, hlibSolver)
tAfter = os.clock()

-------------------------------------------
--  Output
-------------------------------------------
if verbosity >= 1 then
	WriteGridFunctionToVTK(u, "HLib-Solution")
end

if verbosity >= 2 then
	hlibSolver:check_crs_matrix()
end

-- profiling
output = io.open("hlib-profile.txt", "a")

-- get number of processes
numProcs = GetNumProcesses()

assemblePNinit  = GetProfileNode("initLinearSolver")
assemblePNapply = GetProfileNode("applyLinearSolver")

print("\n")
print("dim\tnumPreRefs\tnumRefs\tnumProcs\t\tCPU time (s)\tinitLinearSolver (ms)\tapplyLinearSolver (ms)");
s = string.format("%d\t%d\t\t%d\t\t%d\t\t\t%.2f\t\t%.2f\t\t\t%.2f\n",
		  dim, numPreRefs, numRefs,
		  numProcs,
		  tAfter-tBefore,
		  assemblePNinit:get_avg_total_time_ms(),
		  assemblePNapply:get_avg_total_time_ms())
output:write(s)
print(s)
