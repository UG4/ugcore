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

-- constants
dim = 2

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = "unit_square_tri_four_dirichlet_nodes.ugx"
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
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
	SaveDomain(dom, "refined_grid.ugx")
end

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
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
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_diffusion_tensor(diffusionMatrix)
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

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

print ("Reset initial value")
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
	SaveVectorForConnectionViewer(b, "Rhs.vec")
end

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
	baseConvCheck:set_maximum_steps(500)
	baseConvCheck:set_minimum_defect(1e-8)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose(false)
	-- base = LU()
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
	
	-- Geometric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_debug(dbgWriter)

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
hlibConvCheck:set_verbose(false)
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
