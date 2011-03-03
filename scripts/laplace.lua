----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

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

numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numRefs    = util.GetParamNumber("-numRefs",    3)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

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
--	which automatically creates a parallel refiner if required.
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
SaveDomain(dom, "refined_grid.ugx")

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
elemDisc:set_rhs(rhs)

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
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.mat")

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
	
	-- Gemoetric Multi Grid
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

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
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

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator
print("Init operator (i.e. assemble matrix).")
linOp:init()

-- 2. init solver for linear Operator
print("Init solver for operator.")
solver:init(linOp)

-- 3. apply solver
print("Apply solver.")
solver:apply_return_defect(u,b)

-------------------------------------------
--  Output
-------------------------------------------
WriteGridFunctionToVTK(u, "Solution")