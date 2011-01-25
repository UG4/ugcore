----------------------------------------------------------
--
--   Example - Script to test a 2d reticulum setup
--
--   Author: Markus Breit
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraChooser());

-- Dimension
dim = 2

-- Grid
gridName = "simple_reticulum.ugx"

-- Refinements before distributing grid
numPreRefs = 0

-- Total refinements
numRefs = 2

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
		--return	s*s*(math.sin(s*x) + math.sin(s*y))
		--return -2*y
		return 0;
	end
	
	function ourNeumannBnd2d(x, y, t)
		--local s = 2*math.pi
		--return -s*math.cos(s*x)
		return true, -1.0
	end
	
	function ourDirichletBnd2d(x, y, t)
		local s = 2*math.pi
		--return true, math.sin(s*x) + math.sin(s*y)
		return true, x
		--return true, 2.5
	end

	function membraneDirichletBnd2d(x, y, t)
	local s = 2*math.pi
	--return true, math.sin(s*x) + math.sin(s*y)
	--return true, x*x*y
	return true, 2.5
end
	
--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom = utilCreateDomain(dim)

-- load domain
print("Load Domain from File.")
if utilLoadDomain(dom, gridName) == false then
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

if utilDistributeDomain(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	utilGlobalRefineParallelDomain(dom)
end

-- get subset handler
sh = dom:get_subset_handler()
if sh:num_subsets() ~= 4 then 
print("Domain must have 4 Subsets for this problem.")
exit()
end

-- write grid to file for test purpose
utilSaveDomain(dom, "refined_grid.ugx")

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1FunctionOnSubsets(pattern, "c_er", "er, mem_er", dim)
AddP1FunctionOnSubsets(pattern, "c_cyt", "cyt, mem_er, mem_cyt", dim)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Diffusion Tensor setup
	diffusionMatrix = utilCreateLuaUserMatrix("ourDiffTensor2d", dim)
	--diffusionMatrix = utilCreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
	velocityField = utilCreateLuaUserVector("ourVelocityField2d", dim)
	--velocityField = utilCreateConstUserVector(0.0, dim)

-- Reaction setup
	reaction = utilCreateLuaUserNumber("ourReaction2d", dim)
	--reaction = utilCreateConstUserNumber(0.0, dim)

-- rhs setup
	rhs = utilCreateLuaUserNumber("ourRhs2d", dim)
	--rhs = utilCreateConstUserNumber(0.0, dim)

-- neumann setup
	neumann = utilCreateLuaBoundaryNumber("ourNeumannBnd2d", dim)
	--neumann = utilCreateConstUserNumber(0.0, dim)

-- dirichlet setup
	dirichlet = utilCreateLuaBoundaryNumber("ourDirichletBnd2d", dim)
	--dirichlet = utilCreateConstBoundaryNumber(3.2, dim)
	
-- dirichlet setup
	membraneDirichlet = utilCreateLuaBoundaryNumber("membraneDirichletBnd2d", dim)


-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDiscER = utilCreateFV1ConvDiff(approxSpace, "c_er", "er")
elemDiscER:set_upwind_amount(0.0)
elemDiscER:set_diffusion_tensor(diffusionMatrix)
elemDiscER:set_velocity_field(velocityField)
elemDiscER:set_reaction(reaction)
elemDiscER:set_rhs(rhs)

elemDiscCYT = utilCreateFV1ConvDiff(approxSpace, "c_cyt", "cyt")
elemDiscCYT:set_upwind_amount(0.0)
elemDiscCYT:set_diffusion_tensor(diffusionMatrix)
elemDiscCYT:set_velocity_field(velocityField)
elemDiscCYT:set_reaction(reaction)
elemDiscCYT:set_rhs(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDisc = utilCreateNeumannBoundary(approxSpace, "cyt")
neumannDisc:add_boundary_value(neumann, "c_cyt", "mem_cyt")

neumannDisc2 = utilCreateNeumannBoundary(approxSpace, "er")
neumannDisc2:add_boundary_value(neumann, "c_er", "mem_er")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

--dirichletBND = utilCreateDirichletBoundary(approxSpace)
--dirichletBND:add_boundary_value(dirichlet, "c", "Boundary, MembraneBnd")

--membraneDirichletBND = utilCreateDirichletBoundary(approxSpace)
--membraneDirichletBND:add_boundary_value(membraneDirichlet, "c_membrane", "MembraneBnd")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDiscER)
domainDisc:add_elem_disc(elemDiscCYT)
domainDisc:add_elem_disc(neumannDisc)
domainDisc:add_elem_disc(neumannDisc2)
--domainDisc:add_post_process(dirichletBND)
--domainDisc:add_post_process(membraneDirichletBND)

-------------------------------------------
--  Setup Time Discretization
-------------------------------------------

-- create time discretization
timeDisc = ThetaTimeDiscretization()
timeDisc:set_domain_discretization(domainDisc)
timeDisc:set_theta(0.0) -- 0.0 is backward euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:set_dof_distribution(approxSpace:get_surface_dof_distribution())
op:init()

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()

-- exact Soler
exactSolver = LU()

-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(500)
	baseConvCheck:set_minimum_defect(1e-8)
	baseConvCheck:set_reduction(1e-30)
	baseConvCheck:set_verbose_level(false)
	base = LU()
	--base = LinearSolver()
	--base:set_convergence_check(baseConvCheck)
	--base:set_preconditioner(jac)
	
	-- Transfer and Projection
	transfer = utilCreateP1Prolongation(approxSpace)
	--transfer:set_dirichlet_post_process(dirichletBND)
	projection = utilCreateP1Projection(approxSpace)
	
	-- Gemoetric Multi Grid
	gmg = utilCreateGeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_surface_level(numRefs)
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

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(jac)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(gmg)
bicgstabSolver:set_convergence_check(convCheck)

-------------------------------------------
--  Non linear Solver
-------------------------------------------

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose_level(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(exactSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

newtonSolver:init(op)

-------------------------------------------
--  Solving
-------------------------------------------

-- get grid function
u = approxSpace:create_surface_function()

-- set initial value
u:set(0.0)

-- timestep in seconds
dt = 1
time = 0.0
step = 1
do_steps = 10

out = VTKOutput2d()
out:begin_timeseries("Con", u)
out:print("Con", u, 0, 0.0)

-- Apply Solver
PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, dt, out, "Con", true)

-- Output
out:end_timeseries("Con", u)