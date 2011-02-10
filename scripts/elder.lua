----------------------------------------------------------
--
--   Lua - Script to perform the Elder-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("../scripts/ug_util.lua")

-- choose algebra
algebra = CPUAlgebraChooser()
algebra:set_fixed_blocksize(1)
InitAlgebra(algebra)
-- InitAlgebra also loads all discretization functions and classes

-- constants
dim = 2
gridName = "elder_quads_8x2.ugx"
numPreRefs = 1
numRefs = 4
NumPreTimeSteps = 1
NumTimeSteps = 100

--------------------------------
-- User Data Functions (begin)
--------------------------------

function ConcentrationStart(x, y, t)
	if y == 150 then
		if x > 150 and x < 450 then
		return 1.0
		end
	end
	return 0.0
end

function PressureStart(x, y, t)
	return 9810 * (150 - y)
end

function ConcentrationDirichletBnd(x, y, t)
	if y == 150 then
		if x > 150 and x < 450 then
			return true, 1.0
		end
	end
	if y == 0.0 then
		return true, 0.0
	end

	return false, 0.0
end

function PressureDirichletBnd(x, y, t)
	if y == 150 then
		if x == 0.0 or x == 600 then
			return true, 9810 * (150 - y)
		end
	end
	
	return false, 0.0
end

function Porosity(x,y,t)
	return 0.1
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
if sh:num_subsets() ~= 2 then 
	print("Domain must have 2 Subsets for this problem.")
	exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("Boundary", 1)

-- write grid to file for test purpose
utilSaveDomain(dom, "refined_grid.ugx")


-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1Function(pattern, "c", 2)
AddP1Function(pattern, "p", 2)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-------------------------------------------
--  Setup User Functions
-------------------------------------------

-- dirichlet setup
ConcentrationDirichlet = utilCreateLuaBoundaryNumber("ConcentrationDirichletBnd", dim)
PressureDirichlet = utilCreateLuaBoundaryNumber("PressureDirichletBnd", dim)

-- start setup
ConcentrationStartValue = utilCreateLuaUserNumber("ConcentrationStart", dim)
PressureStartValue = utilCreateLuaUserNumber("PressureStart", dim)

-- Porosity
--porosityValue = utilCreateLuaUserNumber("Porosity", dim)
porosityValue = utilCreateConstUserNumber(0.1, dim)

-- Gravity
gravityValue = utilCreateConstUserVector(0.0, dim)
gravityValue:set_entry(dim-1, -9.81)

-- molecular Diffusion
molDiffusionValue = utilCreateConstDiagUserMatrix( 3.565e-6, dim)

-- Permeability
permeabilityValue = utilCreateConstDiagUserMatrix( 4.845e-13, dim)

-- Density
densityValue = NumberLinker2d();

-- Viscosity
viscosityValue = utilCreateConstUserNumber(1e-3, dim);

-----------------------------------------------------------------
--  Setup FV Element Discretization
-----------------------------------------------------------------

-- create Discretization
domainDisc = DomainDiscretization()

-- create dirichlet boundary for concentration
dirichletBND = utilCreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(ConcentrationDirichlet, "c", "Boundary")
dirichletBND:add_boundary_value(PressureDirichlet, "p", "Boundary")

-- create Finite-Volume Element Discretization for Convection Diffusion Equation
elemDisc = DensityDrivenFlow2d()
elemDisc:set_domain(dom)
elemDisc:set_pattern(pattern)
elemDisc:set_functions("c,p")
elemDisc:set_subsets("Inner")
if elemDisc:set_upwind("no") == false then exit() end
elemDisc:set_consistent_gravity(true)
elemDisc:set_boussinesq_transport(true)
elemDisc:set_boussinesq_flow(true)

elemDisc:set_porosity(porosityValue)
elemDisc:set_gravity(gravityValue)
elemDisc:set_permeability(permeabilityValue)
elemDisc:set_molecular_diffusion(molDiffusionValue)
elemDisc:set_density(densityValue)
elemDisc:set_viscosity(viscosityValue)

-- add Element Discretization to discretization
domainDisc:add_elem_disc(elemDisc)
domainDisc:add_post_process(dirichletBND)

-- create time discretization
timeDisc = ThetaTimeDiscretization()
timeDisc:set_domain_discretization(domainDisc)
timeDisc:set_theta(0.0) -- 0.0 is backward euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:set_dof_distribution(approxSpace:get_surface_dof_distribution())
op:init()

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

-- create GMG
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-10)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose_level(false)

-- base = LapackLUSolver()
base = BiCGStab()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

baseLU = LU()

-- Transfer and Projection
transfer = utilCreateP1Prolongation(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)
projection = utilCreateP1Projection(approxSpace)

gmg = utilCreateGeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_approximation_space(approxSpace)
gmg:set_surface_level(numRefs)
gmg:set_base_level(2)
gmg:set_base_solver(baseLU)
gmg:set_smoother(ilu)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_prolongation(transfer)
gmg:set_projection(projection)

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
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(0.5e-10)
convCheck:set_reduction(1e-8)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(gmg)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(gmg)
bicgstabSolver:set_convergence_check(convCheck)

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose_level(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

newtonSolver:init(op)

-- get grid function
u = approxSpace:create_surface_function()

-- timestep in seconds: 3153600 sec = 0.1 year
dt = 3.1536e6

time = 0.0
step = 1

-- set initial value
InterpolateFunction2d(PressureStartValue, u, "p", time)
InterpolateFunction2d(ConcentrationStartValue, u, "c", time)

-- Apply Solver
out = VTKOutput2d()
out:begin_timeseries("Elder", u)
out:print("Elder", u, 0, 0.0)

-- Perform Time Step
do_steps = NumPreTimeSteps
do_dt = dt/100
PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, do_dt, out, "Elder", true)
step = step + do_steps
time = time + do_dt * do_steps

do_steps = NumTimeSteps - NumPreTimeSteps
if do_steps > 0 then
	do_dt = dt
	PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, do_dt, out, "Elder", true)
	step = step + do_steps
	time = time + do_dt * do_steps
end

-- Output
out:end_timeseries("Elder", u)