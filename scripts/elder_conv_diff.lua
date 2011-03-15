----------------------------------------------------------
--
--   Lua - Script to perform the Elder-Problem with a coupled conv diff eq
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

if dim == 2 then
gridName = "elder_quads_8x2.ugx"
else
	print("Only grid for 2D given.");
	exit();
end

-- choose number of pre-Refinements (before sending grid onto different processes)	
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 2)

-- choose number of time steps
NumPreTimeSteps = util.GetParamNumber("-numPreTimeSteps", 1)
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 100)

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

function TemperatureStart(x, y, t)
	if y > 30 and y < 60 then
		if x > 0 and x < 600 then
			return 1.0
		end
	end

	return 0.0
end

function TemperatureStartBak(x, y, t)
if y > 80 and y < 100 then
if x > 420 and x < 440 then
return 1.0
end
end

if y > 100 and y < 140 then
if x > 230 and x < 250 then
return 1.0
end
end

return 0.0
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
--refiner:refine()
util.GlobalRefineParallelDomain(dom)
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
SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct("p", "Lagrange", 1)
approxSpace:add_fct("T", "Lagrange", 1)
approxSpace:init()
approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------

-- dirichlet setup
ConcentrationDirichlet = util.CreateLuaBoundaryNumber("ConcentrationDirichletBnd", dim)
PressureDirichlet = util.CreateLuaBoundaryNumber("PressureDirichletBnd", dim)

-- start setup
ConcentrationStartValue = util.CreateLuaUserNumber("ConcentrationStart", dim)
PressureStartValue = util.CreateLuaUserNumber("PressureStart", dim)
TemperatureStartValue = util.CreateLuaUserNumber("TemperatureStart", dim)


-- Diffusion Tensor setup
--diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor2d", dim)
diffusionMatrix = util.CreateConstDiagUserMatrix(0.0, dim)

-- Velocity Field setup
function ourVelocityField2d(x, y, t)
	return	0, -0.000001
end

luaVelocityField = util.CreateLuaUserVector("ourVelocityField2d", dim)

constVelocityField = util.CreateConstUserVector(0.0, dim)
constVelocityField:set_entry(1, -0.000001)

-- Porosity
--porosityValue = util.CreateLuaUserNumber("Porosity", dim)
porosityValue = util.CreateConstUserNumber(0.1, dim)

-- Gravity
gravityValue = util.CreateConstUserVector(0.0, dim)
gravityValue:set_entry(dim-1, -9.81)

-- molecular Diffusion
molDiffusionValue = util.CreateConstDiagUserMatrix( 3.565e-6, dim)

-- Permeability
permeabilityValue = util.CreateConstDiagUserMatrix( 4.845e-13, dim)

-- Density
densityValue = ElderDensityLinker2d();

-- Viscosity
viscosityValue = util.CreateConstUserNumber(1e-3, dim);

-----------------------------------------------------------------
--  Setup FV Element Discretization
-----------------------------------------------------------------

-- create Discretization
domainDisc = DomainDiscretization()

-- create dirichlet boundary for concentration
dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(ConcentrationDirichlet, "c", "Boundary")
dirichletBND:add_boundary_value(PressureDirichlet, "p", "Boundary")

-- create Finite-Volume Element Discretization for Convection Diffusion Equation
elemDisc = DensityDrivenFlow2d()
elemDisc:set_approximation_space(approxSpace)
elemDisc:set_functions("c,p")
elemDisc:set_subsets("Inner")
if elemDisc:set_upwind("part") == false then exit() end
elemDisc:set_consistent_gravity(true)
elemDisc:set_boussinesq_transport(true)
elemDisc:set_boussinesq_flow(true)

densityValue:set_input(0, elemDisc:get_brine())

elemDisc:set_porosity(porosityValue)
elemDisc:set_gravity(gravityValue)
elemDisc:set_permeability(permeabilityValue)
elemDisc:set_molecular_diffusion(molDiffusionValue)
elemDisc:set_density(densityValue)
elemDisc:set_viscosity(viscosityValue)

darcyVelocityField = elemDisc:get_darcy_velocity()

-- create Finite-Volume Element Discretization for Convection Diffusion Equation
CDelemDisc = util.CreateFV1ConvDiff(approxSpace, "T", "Inner")
CDelemDisc:set_upwind_amount(1.0)
CDelemDisc:set_diffusion_tensor(diffusionMatrix)
--CDelemDisc:set_velocity_field(constVelocityField)
--CDelemDisc:set_velocity_field(luaVelocityField)
CDelemDisc:set_velocity_field(darcyVelocityField)
CDelemDisc:set_mass_scale(porosityValue)

-- add Element Discretization to discretization
domainDisc:add_elem_disc(elemDisc)
domainDisc:add_elem_disc(CDelemDisc)
domainDisc:add_post_process(dirichletBND)


-----------------------------------------------------------------
--  Setup Time Discretization
-----------------------------------------------------------------

-- create time discretization
timeDisc = ThetaTimeDiscretization()
timeDisc:set_domain_discretization(domainDisc)
timeDisc:set_theta(0.0) -- 0.0 is backward euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:set_dof_distribution(approxSpace:get_surface_dof_distribution())
op:init()

-----------------------------------------------------------------
--  Setup Solver
-----------------------------------------------------------------

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
transfer = util.CreateP1Prolongation(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)
projection = util.CreateP1Projection(approxSpace)

gmg = util.CreateGeometricMultiGrid(approxSpace)
gmg:set_discretization(timeDisc)
gmg:set_approximation_space(approxSpace)
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
step = 0

-- set initial value
InterpolateFunction(PressureStartValue, u, "p", time)
InterpolateFunction(ConcentrationStartValue, u, "c", time)
InterpolateFunction(TemperatureStartValue, u, "T", time)

----------------------------------
-- Time loop
----------------------------------

-- write start solution
print("Writing start values")
out = util.CreateVTKWriter(dim)
out:begin_timeseries("Elder", u)
out:print("Elder", u, step, time)

-- some info output
print( "   numRefs is         " .. numRefs)
print( "   NumTimeSteps is " .. NumTimeSteps   .. ",  NumPreTimeSteps is " .. NumPreTimeSteps )

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1, NumTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- choose time step
	if step <= NumPreTimeSteps then do_dt = dt/100
	else do_dt = dt
	end
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, do_dt)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- update new time
	time = solTimeSeries:time(0) + do_dt
	
	-- plot solution
	out:print("Elder", u, step, time)
	
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()
	
	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
out:end_timeseries("Elder", u)
