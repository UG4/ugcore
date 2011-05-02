----------------------------------------------------------
--
--   Lua - Script to perform a time dependent conv diff problem
--         with a membrane
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 2

if dim == 2 then
	gridName = "unit_square_tri_with_membrane.ugx"
end

numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numRefs    = util.GetParamNumber("-numRefs",    3)

-- choose number of time steps
NumPreTimeSteps = util.GetParamNumber("-numPreTimeSteps", 1)
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

--------------------------------
-- User Data Functions (begin)
--------------------------------

-- Parameters
eps = 1e-1
cx = 0.5
cy = 0.5

ax = 0.25
ay = 0.0

nu = 100
delta = 1e-1

dt = 1e-3

-- exact solution
function exactSolution(x, y, t)
	local xRot = math.cos(nu*t) * (x-cx) - math.sin(nu*t) * (y-cy) 
	local yRot = math.sin(nu*t) * (x-cx) + math.cos(nu*t) * (y-cy) 
	
	local expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t)
	local scale = delta/(delta+4*eps*t)

	return scale * math.exp(expo)
end

-- Functions depending on parameters
function ourDiffTensor2d(x, y, t)
	return	eps, 0, 
			0, eps
end
	
function ourVelocityField2d(x, y, t)
	return	nu*(y - cx), nu*(cy - x)
end
	

function ourDirichletBnd2d(x, y, t)
	return true, exactSolution(x, y, t)
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
	print("numPreRefs must be smaller than numRefs");
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
if sh:num_subsets() ~= 4 then 
print("Domain must have 4 Subsets for this problem.")
exit()
end

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct_on_subset("c_membrane", "Lagrange", 1, "Membrane, MembraneBnd")
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
--diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
diffusionMatrix = util.CreateConstDiagUserMatrix( eps, dim)
identityMatrix = util.CreateConstDiagUserMatrix( 1.0, dim)

-- Velocity Field setup
velocityField = util.CreateLuaUserVector("ourVelocityField"..dim.."d", dim)

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
constDirichlet = util.CreateConstBoundaryNumber(3.2, dim)

-- start value
startValue = util.CreateLuaUserNumber("exactSolution", dim)
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(1.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)

membraneElemDisc = util.CreateFV1ConvDiff(approxSpace, "c_membrane", "Membrane")
membraneElemDisc:set_upwind_amount(0.0)
membraneElemDisc:set_diffusion_tensor(identityMatrix)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(dirichlet, "c", "Boundary, MembraneBnd")

membraneDirichletBND = util.CreateDirichletBoundary(approxSpace)
membraneDirichletBND:add_boundary_value(constDirichlet, "c_membrane", "MembraneBnd")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
domainDisc:add_elem_disc(membraneElemDisc)
domainDisc:add_post_process(dirichletBND)
domainDisc:add_post_process(membraneDirichletBND)

-------------------------------------------
--  Setup Time Discretization
-------------------------------------------

-- create time discretization
timeDisc = ThetaTimeDiscretization()
timeDisc:set_domain_discretization(domainDisc)
timeDisc:set_theta(0.0) -- 0.0 is backward euler

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:set_dof_distribution(approxSpace:get_surface_dof_distribution())
op:init()

-- get grid function
u = approxSpace:create_surface_function()

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
convCheck:set_minimum_defect(1e-9)
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

-- create Exact solver
exactSolver = LU()

-- choose some solver
solver = bicgstabSolver

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose_level(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

newtonSolver:init(op)

-------------------------------------------
--  Apply Solver
-------------------------------------------

-- start
time = 0.0
step = 0

-- set initial value
print("Interpolation start values")
u:set(3.0)
InterpolateFunction(startValue, u, "c", time)

-- filename
filename = "Sol"

-- write start solution
print("Writing start values")
out = util.CreateVTKWriter(dim)
out:print(filename, u, step, time)

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   NumTimeSteps is " .. NumTimeSteps   .. ",  NumPreTimeSteps is " .. NumPreTimeSteps )

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

for step = 1, NumTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

	-- choose time step
	do_dt = dt
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, do_dt)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 

	-- update new time
	time = solTimeSeries:time(0) + do_dt
	
	-- plot solution
	out:print(filename, u, step, time)
	
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)

	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
out:write_time_pvd(filename, u)
