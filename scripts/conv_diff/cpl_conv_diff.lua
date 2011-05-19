--------------------------------------------------------------------------------
--
--   Lua - Script for a system of coupled Conv-Diff-Eq
--
--	This Lua Script sets up a problem of several Convection-Diffusion-Equations.
--	It is intended to test the linker facility.
--
--	Let N be the number of unknown functions c1, ..., cN. For each unknown
--	the following conv-diff eq is discretized:
--
--	\partial_t c_i + \nabla ( \vec{q} ) = 0
--
--	where \vec{q} := 1/N \sum_{i=1}^N ( c_i * \vec{v} - D \nabla c_i ) 
--
--	\vec{v} and D are given user-defined velocity fields and Diffusion tensors.
--	The initial condition for each unknown solution is chosen equally, such that
--	the solved system must have the same solution in each component c_i. But
--	please note, that the resulting linear system is fully-coupled between the
--	components, i.e. with the FE pattern there is a non-zero entry for the 
--  whole dense block matrix of size NxN.
--
--   Author: Andreas Vogel
--
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 2

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	gridName = "unit_square_01/unit_square_01_quads_8x8.ugx"
end

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    4)
numSys     = util.GetParamNumber("-numSys",     1)

-- choose number of time steps
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 10)

print(" Choosen Parater:")
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("    grid       = " .. gridName)

--------------------------------
-- User Data Functions (begin)
--------------------------------

-- This scales the amount of diffusion of the problem
eps = 1e-1

-- The coordinates (cx, cy) specify the rotation center of the cone
cx = 0.5
cy = 0.5

-- The coordinates (ax, ay) specify the position of the highest point of the
-- cone at start time t=0.0
ax = 0.25
ay = 0.0

-- The parameter nu specifies the rotation velocity
nu = 100

-- The parameter delta is a scaling factor influencing the steepness of the cone
delta = 1e-1

-- This is the choosen time step
dt = 1e-3

-- This is the exact solution for our problem
function exactSolution(x, y, t)
	local xRot = math.cos(nu*t) * (x-cx) - math.sin(nu*t) * (y-cy) 
	local yRot = math.sin(nu*t) * (x-cx) + math.cos(nu*t) * (y-cy) 
	
	local expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t)
	local scale = delta/(delta+4*eps*t)

	return scale * math.exp(expo)
end

-- The Diffusion tensor
function ourDiffTensor2d(x, y, t)
	return	eps, 0, 
			0, eps
end
	
-- The velocity field
function ourVelocityField2d(x, y, t)
	return	nu*(y - cx) / numSys , nu*(cy - x) / numSys
end
	
-- The dirichlet condition
function ourDirichletBnd2d(x, y, t)
	return true, exactSolution(x, y, t)
end

--------------------------------
-- User Data Functions (end)
--------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
for i = 1, numSys do
	approxSpace:add_fct("c"..i, "Lagrange", 1)
end
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
--diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
diffusionMatrix = util.CreateConstDiagUserMatrix( eps, dim)

-- Velocity Field setup
velocityField = util.CreateLuaUserVector("ourVelocityField"..dim.."d", dim)

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)

-- start value
startValue = util.CreateLuaUserNumber("exactSolution", dim)
	
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------
Flux = ScaleAddLinkerVector2d()
scalarDiff = util.CreateConstUserNumber(-eps / numSys , dim)

elemDisc = {}
for i=1, numSys do
	elemDisc[i] = FV1ConstantEquation2d()
	elemDisc[i]:set_approximation_space(approxSpace)
	elemDisc[i]:set_functions("c"..i)
	elemDisc[i]:set_subsets("Inner")

	Flux:add(elemDisc[i]:get_concentration(), velocityField)
	Flux:add(scalarDiff, elemDisc[i]:get_concentration_grad())

	elemDisc[i]:set_mass_scale(elemDisc[i]:get_concentration())
	elemDisc[i]:set_velocity(Flux)
end
-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
for i=1, numSys do
	dirichletBND:add_boundary_value(dirichlet, "c"..i, "Boundary")
end

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
for i=1, numSys do
	domainDisc:add_elem_disc(elemDisc[i])
end
domainDisc:add_post_process(dirichletBND)

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
	base = LU()
--	base = LinearSolver()
--	base:set_convergence_check(baseConvCheck)
--	base:set_preconditioner(jac)
	
	-- Transfer and Projection
	transfer = util.CreateP1Prolongation(approxSpace)
	transfer:set_dirichlet_post_process(dirichletBND)
	projection = util.CreateP1Projection(approxSpace)
	
	-- Gemoetric Multi Grid
	gmg = util.CreateGeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(ilu)
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
solver = linSolver

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
--newtonSolver:set_debug(dbgWriter)
newtonSolver:init(op)

-------------------------------------------
--  Apply Solver
-------------------------------------------

-- start
time = 0.0
step = 0

-- set initial value
print("Interpolation start values")
for i=1, numSys do
	InterpolateFunction(startValue, u, "c"..i, time)
end

-- filename
filename = "Sol"

-- write start solution
print("Writing start values")
out = util.CreateVTKWriter(dim)
out:print(filename, u, step, time)

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   NumTimeSteps is " .. NumTimeSteps)

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

-- check if profiler is available
if GetProfilerAvailable() == true then
    -- get node
    pn = GetProfileNode("Ass_TJac")
    -- check if node is valid
    if pn:is_valid() then
        print("Called Ass_TJac "..pn:get_avg_entry_count().." times.")
        print("Spend "..pn:get_avg_self_time_ms().." ms for Ass_TJac alone (without child nodes).")
        print("Spend "..pn:get_avg_total_time_ms().." ms in Ass_TJac altogether.\n")
        
        print("Ass_TJac call list total time sorted")
        print(pn:total_time_sorted())
    else
        print("Ass_TJac is not known to the profiler.")
    end
else
    print("Profiler not available.")
end 

-- check if profiler is available
if GetProfilerAvailable() == true then
    -- get node
    pn = GetProfileNode("Ass_TDef")
    -- check if node is valid
    if pn:is_valid() then
        print("Called Ass_TDef "..pn:get_avg_entry_count().." times.")
        print("Spend "..pn:get_avg_self_time_ms().." ms for Ass_TDef alone (without child nodes).")
        print("Spend "..pn:get_avg_total_time_ms().." ms in Ass_TDef altogether.\n")
        
        print("Ass_TDef call list total time sorted")
        print(pn:total_time_sorted())
    else
        print("Ass_TDef is not known to the profiler.")
    end
else
    print("Profiler not available.")
end 
