----------------------------------------------------------
--
--   Example - Script to test a 2d reticulum setup
--
--   Author: Markus Breit
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- Dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Grid
gridName = "rc19_amp.ugx"
--gridName = "RC19amp_ug4_finished.ugx"
--gridName = "simple_reticulum_3d.ugx"

-- Refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- Total refinements
numRefs    = util.GetParamNumber("-numRefs",    0)

-- choose number of time steps
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)

--------------------------------
-- User Data Functions (begin)
--------------------------------
function CaCytStart(x, y, z, t)
    return 7.5e-8
end

function CaERStart(x, y, z, t)
    return 1.0e-4
end

function IP3Start(x, y, z, t)
    return 1.0e-6
end

function ourDiffTensorCA(x, y, z, t)
    return	40, 0, 0,
            0, 40, 0,
            0, 0, 40
end

function ourDiffTensorIP3(x, y, z, t)
    return	280, 0, 0,
            0, 280, 0,
            0, 0, 280
end

function ourRhs(x, y, z, t)
    return 0;
end


-- firing pattern
syns = {}
for i=6,15 do
	syns["start"..i] = 0.005*(i-6)
	syns["end"..i] = 0.005*(i-6)+1
end

function ourNeumannBndCA(x, y, z, t, si)
	-- burst for active synapses
	if 	(si>=6 and si<=15 and syns["start"..si]<=t and t<syns["end"..si])
	then influx = 5e-6 * 11.0/16.0*(1.0+5.0/((10.0*(t-syns["start"..si])+1)*(10.0*(t-syns["start"..si])+1)))
	else influx = 0.0
	end
	
    return true, influx
end

ip3EntryDelay = 0.05;
function ourNeumannBndIP3(x, y, z, t, si)
	-- burst for active synapses
	if 	(si>=6 and si<=15 and syns["start"..si]+ip3EntryDelay<=t and t<syns["end"..si]+ip3EntryDelay)
	then influx = 0.0001
	else influx = 0.0
	end
	
    return true, influx
end

	
--------------------------------
-- User Data Functions (end)
--------------------------------

-- create, load, refine and distribute domain
print("create, refine and distribute domain")
neededSubsets = {}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create approximation space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
innerDomain = "er, mem_er"
--outerDomain = "cyt, mem_cyt, mem_er"
outerDomain = "cyt, nuc, mem_cyt, mem_er, mem_nuc"
synapses = ""
---[[
for i=1,9 do
	synapses = synapses .. ", syn" .. i
end
outerDomain = outerDomain .. synapses
--]]
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Start value function setup
    CaCytStartValue = LuaUserNumber3d("CaCytStart")
    CaERStartValue = LuaUserNumber3d("CaERStart")
    IP3StartValue = LuaUserNumber3d("IP3Start")

-- Diffusion Tensor setup
	diffusionMatrixCA = LuaUserMatrix3d("ourDiffTensorCA")
	diffusionMatrixIP3 = LuaUserMatrix3d("ourDiffTensorIP3")

-- rhs setup
	rhs = LuaUserNumber3d("ourRhs")
	--rhs = ConstUserNumber(0.0)

-- neumann setup
	neumannCA = LuaCondUserNumber3d("ourNeumannBndCA")
	--neumannCA = ConstUserNumber(0.0)
	neumannIP3 = LuaCondUserNumber3d("ourNeumannBndIP3")

--[[
-- dirichlet setup
	dirichlet = LuaCondUserNumber3d("ourDirichletBnd")
	
-- dirichlet setup
	membraneDirichlet = LuaCondUserNumber3d("membraneDirichletBnd")
--]]

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them

if dim == 2 then 
    upwind = NoUpwind2d()
elseif dim == 3 then 
    upwind = NoUpwind3d()
end

elemDiscER = ConvectionDiffusion("ca_er", "er") 
elemDiscER:set_disc_scheme("fv1")
elemDiscER:set_diffusion(diffusionMatrixCA)
elemDiscER:set_source(rhs)
elemDiscER:set_upwind(upwind)

elemDiscCYT = ConvectionDiffusion("ca_cyt", "cyt, nuc")
elemDiscCYT:set_disc_scheme("fv1")
elemDiscCYT:set_diffusion(diffusionMatrixCA)
elemDiscCYT:set_source(rhs)
elemDiscCYT:set_upwind(upwind)

elemDiscIP3 = ConvectionDiffusion("ip3", "cyt, nuc")
elemDiscIP3:set_disc_scheme("fv1")
elemDiscIP3:set_diffusion(diffusionMatrixIP3)
elemDiscIP3:set_source(rhs)
elemDiscIP3:set_upwind(upwind)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDiscCA = NeumannBoundary("cyt")
neumannDiscCA:add(neumannCA, "ca_cyt", "mem_cyt" .. synapses)
neumannDiscIP3 = NeumannBoundary("cyt")
neumannDiscIP3:add(neumannIP3, "ip3", "mem_cyt" .. synapses)

-----------------------------------------------------------------
--  Setup inner boundary (channels on ER membrane)
-----------------------------------------------------------------

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
innerDisc = FV1InnerBoundaryCalciumER("ca_cyt, ca_er, ip3", "mem_er")

--[[
-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary, MembraneBnd")

membraneDirichletBND = DirichletBoundary()
membraneDirichletBND:add(membraneDirichlet, "c_membrane", "MembraneBnd")

--]]

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(neumannDiscCA)
domainDisc:add(innerDisc)
--domainDisc:add(dirichletBND)
--domainDisc:add(membraneDirichletBND)

-------------------------------------------
--  Setup Time Discretization
-------------------------------------------

-- create time discretization
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
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

-- exact Solver
exactSolver = LU()


-- create GMG ---
-----------------

	-- Base Solver
	baseConvCheck = StandardConvergenceCheck()
	baseConvCheck:set_maximum_steps(10000)
	baseConvCheck:set_minimum_defect(1e-28)
	baseConvCheck:set_reduction(1e-2)
	baseConvCheck:set_verbose(false)
	--base = LU()
	---[[
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(gs)
	--]]
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	--gmg:set_surface_level(numRefs)
	gmg:set_base_level(0)
	gmg:set_base_solver(base)
	gmg:set_smoother(gs)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)

--[[
-- create AMG ---
-----------------

	amg = RSAMGPreconditioner()
	amg:set_nu1(2)
	amg:set_nu2(2)
	amg:set_gamma(1)
	amg:set_presmoother(gs)
	amg:set_postsmoother(gs)
	amg:set_base_solver(base)
	--amg:set_debug(u)
--]]

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(50)
convCheck:set_minimum_defect(1e-21)
convCheck:set_reduction(1e-4)
convCheck:set_verbose(true)

--[[
-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(jac)
cgSolver:set_convergence_check(convCheck)
--]]

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(gmg)
bicgstabSolver:set_convergence_check(convCheck)

-------------------------------------------
--  Non linear Solver
-------------------------------------------

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(20)
newtonConvCheck:set_minimum_defect(1e-21)
newtonConvCheck:set_reduction(1e-8)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)

newtonSolver:init(op)

-------------------------------------------
--  Solving
-------------------------------------------

-- get grid function
u = GridFunction(approxSpace)

-- set initial value
Interpolate(CaCytStartValue, u, "ca_cyt", 0.0)
Interpolate(CaERStartValue, u, "ca_er", 0.0)
Interpolate(IP3StartValue, u, "ip3", 0.0)

SaveVectorForConnectionViewer(u,"debug_vector.vec")

-- timestep in seconds
dt = 0.001
time = 0.0
step = 0

-- filename
filename = "/Users/markus/Developing/ug4/trunk/bin/retic/result"

-- write start solution
print("Writing start values")
out = VTKOutput()
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
	
	-- choose time step (currently constant)
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
