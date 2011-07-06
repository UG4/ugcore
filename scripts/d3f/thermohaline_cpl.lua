--------------------------------------------------------------------------------
--
--   Lua - Script to perform the parcel - problem
--
--   Author: Andreas Vogel
--
--------------------------------------------------------------------------------

ug_load_script("../scripts/ug_util.lua")

-- choose algebra
algebra = CPUAlgebraSelector()
algebra:set_fixed_blocksize(1)
InitAlgebra(algebra)

-- choose dimension
dim = 2

-- choose grid
if dim == 2 then
	gridName = "grids/parcel_square.ugx"
end

-- choose number of pre-Refinements (before sending grid onto different processes)	
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 2)

-- choose number of time steps
NumPreTimeSteps = util.GetParamNumber("-numPreTimeSteps", 1)
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 100)


--------------------------------------------------------------------------------
-- Grid and Approximation Space
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary", "Top", "Bottom"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct("p", "Lagrange", 1)
approxSpace:add_fct("T", "Lagrange", 1)
approxSpace:init()
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee
if OrderCuthillMcKee(approxSpace, false) == false then
	print("ERROR when ordering Cuthill-McKee"); exit();
end

--------------------------------------------------------------------------------
--  Constants
--------------------------------------------------------------------------------

local phi = 0.1
local permeability = 5e-12
local molDiffusion = 1e-8
local thermalConductivity = 1.8

local C_s = 1e3
local C_f = 4.184e3

local rho_ref = 1e3
local rho_s = 2.650e3

local ref_rho_pw = 800
local ref_rho_pb = 1000

local alpha_w = 1.45e-3
local alpha_b = 8.13e-4

local theta_w = 523.15
local theta_b = 563.15

local mass_scale = phi * rho_ref * C_f + (1-phi)* rho_s * C_s

wQ = 0.0
tQ = 473.15

dir = "posBouyancy"
if dir == "posBouyancy" then
wP = 0.47
tP = 573.15
else
wP = 0.55
tP = 523.15
end

-- Size of Domain [m]
local L = 2500
local l = 500

--------------------------------------------------------------------------------
--  User Data for Linking Process
--------------------------------------------------------------------------------

-- Porosity
Porosity = util.CreateConstUserNumber(phi, dim)

-- Diffusion
--MolecularDiffusion = util.CreateConstDiagUserMatrix(molDiffusion, dim)
DiffDisp = util.CreateConstDiagUserMatrix(molDiffusion, dim)

-- Mass Scale
TempMassScale = util.CreateConstUserNumber(mass_scale, dim);

--------------------------------------------------------------------------------
--  Density
--------------------------------------------------------------------------------

function DensityFct(c, T)
	local rho_pw = ref_rho_pw * math.exp(-alpha_w * (T-theta_w))
	local rho_pb = ref_rho_pb * math.exp(-alpha_b * (T-theta_b))

	return rho_pw * rho_pb / (rho_pb - (rho_pb - rho_pw) * c)
end

function DDensityFct_c(c, T)
	local rho_pw = ref_rho_pw * math.exp(-alpha_w * (T-theta_w))
	local rho_pb = ref_rho_pb * math.exp(-alpha_b * (T-theta_b))

	local denom = 1 / (rho_pb - (rho_pb - rho_pw) * c)

	return (rho_pb - rho_pw) * rho_pw * rho_pb * (denom*denom)
end

function DDensityFct_T(c, T)
	local rho_pw = ref_rho_pw * math.exp(-alpha_w * (T-theta_w))
	local rho_pb = ref_rho_pb * math.exp(-alpha_b * (T-theta_b))

	local D_rho_pw = - alpha_w * rho_pw
	local D_rho_pb = - alpha_b * rho_pb
	
	local denom = 1 / (rho_pb - (rho_pb - rho_pw) * c)

	return ((D_rho_pw * rho_pb + rho_pw * D_rho_pb) * denom
			- (rho_pw*rho_pb) * (denom*denom) * D_rho_pb
			+ (rho_pw*rho_pb) * (denom*denom) * D_rho_pb * c
			- (rho_pw*rho_pb) * (denom*denom) * D_rho_pw * c)
end

if 		dim == 2 then Density = LuaUserFunctionNumber2d();
elseif 	dim == 3 then Density = LuaUserFunctionNumber3d();
else print("Dimension not supported"); exit() end
	
Density:set_lua_value_callback("DensityFct", 2);
Density:set_lua_deriv_callback(0, "DDensityFct_c");
Density:set_lua_deriv_callback(1, "DDensityFct_T");

--Density = ElderDensityLinker2d()

--------------------------------------------------------------------------------
--  Darcy Velocity
--------------------------------------------------------------------------------

-- Permeability
Permeability = util.CreateConstDiagUserMatrix(permeability, dim)

-- Viscosity
Viscosity = util.CreateConstUserNumber(1e-3, dim);

-- Gravity
Gravity = util.CreateConstUserVector(0.0, dim); Gravity:set_entry(dim-1, -9.81)

-- Darcy Velocity
DarcyVelocity = util.CreateDarcyVelocityLinker(dim)
DarcyVelocity:set_permeability(Permeability)
DarcyVelocity:set_viscosity(Viscosity)
DarcyVelocity:set_density(Density)
DarcyVelocity:set_gravity(Gravity)

-- set the product Density * Porosity (only Porosity for Bussinesq)
rhophi = Porosity
--rhophi = ScaleAddLinkerNumber2d()
--rhophi:add(Porosity, RefDensity)

--------------------------------------------------------------------------------
--  Coefficients for Heat Equation
--------------------------------------------------------------------------------

-- Diffusion
ThermoDisp = util.CreateConstDiagUserMatrix(thermalConductivity, dim)

C_fValue = util.CreateConstUserNumber(C_f * rho_ref, dim)
C_f_DarcyVel = ScaleAddLinkerVector2d()
C_f_DarcyVel:add(C_fValue, DarcyVelocity)

--------------------------------------------------------------------------------
--  Setup FV Element Discretization
--------------------------------------------------------------------------------

FlowEq = util.CreateFV1ConstEq(approxSpace, "p", "Inner")
FlowEq:set_mass_scale(rhophi)
FlowEq:set_velocity(DarcyVelocity)
print("Flow Equation created.")

TransportEq = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
upwind = util.CreateUpwind("part", dim)
TransportEq:set_upwind(upwind)
TransportEq:set_mass_scale(rhophi)
TransportEq:set_velocity_field(DarcyVelocity)
TransportEq:set_diffusion_tensor(DiffDisp)
print("Transport Equation created.")

EnergyEq = util.CreateFV1ConvDiff(approxSpace, "T", "Inner")
upwindTemp = util.CreateUpwind("part", dim)
EnergyEq:set_upwind(upwindTemp)
EnergyEq:set_mass_scale(TempMassScale)
EnergyEq:set_velocity_field(C_f_DarcyVel)
EnergyEq:set_diffusion_tensor(ThermoDisp)
print("Transport Equation created.")


Density:set_input(0, TransportEq:get_concentration())
Density:set_input(1, EnergyEq:get_concentration())
DarcyVelocity:set_pressure_gradient(FlowEq:get_concentration_grad())

--------------------------------------------------------------------------------
--  Boundary Conditions
--------------------------------------------------------------------------------

function PressureDirichletBnd(x, y, t)
	if y == L then
		if x == 0.0 or x == L then
			return true, 0.0
		end
	end
	
	return false, 0.0
end

-- dirichlet setup
PressureDirichlet = util.CreateLuaBoundaryNumber("PressureDirichletBnd", dim)
ZeroDirichlet = util.CreateConstBoundaryNumber(0.0, dim)
TempDirichlet = util.CreateConstBoundaryNumber(tQ, dim)

-- create dirichlet boundary for concentration
dirichletBND = util.CreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(ZeroDirichlet, "c", "Top, Bottom")
dirichletBND:add_boundary_value(TempDirichlet, "T", "Top, Bottom")
dirichletBND:add_boundary_value(PressureDirichlet, "p", "Top")

--------------------------------------------------------------------------------
--  Spacial Problem
--------------------------------------------------------------------------------

-- create Discretization
domainDisc = DomainDiscretization()

-- add Element Discretization to discretization
domainDisc:add_elem_disc(TransportEq)
domainDisc:add_elem_disc(FlowEq)
domainDisc:add_elem_disc(EnergyEq)
domainDisc:add_post_process(dirichletBND)

--------------------------------------------------------------------------------
--  Time Discretization
--------------------------------------------------------------------------------

-- create time discretization
timeDisc = ThetaTimeDiscretization()
timeDisc:set_domain_discretization(domainDisc)
timeDisc:set_theta(0.0) -- 0.0 is backward euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:set_dof_distribution(approxSpace:get_surface_dof_distribution())
op:init()

--------------------------------------------------------------------------------
--  Solver Setup
--------------------------------------------------------------------------------

-- get grid function
u = approxSpace:create_surface_function()

-- debug writer
dbgWriter = util.CreateGridFunctionDebugWriter(dim)
dbgWriter:set_reference_grid_function(u)
dbgWriter:set_vtk_output(false)

-- create algebraic Preconditioner
print("Creating Preconditioner.")
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilu:set_debug(dbgWriter)
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
gmg:set_base_level(0)
gmg:set_base_solver(baseLU)
gmg:set_smoother(ilu)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_prolongation(transfer)
gmg:set_projection(projection)

amg = RSAMGPreconditioner()
--amg:enable_aggressive_coarsening_A(2)
amg:set_num_presmooth(3)
amg:set_num_postsmooth(3)
amg:set_cycle_type(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
amg:set_max_levels(5)
amg:set_max_nodes_for_base(300)
amg:set_max_fill_before_base(0.7)
amg:set_fsmoothing(true)

-- create Convergence Check
print("Creating Solver.")
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(0.5e-10)
convCheck:set_reduction(1e-20)

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
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

--------------------------------------------------------------------------------
-- Interpolate Start Values
--------------------------------------------------------------------------------

local lowerParcel = (L-l)/2
local upperParcel = L - ((L-l)/2)

function ConcentrationStart(x, y, t)
	if y > lowerParcel and y < upperParcel then
		if x > lowerParcel and x < upperParcel then
			return wP
		end
	end
	return wQ
end

function PressureStart(x, y, t)
	return 0.0
end

function TemperatureStart(x, y, t)
	if y > lowerParcel and y < upperParcel then
		if x > lowerParcel and x < upperParcel then
			return tP
		end
	end
	return tQ
end

ConcentrationStartValue = util.CreateLuaUserNumber("ConcentrationStart", dim)
PressureStartValue = util.CreateLuaUserNumber("PressureStart", dim)
TemperatureStartValue = util.CreateLuaUserNumber("TemperatureStart", dim)

time = 0.0
InterpolateFunction(PressureStartValue, u, "p", time)
InterpolateFunction(ConcentrationStartValue, u, "c", time)
InterpolateFunction(TemperatureStartValue, u, "T", time)

--------------------------------------------------------------------------------
-- Time loop
--------------------------------------------------------------------------------

-- timestep in seconds: 3.1536e6 sec = 0.1 year
dt = 3.1536e6

step = 0

-- filename
filename = "Elder"

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   NumTimeSteps is " .. NumTimeSteps   .. ",  NumPreTimeSteps is " .. NumPreTimeSteps )

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

-- write start solution
print("Writing start values")
out = util.CreateVTKWriter(dim)
out:print(filename, u, step, time)

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
