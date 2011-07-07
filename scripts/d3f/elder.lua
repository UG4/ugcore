----------------------------------------------------------
--
--   Lua - Script to perform the Elder-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

ug_load_script("../scripts/ug_util.lua")

-- choose algebra
algebra = CPUAlgebraSelector()
algebra:set_fixed_blocksize(1)
InitAlgebra(algebra)
-- InitAlgebra also loads all discretization functions and classes

-- choose dimension
dim = 2

-- choose grid
if dim == 2 then
	gridName = "grids/elder_quads_8x2.ugx"
else
	gridName = "grids/elder_hex_8x8x2.ugx"
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

if dim == 2 then 
--------------------------------
-- User Data Functions 2D
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

else 
--------------------------------
-- User Data Functions 3D
--------------------------------
	
function ConcentrationStart(x, y, z, t)
	if z == 150 then
		if y > 150 and y < 450 then
			if x > 150 and x < 450 then
				return 1.0
			end
		end
	end
	return 0.0
end

function PressureStart(x, y, z, t)
	return 9810 * (150 - z)
end

function ConcentrationDirichletBnd(x, y, z, t)
	if z == 150 then
		if y > 150 and y < 450 then
			if x > 150 and x < 450 then
				return true, 1.0
			end
		end
	end
	if z == 0.0 then
		return true, 0.0
	end
	
	return false, 0.0
end

function PressureDirichletBnd(x, y, z, t)
	if z == 150 then
		if y == 0.0 or y == 600 then
			if x == 0.0 or x == 600 then
				return true, 9810 * (150 - z)
			end
		end
	end
	
	return false, 0.0
end

function Porosity(x,y,z,t)
	return 0.1
end

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
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct("p", "Lagrange", 1)
approxSpace:init()
approxSpace:print_statistic()
approxSpace:print_layout_statistic()

-- lets order indices using Cuthill-McKee
if OrderCuthillMcKee(approxSpace, false) == false then
	print("ERROR when ordering Cuthill-McKee"); exit();
end

-------------------------------------------
--  Setup User Functions
-------------------------------------------

-- dirichlet setup
ConcentrationDirichlet = util.CreateLuaBoundaryNumber("ConcentrationDirichletBnd", dim)
PressureDirichlet = util.CreateLuaBoundaryNumber("PressureDirichletBnd", dim)
-- COMMENT IN BELOW TO USE C++ HARD CODED STANDARD ELDER BOUNDARY DATA
--ConcentrationDirichlet = ElderConcentrationBoundaryData2d()
--PressureDirichlet = ElderPressureBoundaryData2d()

-- start setup
ConcentrationStartValue = util.CreateLuaUserNumber("ConcentrationStart", dim)
PressureStartValue = util.CreateLuaUserNumber("PressureStart", dim)

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

-- Viscosity
viscosityValue = util.CreateConstUserNumber(1e-3, dim);

-- Density
function DensityFct(c)

	return 1000 + 200 * c
end

function DDensityFct_c(c)
	return 200
end


if dim == 2 then 
--	densityValue = ElderDensityLinker2d();
	densityValue = LuaUserFunctionNumber2d();
	densityValue:set_lua_value_callback("DensityFct", 1);
	densityValue:set_lua_deriv_callback(0, "DDensityFct_c");
else 
--	densityValue = ElderDensityLinker3d();
	densityValue = LuaUserFunctionNumber3d();
	densityValue:set_lua_value_callback("DensityFct", 1);
	densityValue:set_lua_deriv_callback(0, "DDensityFct_c");
end

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
if dim == 2 then elemDisc = DensityDrivenFlow2d()
else  elemDisc = DensityDrivenFlow3d() end

elemDisc:set_approximation_space(approxSpace)
elemDisc:set_functions("c,p")
elemDisc:set_subsets("Inner")
elemDisc:set_consistent_gravity(false)
elemDisc:set_boussinesq_transport(true)
elemDisc:set_boussinesq_flow(true)


-- Select upwind
upwind = util.CreateUpwind("part", dim)
if elemDisc:set_upwind(upwind) == false then exit() end


densityValue:set_input(0, elemDisc:get_brine())

print("Setting Porosity.")
elemDisc:set_porosity(porosityValue)
print("Setting Gravity.")
elemDisc:set_gravity(gravityValue)
print("Setting Permeability.")
elemDisc:set_permeability(permeabilityValue)
print("Setting mol. Diffusion.")
elemDisc:set_molecular_diffusion(molDiffusionValue)
print("Setting Density User func.")
elemDisc:set_density(densityValue)
print("Setting Viscosity.")
elemDisc:set_viscosity(viscosityValue)

-- add Element Discretization to discretization
print("Adding elem disc to global problem.")
domainDisc:add_elem_disc(elemDisc)
print("Adding bnd conds to global problem.")
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
--ilu:set_debug(dbgWriter)
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
baseJac = Jacobi();
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(baseJac)

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
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_prolongation(transfer)
gmg:set_projection(projection)

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

-- create Convergence Check
print("Creating Solver.")
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(0.5e-8)
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
bicgstabSolver:set_preconditioner(ilu)
bicgstabSolver:set_convergence_check(convCheck)

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-6)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose_level(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

-- timestep in seconds: 3153600 sec = 0.1 year
dt = 3.1536e6

time = 0.0
step = 0

-- set initial value
print("Interpolation start values")
InterpolateFunction(PressureStartValue, u, "p", time)
InterpolateFunction(ConcentrationStartValue, u, "c", time)

----------------------------------
-- Time loop
----------------------------------
-- filename
filename = "Elder"

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
