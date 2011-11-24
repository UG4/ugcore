-------------------------------------------------------------------------------
--
--   Lua - Script to test the time stepping
--
--	 This script uses the ConvectionDiffusion element disc but with zero
--	 convection, diffusion and reaction. Thus, the computed equation is
--		\partial_t u = f(t)
--
--   Author: Andreas Vogel
--
-------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

-- constants
dim = 2

gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")

numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numRefs    = util.GetParamNumber("-numRefs",    2)

-- choose number of time steps
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 10)
dt =  util.GetParamNumber("-dt", 1)
timeScheme = util.GetParam("-scheme", "Alexander")
maxOrder = util.GetParamNumber("-maxOrder", 2)

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-----------------------------------------------------------------
--  Discretization Setup
-----------------------------------------------------------------

-- This is the exact solution for our problem
function exactSolution(x, y, t)
	return math.exp(-.5*t)
end

-- The dirichlet condition
function dirichletBnd2d(x, y, t)
	return true, exactSolution(x, y, t)
end

-- the rhs
function rhs(x, y, t)
	return 1.0
end

-- the source
function reaction(x, y, t)
	return 0.5
end

dirichlet = LuaBoundaryNumber("dirichletBnd"..dim.."d")
startValue = LuaUserNumber("exactSolution")
sourceValue = LuaUserNumber("rhs")
reactionValue = LuaUserNumber("reaction")
	
elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_upwind(NoUpwind())
elemDisc:set_disc_scheme("fv1")
elemDisc:set_reaction(reactionValue)
--elemDisc:set_source(sourceValue)

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)

if timeScheme == "Theta" then 
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler

elseif timeScheme == "Alexander" then 
timeDisc = ThetaTimeStep(domainDisc, "Alexander")

elseif timeScheme == "FracStep" then 
timeDisc = ThetaTimeStep(domainDisc, "FracStep")

elseif timeScheme == "bdf" then
timeDisc = BDF(domainDisc)
timeDisc:set_order(1)

else
write("No time stepping specified"); exit();
end  

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:set_dof_distribution(approxSpace:surface_dof_distribution())
op:init()

-- get grid function
u = GridFunction(approxSpace)

ilu = ILU()

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-9)
convCheck:set_reduction(1e-12)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- choose some solver
solver = cgSolver

solver = LU()

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10000)
newtonConvCheck:set_verbose_level(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:init(op)

-------------------------------------------------------------------------------
--  Apply Solver
-------------------------------------------------------------------------------

-- start
time = 0.0
step = 0

-- set initial value
print("Interpolation start values")
InterpolateFunction(startValue, u, "c", time)

-- filename
filename = "Sol"

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

for step = 1, NumTimeSteps + maxOrder - 1 do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")

	-- choose time step
	if step <= maxOrder then
		do_dt = dt/maxOrder
	else
		do_dt = dt
	end
	
	for stage = 1, timeDisc:num_stages() do
		print("      +++ STAGE " .. stage .. " BEGIN ++++++")

		timeDisc:set_stage(stage)

		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
		
		-- prepare newton solver
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 
		
		-- apply newton solver
		if newtonSolver:apply(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 
	
		-- update new time
		time = timeDisc:future_time()
		
		-- push oldest solutions with new values to front, oldest sol pointer is poped from end	
		if timeScheme == "bdf" and step < maxOrder then
			print("Increasing order to "..step+1)
			timeDisc:set_order(step+1)
			uNew = u:clone()
			solTimeSeries:push(uNew, time)
		else 
			oldestSol = solTimeSeries:oldest()
			VecScaleAssign(oldestSol, 1.0, u)
			solTimeSeries:push_discard_oldest(oldestSol, time)
		end

		print("      +++ STAGE " .. stage .. " END   ++++++")
	end
	
	-- plot solution
	out:print(filename, u, step, time)

	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
out:write_time_pvd(filename, u)
