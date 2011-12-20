----------------------------------------------------------
--
--   Example - Script to test a 2d reticulum setup
--
--   Author: Markus Breit
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- Dimension
dim = 2

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Grid
gridName = "simple_reticulum.ugx"

-- Refinements before distributing grid
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- Total refinements
numRefs    = util.GetParamNumber("-numRefs",    2)

-- choose number of time steps
NumTimeSteps =  util.GetParamNumber("-numTimeSteps", 5)

--------------------------------
-- User Data Functions (begin)
--------------------------------
    function CaCytStart(x, y, t)
        if 4.5 < x and x < 4.7 and 0.55 < y and y < 0.75
        then return 7.5e-4
        end
        return 7.5e-8
    end

    function CaERStart(x, y, t)
        return 1.0e-4
    end

    function IP3Start(x, y, t)
        return 1.0e-6
    end

	function ourDiffTensor2dCA(x, y, t)
		return	40, 0, 
				0, 40
	end
	
	function ourDiffTensor2dIP3(x, y, t)
		return	280, 0, 
				0, 280
	end
	
	function ourRhs2d(x, y, t)
		return 0;
	end
	
	function ourNeumannBnd2d(x, y, t)
		return true, 0.0
	end
	
	function ourDirichletBnd2d(x, y, t)
		return true, 0.0
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
dom = Domain()

-- load domain
print("Load Domain from File.")
if LoadDomain(dom, gridName) == false then
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
sh = dom:subset_handler()
if sh:num_subsets() ~= 4 then 
print("Domain must have 4 Subsets for this problem.")
exit()
end

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, "cyt, mem_er, mem_cyt")
approxSpace:add_fct("ca_er", "Lagrange", 1, "er, mem_er")
approxSpace:add_fct("ip3", "Lagrange", 1, "cyt, mem_er, mem_cyt")

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Start value function setup
    CaCytStartValue = LuaUserNumber("CaCytStart")
    CaERStartValue = LuaUserNumber("CaERStart")
    IP3StartValue = LuaUserNumber("IP3Start")

-- Diffusion Tensor setup
	diffusionMatrixCA = LuaUserMatrix("ourDiffTensor2dCA")
	diffusionMatrixIP3 = LuaUserMatrix("ourDiffTensor2dIP3")

-- rhs setup
	rhs = LuaUserNumber("ourRhs2d")
	--rhs = ConstUserNumber(0.0)

-- neumann setup
	neumann = LuaBoundaryNumber("ourNeumannBnd2d")
	--neumann = ConstUserNumber(0.0)

-- dirichlet setup
	dirichlet = LuaBoundaryNumber("ourDirichletBnd2d")
	--dirichlet = ConstBoundaryNumber(3.2)
	
-- dirichlet setup
	membraneDirichlet = LuaBoundaryNumber("membraneDirichletBnd2d")


-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them

 -- muss hier im letzten Arg (subsets) nicht auch ein mem_er stehen?
 -- Antwort: Nein, es geht hier im die 2d elemente von "er", der Rand von "er"
 --          spielt fuer das assemblieren keine Rolle. Randwerte kommen dann 
 --          spaeter extra. 

if dim == 2 then 
    upwind = NoUpwind2d()
elseif dim == 3 then 
    upwind = NoUpwind3d()
end

elemDiscER = ConvectionDiffusion("ca_er", "er") 
elemDiscER:set_disc_scheme("fv1")
elemDiscER:set_diffusion_tensor(diffusionMatrixCA)
elemDiscER:set_source(rhs)
if elemDiscER:set_upwind(upwind) == false then exit() end

elemDiscCYT = ConvectionDiffusion("ca_cyt", "cyt")
elemDiscCYT:set_disc_scheme("fv1")
elemDiscCYT:set_diffusion_tensor(diffusionMatrixCA)
elemDiscCYT:set_source(rhs)
if elemDiscCYT:set_upwind(upwind) == false then exit() end

elemDiscIP3 = ConvectionDiffusion("ip3", "cyt")
elemDiscIP3:set_disc_scheme("fv1")
elemDiscIP3:set_diffusion_tensor(diffusionMatrixIP3)
elemDiscIP3:set_source(rhs)
if elemDiscIP3:set_upwind(upwind) == false then exit() end

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDiscCYT = FV1NeumannBoundary("cyt")
neumannDiscCYT:add(neumann, "ca_cyt", "mem_cyt")

-- we pass here the function needed to evaluate the flux function.
-- The order in which the discrete fct are passed is crucial!
innerDisc = FV1InnerBoundaryCalciumER(3, "ca_cyt, ca_er, ip3", "mem_er")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

--dirichletBND = DirichletBoundary()
--dirichletBND:add(dirichlet, "c", "Boundary, MembraneBnd")

--membraneDirichletBND = DirichletBoundary()
--membraneDirichletBND:add(membraneDirichlet, "c_membrane", "MembraneBnd")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(neumannDiscCYT)
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
op:set_dof_distribution(approxSpace:surface_dof_distribution())
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
	transfer = P1Prolongation(approxSpace)
	--transfer:set_dirichlet_post_process(dirichletBND)
	projection = P1Projection(approxSpace)
	
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	--gmg:set_surface_level(numRefs)
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
newtonConvCheck:set_minimum_defect(5e-12)
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
u = GridFunction(approxSpace)

-- set initial value
InterpolateFunction(CaCytStartValue, u, "ca_cyt", 0.0)
InterpolateFunction(CaERStartValue, u, "ca_er", 0.0)
InterpolateFunction(IP3StartValue, u, "ip3", 0.0)

-- timestep in seconds
dt = 0.001
time = 0.0
step = 0

-- filename
filename = "Con"

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
