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
numPreRefs = GetParamNumber("-numPreRefs", 0)

-- Total refinements
numRefs    = GetParamNumber("-numRefs",    2)

-- choose number of time steps
NumTimeSteps =  GetParamNumber("-numTimeSteps", 5)

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
		return	500, 0, 
				0, 500
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
dom = utilCreateDomain(dim)

-- load domain
print("Load Domain from File.")
if utilLoadDomain(dom, gridName) == false then
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
sh = dom:get_subset_handler()
if sh:num_subsets() ~= 4 then 
print("Domain must have 4 Subsets for this problem.")
exit()
end

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1FunctionOnSubsets(pattern, "ca_cyt", "cyt, mem_er, mem_cyt", dim)
AddP1FunctionOnSubsets(pattern, "ca_er", "er, mem_er", dim)
AddP1FunctionOnSubsets(pattern, "ip3", "cyt, mem_er, mem_cyt", dim)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Start value function setup
    CaCytStartValue = utilCreateLuaUserNumber("CaCytStart", dim)
    CaERStartValue = utilCreateLuaUserNumber("CaERStart", dim)
    IP3StartValue = utilCreateLuaUserNumber("IP3Start", dim)

-- Diffusion Tensor setup
	diffusionMatrixCA = utilCreateLuaUserMatrix("ourDiffTensor2dCA", dim)
	diffusionMatrixIP3 = utilCreateLuaUserMatrix("ourDiffTensor2dIP3", dim)

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

-- Note: No VelocityField and Reaction is set. The assembling assumes default
--       zero values for them

 -- muss hier im letzten Arg (subsets) nicht auch ein mem_er stehen?
 -- Antwort: Nein, es geht hier im die 2d elemente von "er", der Rand von "er"
 --          spielt für das assemblieren keine Rolle. Randwerte kommen dann 
 --          später extra. 
elemDiscER = utilCreateFV1ConvDiff(approxSpace, "ca_er", "er") 
elemDiscER:set_upwind_amount(0.0)
elemDiscER:set_diffusion_tensor(diffusionMatrixCA)
elemDiscER:set_rhs(rhs)

elemDiscCYT = utilCreateFV1ConvDiff(approxSpace, "ca_cyt", "cyt")
elemDiscCYT:set_upwind_amount(0.0)
elemDiscCYT:set_diffusion_tensor(diffusionMatrixCA)
elemDiscCYT:set_rhs(rhs)

elemDiscIP3 = utilCreateFV1ConvDiff(approxSpace, "ip3", "cyt")
elemDiscIP3:set_upwind_amount(0.0)
elemDiscIP3:set_diffusion_tensor(diffusionMatrixIP3)
elemDiscIP3:set_rhs(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

neumannDiscCYT = utilCreateNeumannBoundary(approxSpace, "cyt")
neumannDiscCYT:add_boundary_value(neumann, "ca_cyt", "mem_cyt")

-- we pass here the function needed to evaluate the flux function. The order in 
-- which the discrete fct are passed is crutial!
innerDisc = utilCreateInnerBoundary(approxSpace, "ca_cyt, ca_er, ip3", "mem_er")
--innerDisc:add_flux("ca_cyt, ip3", "mem_er")

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
domainDisc:add_elem_disc(elemDiscIP3)
domainDisc:add_elem_disc(neumannDiscCYT)
domainDisc:add_elem_disc(innerDisc)
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
u = approxSpace:create_surface_function()

-- set initial value
InterpolateFunction2d(CaCytStartValue, u, "ca_cyt", 0.0)
InterpolateFunction2d(CaERStartValue, u, "ca_er", 0.0)
InterpolateFunction2d(IP3StartValue, u, "ip3", 0.0)

-- timestep in seconds
dt = 0.01
time = 0.0
step = 1

-- write start solution
print("Writing start values")
out = utilCreateVTKWriter(dim)
out:begin_timeseries("Con", u)
out:print("Con", u, 0, time)

-- some info output
print( "   numPreRefs is   " .. numPreRefs ..     ",  numRefs is         " .. numRefs)
print( "   NumTimeSteps is " .. NumTimeSteps)

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
prevSol = PreviousSolutions()
prevSol:push(uOld, time)

for step = 1, NumTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	-- choose time step (currently constant)
	do_dt = dt
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(prevSol, do_dt)
	
	-- prepare newton solver
	if newtonSolver:prepare(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- apply newton solver
	if newtonSolver:apply(u) == false then print ("Newton solver failed at step "..step.."."); exit(); end 
	
	-- update new time
	time = prevSol:time(0) + do_dt
	
	-- plot solution
	out:print("Con", u, step, time)
	
	-- get oldest solution
	oldestSol = prevSol:oldest_solution()
	
	-- copy values into oldest solution (we reuse the memory here)
	VecScaleAssign(oldestSol, 1.0, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	prevSol:push_discard_oldest(oldestSol, time)
	
	print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- end timeseries, produce gathering file
out:end_timeseries("Con", u)