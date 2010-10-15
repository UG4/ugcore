----------------------------------------------------------
--
--   Lua - Script to perform the Elder-Problem
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom2d = Domain2d()

-- load domain
print("Load Domain from File.")
if LoadDomain2d(dom2d, "elder_quads_8x2.ugx") == false then 
   print("Loading Domain failed.")
   exit()
end

-- get subset handler
sh = dom2d:get_subset_handler()

-- create Refiner
print("Create Refiner")
refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom2d:get_grid())
refiner:refine()
refiner:refine()

-- write grid to file for test purpose
SaveDomain2d(dom2d, "refined_grid.ugx")

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(dom2d:get_subset_handler())

-- add function to pattern
print("Add one function to function pattern")
AddP1Function(pattern, "c", 2)
AddP1Function(pattern, "p", 2)

-- lock pattern
print("Lock function pattern")
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace2d()

-- assign function pattern
print("Assign Domain and FunctionPattern to ApproximationSpace")
approxSpace:assign_domain(dom2d)
approxSpace:assign_function_pattern(pattern)

-- name function
_C_ = 0
_P_ = 1

-- name subsets
if sh:num_subsets() ~= 2 then 
   print("Domain must have 2 Subsets for this problem.")
   exit()
end
_INNER_ = 0
_BND_ = 1
sh:set_subset_name("Inner", 0)
sh:set_subset_name("Boundary", 1)

-- create Discretization
domainDisc = DomainDiscretization()

-- create dirichlet boundary for concentration
cBndFct = cElderDirichletBoundaryFunction2d()
dirichletBND1 = DirichletBND2d()
dirichletBND1:set_domain(dom2d)
dirichletBND1:set_dirichlet_function(cBndFct)
dirichletBND1:set_function(_C_)

-- create dirichlet boundary for pressure
pBndFct = pElderDirichletBoundaryFunction2d()
dirichletBND2 = DirichletBND2d()
dirichletBND2:set_domain(dom2d)
dirichletBND2:set_dirichlet_function(pBndFct)
dirichletBND2:set_function(_P_)

-- add Dirichlet boundary to discretization
domainDisc:add_dirichlet_bnd(dirichletBND1, _BND_)
domainDisc:add_dirichlet_bnd(dirichletBND2, _BND_)

-- create Finite-Volume Element Discretization for Convection Diffusion Equation
elderElemFct = ElderUserFunction2d()
elemDisc = FV1DensityDrivenFlowElemDisc2d()
elemDisc:set_domain(dom2d)
elemDisc:set_upwind_amount(0.0)
elemDisc:set_user_functions(elderElemFct)

-- add Element Discretization to discretization
domainDisc:add(elemDisc, pattern, "c,p", "Inner")

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
jac = JacobiPreconditioner()
jac:set_damp(0.8)
gs = GSPreconditioner()
sgs = SGSPreconditioner()
bgs = BGSPreconditioner()
ilu = ILUPreconditioner()
ilut = ILUTPreconditioner()

-- create GMG
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-8)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose_level(false)

-- base = LapackLUSolver()
base = BiCGStabSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

gmg = GeometricMultiGridPreconditioner2d()
gmg:set_discretization(timeDisc)
gmg:set_approximation_space(approxSpace)
gmg:set_surface_level(2)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)

-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(1000)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheck)

-- create CG Solver
cgSolver = CGSolver()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheck)

-- create BiCGStab Solver
bicgstabSolver = BiCGStabSolver()
bicgstabSolver:set_preconditioner(ilu)
bicgstabSolver:set_convergence_check(convCheck)

-- convergence check
newtonConvCheck = StandardConvergenceCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(1e-6)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose_level(true)

newtonLineSearch = StandardLineSearch()

-- create Newton Solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)

newtonSolver:init(op)



-- get grid function
u = approxSpace:create_surface_function("u", true)

-- set initial value
interpol = InterpolateElder()
interpol:invoke(u)

-- Apply Solver
out = VTKOutput2d()
out:begin_timeseries("Elder", u)
out:print("Elder", u, 0, 0.0)

-- timestep in seconds: 3153600 sec = 0.1 year
dt = 3.1536e6

time = 0.0
step = 1

-- Perform Time Step
do_steps = 5
PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, dt/100, out, "Elder")
step = step + do_steps
time = time + dt/100 * do_steps

do_steps = 15
PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, dt/10, out, "Elder")
step = step + do_steps
time = time + dt/10 * do_steps

do_steps = 10
PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, dt/5, out, "Elder")
step = step + do_steps
time = time + dt/5 * do_steps

do_steps = 100
PerformTimeStep2d(newtonSolver, u, timeDisc, do_steps, step, time, dt, out, "Elder")
step = step + do_steps
time = time + dt * do_steps

-- Output
out:end_timeseries("Elder", u)