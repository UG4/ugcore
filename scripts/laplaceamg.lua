----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem with AMG
--
--   Author: Martin Rupp, Andreas Vogel
--
----------------------------------------------------------

-- create Instance of a Domain
print("Create Domain.")
dom2d = Domain2d()

-- load domain
print("Load Domain from File.")
LoadDomain2d(dom2d, "../data/unit_square_tri.ugx")

-- create Refiner
print("Create Refiner")
refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom2d:get_grid())

REFINEMENTS = 5

for i= 1, REFINEMENTS do
	refiner:refine()
end


-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()

-- add function to pattern
print("Add one function to function pattern")
AddP1Function(pattern, "c", 2)

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

-- name subsets
_INNER_ = 0
_BND_ = 1

-- create Discretization
domainDisc = DomainDiscretization()

-- create dirichlet boundary
bndFct = SinusDirichletBoundaryFunction2d()
dirichletBND = DirichletBND2d()
dirichletBND:set_domain(dom2d)
dirichletBND:set_dirichlet_function(bndFct)
dirichletBND:set_function(_C_)

-- add Dirichlet boundary to discretization
domainDisc:add_dirichlet_bnd(dirichletBND, _BND_)

-- create Finite-Volume Element Discretization for Convection Diffusion Equation
sinusElemFct = SinusConvDiffUserFunction2d()
elemDisc = FV1ConvectionDiffusionElemDisc()
elemDisc:set_domain(dom2d)
elemDisc:set_upwind_amount(0.0)
elemDisc:set_user_functions(sinusElemFct)

-- add Element Discretization to discretization
fctGroup = FunctionGroup()
fctGroup:clear()
fctGroup:add_function( _C_ )
domainDisc:add(elemDisc, fctGroup, _INNER_, 2)

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:export_rhs(true)
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function("u", true)
b = approxSpace:create_surface_function("b", true)

-- set initial value
u:set(1.0)

-- init Operator
linOp:init()

-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")

-- set dirichlet values in start iterate
linOp:set_dirichlet_values(u)
b:assign(linOp:get_rhs())

-- create algebraic Preconditioner
jac = JacobiPreconditioner()
jac:set_damp(0.8)
gs = GSPreconditioner()
sgs = SGSPreconditioner()
bgs = BGSPreconditioner()
ilu = ILUPreconditioner()
ilut = ILUTPreconditioner()

-- create AMG
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-8)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose_level(false)

amg = AMGPreconditioner()


-- base = LapackLUSolver()
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(jac)

amg:set_nu1(2)
amg:set_nu2(2)
amg:set_gamma(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_debug(u)


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(30)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-12)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(amg)
linSolver:set_convergence_check(convCheck)

-- Apply Solver
ApplyLinearSolver(linOp, u, b, linSolver)

-- Output
WriteGridFunctionToVTK(u, "Solution")