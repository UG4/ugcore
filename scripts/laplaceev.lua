----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-EV-Problem
--	 
--   Author: Andreas Vogel / Martin Rupp
--
----------------------------------------------------------

-- make sure that ug_util is in the right path.
-- currently only the path in which you start your application is valid.
dofile("../scripts/ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraChooser());

-- constants
dim = 2
gridName = "unit_square_tri.ugx"
numPreRefs = 0
numRefs = 4



-- create Instance of a Domain
print("Create Domain.")
dom = utilCreateDomain(dim)

-- load domain
print("Load Domain from File.")
if utilLoadDomain(dom, gridName) == false then
print("Loading Domain failed.")
exit()
end

-- get subset handler
sh = dom:get_subset_handler()
--if sh:num_subsets() ~= 2 then 
--	print("Domain must have 2 Subsets for this problem.")
--	exit()
--end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- create Refiner
print("Create Refiner")
if numPreRefs > numRefs then
print("numPreRefs must be smaller/equal than numRefs");
exit();
end

refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom:get_grid())
for i=1,numPreRefs do
refiner:refine()
end

if utilDistributeDomain(dom) == false then
print("Error while Distributing Grid.")
exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
GlobalRefineParallelDomain2d(dom)
end

-- write grid to file for test purpose
utilSaveDomain(dom, "refined_grid.ugx")

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1Function(pattern, "c", dim)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Diffusion Tensor setup
diffusionMatrix = utilCreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
velocityField = utilCreateConstUserVector(0.0, dim)

-- Reaction setup
reaction = utilCreateConstUserNumber(0.0, dim)

-- rhs setup
rhs = utilCreateConstUserNumber(0.0, dim)

-- dirichlet setup
dirichlet = utilCreateConstBoundaryNumber(0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = utilCreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_rhs(rhs)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = utilCreateDirichletBoundary(approxSpace)
dirichletBND:add_boundary_value(dirichlet, "c", "DirichletBoundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization()
domainDisc:add_elem_disc(elemDisc)
domainDisc:add_post_process(dirichletBND)

-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:export_rhs(false)
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- init Operator
print ("Assemble Operator ... ")
linOp:init()
print ("done")


-- get grid function
n_ev = 2
x0 = approxSpace:create_surface_function("x0", true)
x0:set_random(0.0, 1.0);
linOp:set_dirichlet_values(x0)
x1 = approxSpace:create_surface_function("x1", true)
x1:set_random(0.0, 1.0);
linOp:set_dirichlet_values(x1)

-- write matrix for test purpose
SaveMatrixForConnectionViewer(x0, linOp, "Stiffness.mat")
-- SaveVectorForConnectionViewer(b, "Rhs.mat")

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
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
transfer = utilCreateP1Prolongation(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)
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


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

eig = EigenSolver()
eig:add_vector(x0)
eig:add_vector(x1)
eig:set_preconditioner(gmg)
eig:set_linear_operator_A(linOp);
eig:apply()


-- Output
WriteGridFunctionToVTK(u, "Solution")