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
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 2
gridName = "unit_square_01/unit_square_01_tri_2x2.ugx"
numPreRefs = 0
numRefs = 4



-- create Instance of a Domain
print("Create Domain.")
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
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

if util.DistributeDomain(dom) == false then
print("Error while Distributing Grid.")
exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
GlobalRefineParallelDomain2d(dom)
end

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- Diffusion Tensor setup
diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
velocityField = util.CreateConstUserVector(0.0, dim)

-- Reaction setup
reaction = util.CreateConstUserNumber(0.0, dim)

-- rhs setup
rhs = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
dirichlet = util.CreateConstBoundaryNumber(0, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = util.CreateDirichletBoundary(approxSpace)
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
transfer = util.CreateP1Prolongation(approxSpace)
transfer:set_dirichlet_post_process(dirichletBND)
projection = util.CreateP1Projection(approxSpace)

-- Gemoetric Multi Grid
gmg = util.CreateGeometricMultiGrid(approxSpace)
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