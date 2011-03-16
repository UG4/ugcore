----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Martin Rupp / Andreas Vogel
--
----------------------------------------------------------

ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraChooser());

-- constants
if util.HasParamOption("-3d") then
	dim = 3
else
	dim = 2
end

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_tri.ugx")
	-- gridName = "unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_cube_hex.ugx")
	-- gridName = "unit_cube_tets_regular.ugx"
end

-- choose number of pre-Refinements (before sending grid onto different processes)	
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 5)


--------------------------------
-- User Data Functions (begin)
--------------------------------
	function ourDiffTensor2d(x, y, t)
		return	1, 0, 
				0, 1
	end
	
	function ourVelocityField2d(x, y, t)
		return	0, 0
	end
	
	function ourReaction2d(x, y, t)
		return	0
	end
	
	function ourRhs2d(x, y, t)
		local s = 2*math.pi
		-- return	s*s*(math.sin(s*x) + math.sin(s*y))
		--return -2*y
		return 0
	end
	
	function ourNeumannBnd2d(x, y, t)
		--local s = 2*math.pi
		--return -s*math.cos(s*x)
		return true, -2*x*y
	end
	
	function ourDirichletBnd2d(x, y, t)
		local s = 2*math.pi
		return true, math.sin(s*x) + math.sin(s*y)
		--return true, x*x*y
		--return true, 2.5
	end

	function ourDiffTensor3d(x, y, z, t)
		return	1, 0, 0,
				0, 1, 0,
				0, 0, 1
	end
	
	function ourVelocityField3d(x, y, z, t)
		return	0, 0, 0
	end
	
	function ourReaction3d(x, y, z, t)
		return	0
	end
	
	function ourRhs3d(x, y, z, t)
		--local s = 2*math.pi
		--return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
		return 0;
	end
	
	function ourNeumannBnd3d(x, y, t)
		--local s = 2*math.pi
		--return -s*math.cos(s*x)
		return true, -2*x*y*z
	end
	
	function ourDirichletBnd3d(x, y, z, t)
		--local s = 2*math.pi
		--return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
		return true, x
	end
--------------------------------
-- User Data Functions (end)
--------------------------------


-- create Instance of a Domain
print("Create Domain.")
dom = util.CreateDomain(dim)

-- load domain
print("Load Domain from File.")
if util.LoadDomain(dom, gridName) == false then
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
-- which automatically creates a parallel refiner if required.
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
if sh:num_subsets() ~= 2 then 
print("Domain must have 2 Subsets for this problem.")
exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- write grid to file for test purpose
SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()
-- approxSpace:print_layout_statistic()
-- approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")

-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDiffTensor2d")
-- Diffusion Tensor setup
diffusionMatrix = util.CreateLuaUserMatrix("ourDiffTensor"..dim.."d", dim)
--diffusionMatrix = util.CreateConstDiagUserMatrix(1.0, dim)

-- Velocity Field setup
velocityField = util.CreateLuaUserVector("ourVelocityField"..dim.."d", dim)
--velocityField = util.CreateConstUserVector(0.0, dim)

-- Reaction setup
reaction = util.CreateLuaUserNumber("ourReaction"..dim.."d", dim)
--reaction = util.CreateConstUserNumber(0.0, dim)

-- rhs setup
rhs = util.CreateLuaUserNumber("ourRhs"..dim.."d", dim)
--rhs = util.CreateConstUserNumber(0.0, dim)

-- neumann setup
neumann = util.CreateLuaBoundaryNumber("ourNeumannBnd"..dim.."d", dim)
--neumann = util.CreateConstUserNumber(0.0, dim)

-- dirichlet setup
dirichlet = util.CreateLuaBoundaryNumber("ourDirichletBnd"..dim.."d", dim)
--dirichlet = util.CreateConstBoundaryNumber(3.2, dim)

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = util.CreateFV1ConvDiff(approxSpace, "c", "Inner")
elemDisc:set_upwind_amount(0.0)
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_rhs(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = util.CreateNeumannBoundary(approxSpace, "Inner")
--neumannDisc:add_boundary_value(neumann, "c", "NeumannBoundary")

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
--domainDisc:add_elem_disc(neumannDisc)
domainDisc:add_post_process(dirichletBND)


-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:export_rhs(true)
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = approxSpace:create_surface_function()
b = approxSpace:create_surface_function()

-- set initial value
u:set(1.0)

-- init Operator
print ("Assemble Operator ... ")
linOp:init()
print ("done")

-- set dirichlet values in start iterate
print ("set dirichlet... ")
linOp:set_dirichlet_values(u)
b:assign(linOp:get_rhs())
print ("done")
-- write matrix for test purpose
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
-- SaveVectorForConnectionViewer(b, "Rhs.mat")

-- create algebraic Preconditioners

print ("create preconditioners... ")
jac = Jacobi()
jac:set_damp(0.8)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()


-- create Base Solver
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-16)
baseConvCheck:set_reduction(1e-16)
baseConvCheck:set_verbose_level(false)

if false then
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
else
	base = LU()
end

-- Testvectors for AMG ---
--------------------------

function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
	local amgTestvector;
	if dim == 1 then
		amgTestvector = GridFunctionVectorWriter1d()
	elseif dim == 2 then
		amgTestvector = GridFunctionVectorWriter2d()
	elseif dim == 3 then
		amgTestvector = GridFunctionVectorWriter3d()
	end
	amgTestvector:set_reference_grid_function(gridfunction)
	amgTestvector:set_user_data(util.CreateLuaUserNumber(luaCallbackName, dim))
	return amgTestvector	
end


function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	local amgDirichlet0 = GridFunctionVectorWriterDirichlet02d()
	amgDirichlet0:init(dirichletBND, approxSpace)
	return amgDirichlet0
end


function ourTestvector2d_0_0(x, y, t)
	return 0
end

function ourTestvector2d_1_1(x, y, t)
	return math.sin(math.pi*x)*math.sin(math.pi*y)
end

function ourTestvector2d_2_1(x, y, t)
	return math.sin(2*math.pi*x)*math.sin(math.pi*y)
end


function ourTestvector2d_1_2(x, y, t)
	return math.sin(math.pi*x)*math.sin(2*math.pi*y)
end


function ourTestvector2d_2_2(x, y, t)
	return math.sin(2*math.pi*x)*math.sin(2*math.pi*y)
end



-- create AMG ---
-----------------
bUseFAMG = 1
if bUseFAMG == 1 then
	amg = FAMGPreconditioner()	
	amg:set_delta(0.5)
	amg:set_theta(0.95)
	amg:set_aggressive_coarsening(true)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.8)
	amg:set_testvector_damps(5)
	
	
	-- add testvectors with lua callbacks (see ourTestvector2d_1_1)
	-- amg:add_vector_writer(CreateAMGTestvector(u, "ourTestvector2d_0_0", dim), 1.0)
	-- amg:add_vector_writer(CreateAMGTestvector(u, "ourTestvector2d_1_1", dim), 1.0)
	-- amg:add_vector_writer(CreateAMGTestvector(u, "ourTestvector2d_1_2", dim), 1.0)
	amg:add_vector_writer(CreateAMGTestvector(u, "ourTestvector2d_2_1", dim), 1.0)
	-- amg:add_vector_writer(CreateAMGTestvector(u, "ourTestvector2d_2_2", dim), 1.0)
	
	-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
	-- amg:add_vector_writer(CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace), 1.0)

else
	amg = AMGPreconditioner()
	amg:enable_aggressive_coarsening_A(2)
end


vectorWriter = GridFunctionPositionProvider2d()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider2d(vectorWriter)
 
amg:set_num_presmooth(3)
amg:set_num_postsmooth(3)
amg:set_cycle_type(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
amg:set_max_levels(2)
amg:set_matrix_write_path("/Users/mrupp/matrices/")
amg:set_max_nodes_for_base(5)
amg:set_max_fill_before_base(0.7)
amg:set_fsmoothing(0.0)

amg:tostring()




-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(30)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

print("done.")
-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(amg)
linSolver:set_convergence_check(convCheck)

-- Apply Solver
ApplyLinearSolver(linOp, u, b, linSolver)

-- amg:check(u, b);

-- SaveVectorForConnectionViewer(u, "u.mat")

-- Output
WriteGridFunctionToVTK(u, "Solution")
