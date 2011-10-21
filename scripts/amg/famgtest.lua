----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Martin Rupp / Andreas Vogel
--
----------------------------------------------------------
-- SetOutputProcessRank(0)
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")

-- constants
if util.HasParamOption("-3d") then
	dim = 3
else
	dim = 2
end

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	-- gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	-- gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 5)

-- choose number of pre-Refinements (before sending grid onto different processes)	
numPreRefs = util.GetParamNumber("-numPreRefs", math.min(5, numRefs-2))

maxBase = util.GetParamNumber("-maxBase", 3000)

RAepsilon = util.GetParamNumber("-RAepsilon", 1)
RAalpha = util.GetParamNumber("-RAalpha", 0)

bFileOutput = true
bOutput = true
bUseNestedAMG = true
print("Parameters: ")
print("    numPreRefs = "..numPreRefs)
print("    numRefs = "..numRefs)
print("    maxBase = "..maxBase)
print("    dim = "..dim)
print("    gridName = "..gridName)
print("    RAepsilon = "..RAepsilon)
print("    RAalpha = "..RAalpha.." degree")
RAalpha = RAalpha * (2*math.pi/360)
print("    RAalpha = "..RAalpha.." grad")


function my_assert(condition, text)
	if condition == false then
		error(text)
	end
end

--------------------------------
-- User Data Functions (begin)
--------------------------------
	function ourDiffTensor2d(x, y, t)
		return	1, 0, 
				0, 1
	end
	
	-- discontinuous coefficient
	function jumpDiffTensor2d(x, y, t)
		if (math.abs (x-0.5)<1/4 and math.abs(y-0.5)<1/4) then
			return	RAepsilon, 0, 0, 1
		else 
			return	1, 0, 0, RAepsilon
		end
	end
	
	function CreateRotatedAnisotropyMatrix2d(alpha, epsilon)
		local sinalpha = math.sin(alpha)
		local cosalpha = math.cos(alpha)
		RAmat = ConstUserMatrix2d()
		print((sinalpha*sinalpha + epsilon*cosalpha*cosalpha)..", "..(1-epsilon)*sinalpha*cosalpha)
		print((1-epsilon)*sinalpha*cosalpha..", "..epsilon*sinalpha*sinalpha + cosalpha*cosalpha)
		RAmat:set_entry(0, 0, sinalpha*sinalpha + epsilon*cosalpha*cosalpha)
		RAmat:set_entry(0, 1, (1-epsilon)*sinalpha*cosalpha)
		RAmat:set_entry(1, 0, (1-epsilon)*sinalpha*cosalpha)
		RAmat:set_entry(1, 1, epsilon*sinalpha*sinalpha + cosalpha*cosalpha)
		return RAmat
	end		
	
	function ourVelocityField2d(x, y, t)
		return	0, 0 -- 150 ist grenzwertig
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

	
	function ourDirichletBnd2d(x, y, t)
		local s = 2*math.pi
		--return true, math.sin(s*x) + math.sin(s*y)
		return true, 0
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
		return true, 0
	end
--------------------------------
-- User Data Functions (end)
--------------------------------

-- create Instance of a Domain
dom = Domain()
if LoadDomain(dom, gridName) == false then
print("Loading Domain failed.")
exit()
end

-- create Refiner
my_assert(numPreRefs < numRefs, "numPreRefs must be smaller than numRefs");

refiner = GlobalDomainRefiner(dom)
for i=1,numPreRefs do
refiner:refine()
end

-- Distribute the domain to all involved processes
my_assert(DistributeDomain(dom) == true, "Error while Distributing Grid.")

-- Perform post-refine
for i=numPreRefs+1,numRefs do
refiner:refine()
end

-- get subset handler
sh = dom:get_subset_handler()
my_assert(sh:num_subsets() == 2, "Domain must have 2 Subsets for this problem.")
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- write grid to file for test purpose
-- SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-------------------------------------------
--  Setup User Functions
-------------------------------------------
if RAepsilon ~= 1.0 then
	diffusionMatrix = CreateRotatedAnisotropyMatrix2d(RAalpha, RAepsilon)
	-- diffusionMatrix = LuaUserMatrix("jumpDiffTensor"..dim.."d")
else
	diffusionMatrix = LuaUserMatrix("ourDiffTensor"..dim.."d")		
end


-- diffusionMatrix = ConstUserMatrix(1.0)

-- Velocity Field setup
velocityField = LuaUserVector("ourVelocityField"..dim.."d")
reaction = LuaUserNumber("ourReaction"..dim.."d")
rhs = LuaUserNumber("ourRhs"..dim.."d")
neumann = LuaBoundaryNumber("ourNeumannBnd"..dim.."d")
dirichlet = LuaBoundaryNumber("ourDirichletBnd"..dim.."d")

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------
if dim == 2 then
upwind = WeightedUpwind2d(); 
else
upwind = WeightedUpwind3d();
end
upwind:set_weight(0.0)
elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
my_assert(elemDisc:set_upwind(upwind), "could not set upwind")
elemDisc:set_diffusion_tensor(diffusionMatrix)
elemDisc:set_velocity_field(velocityField)
elemDisc:set_reaction(reaction)
elemDisc:set_source(rhs)

-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

--neumannDisc = FV1NeumannBoundary("Inner")
--neumannDisc:add(neumann, "c", "NeumannBoundary")

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(dirichlet, "c", "DirichletBoundary")

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
--domainDisc:add(neumannDisc)
domainDisc:add(dirichletBND)



-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)
linOp:set_dof_distribution(approxSpace:get_surface_dof_distribution())

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- set initial value
u:set_random(-1.0, 1.0)

-- init Operator and set dirichlet values in start iterate
print ("Assemble Operator and dirichlet ... ")
linOp:init_op_and_rhs(b)
linOp:set_dirichlet_values(u)
print ("done")

-- write matrix for test purpose
if bOutput then
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.vec")
end

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

-- create AMG ---
-----------------

bUseFAMG = 1
if bUseFAMG == 1 then
	print ("create FAMG... ")
	-- Testvectors for FAMG ---
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
		amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
		return amgTestvector	
	end
		
	function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
		local amgDirichlet0 = GridFunctionVectorWriterDirichlet02d()
		amgDirichlet0:init(dirichletBND, approxSpace)
		return amgDirichlet0
	end
	
	amg = FAMGPreconditioner()	
	amg:set_delta(0.5)
	amg:set_theta(0.95)
	amg:set_aggressive_coarsening(false)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.8)	
		
	-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
	testvectorwriter = CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	testvector = GridFunction(approxSpace)
	testvectorwriter:update(testvector)	
	amg:add_vector_writer(testvectorwriter, 1.0)
	amg:set_testvector_damps(1)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.8)
		
	-- amg:set_debug_level_get_ratings(4)
	-- amg:set_debug_level_coloring(4)
	-- amg:set_debug_level_communicate_prolongation(4)
	-- amg:set_debug_level_overlap(4,4)
	-- amg:set_debug_level_precalculate_coarsening(4)
	-- amg:set_debug_level_calculate_parent_pairs(4)
else
	print ("create AMG... ")
	amg = RSAMGPreconditioner()
	amg:enable_aggressive_coarsening_A(2)
end


vectorWriter = GridFunctionPositionProvider2d()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider2d(vectorWriter)
if bOutput then
amg:set_matrix_write_path("/Users/mrupp/matrices/")
end

amg:set_num_presmooth(2)
amg:set_num_postsmooth(2)
amg:set_cycle_type(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
if bUseNestedAMG then
	amg2 = RSAMGPreconditioner()
	amg2:set_base_solver(base)
	amg2:set_num_presmooth(3)
	amg2:set_num_postsmooth(3)
	amg2:set_cycle_type(1)
	amg2:set_presmoother(jac)
	amg2:set_postsmoother(jac)
	amg2:set_max_nodes_for_base(1000)
	amg2:set_max_fill_before_base(0.7)
	amg2:set_fsmoothing(true)
	amg2:set_epsilon_truncation(0)
	-- amg2:set_matrix_write_path("/Users/mrupp/matrices/2/")
	convCheck2 = StandardConvergenceCheck()
	convCheck2:set_maximum_steps(30)
	convCheck2:set_minimum_defect(1e-11)
	convCheck2:set_reduction(1e-12)
	convCheck2:set_verbose_level(false)
	linSolver2 = LinearSolver()
	linSolver2:set_preconditioner(amg2)
	linSolver2:set_convergence_check(convCheck2)
	amg:set_base_solver(linSolver2)
else
	amg:set_base_solver(base)
end

amg:set_max_levels(2)

amg:set_min_nodes_on_one_processor(10000)
amg:set_preferred_nodes_on_one_processor(10000)
amg:set_max_nodes_for_base(maxBase)
amg:set_max_fill_before_base(0.7)
amg:set_fsmoothing(true)
amg:set_epsilon_truncation(0)
amg:tostring()


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(9)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

print("done.")
-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(amg)
linSolver:set_convergence_check(convCheck)

-- Apply Solver

b2 = GridFunction(approxSpace)
b2:assign(b)

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator
print("Init operator (i.e. assemble matrix).")
linOp:init()

-- 2. init solver for linear Operator
print("Init solver for operator.")
linSolver:init(linOp)

-- 3. apply solver
if true then
print("Apply solver.")
tBefore = os.clock()
linSolver:apply_return_defect(u,b)
tSolve = os.clock()-tBefore
-- WriteGridFunctionToVTK(u, "Solution")

print("done")
else
for i = 1, 9 do
convCheck:set_maximum_steps(i)
srand(0)
u:set_random(-1.0, 1.0)
linOp:set_dirichlet_values(u)
b:assign(b2)
linSolver:apply_return_defect(u,b)
WriteGridFunctionToVTK(u, "Solution"..i)
end
end

if bOutput then
WriteGridFunctionToVTK(u, "Solution")
end

if false then
	print("Testing the AMG. Restart, then 6 Iteration steps:")
	b:assign(b2)
	u:set_random(-1.0, 1.0)
	linOp:set_dirichlet_values(u)
	convCheck:set_maximum_steps(6)
	linSolver:apply_return_defect(u,b)
	print("done. Now check:")
	amg:check(u, b)
	amg:check2(u, b)
end

printf = function(s,...)
	print(s:format(...))
end -- function


function fsize (file)
	local current = file:seek()      -- get current position
    local size = file:seek("end")    -- get file size
    file:seek("set", current)        -- restore position
    return size
end

if bFileOutput and GetProcessRank() == 0 then
	output = io.open("output_"..os.date("y%Ym%md%d")..".txt", "a")
	if fsize(output) == 0 then 
		output:write("procs")
		output:write("\tnumRefs")
		output:write("\tndofs")
		output:write("\tsteps")
		output:write("\tlastReduction")
		output:write("\ttSetupAmg [ms]")
		output:write("\tc_A")
		output:write("\tc_G")
		output:write("\tused Levels")
		output:write("\ttSolve [s]")
		output:write("\n")
	end
	output:write(GetNumProcesses())
	output:write("\t"..numRefs)
	output:write("\t"..amg:get_level_information(0):get_nr_of_nodes())
	output:write("\t"..convCheck:step())
	output:write("\t"..convCheck:defect()/convCheck:previous_defect())
	output:write("\t"..amg:get_timing_whole_setup_ms())
	output:write("\t"..amg:get_operator_complexity())
	output:write("\t"..amg:get_grid_complexity())
	output:write("\t"..amg:get_used_levels())
	output:write("\t"..tSolve)
	output:write("\n")
	print(s)
	-- else
end
