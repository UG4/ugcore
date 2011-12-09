----------------------------------------------------------
--
--   Lua - Script to perform the Laplace-Problem
--
--   Author: Martin Rupp / Andreas Vogel
--
----------------------------------------------------------
-- SetOutputProcessRank(0)
SetOutputProfileStats(false)

GetLogAssistant():enable_file_output(true, "ug_log"..GetProcessRank())
SetOutputProcessRank(-1)

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
	-- gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")	
	gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	-- gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 1)

-- choose number of pre-Refinements (before sending grid onto different processes)	
numPreRefs = util.GetParamNumber("-numPreRefs", math.min(5, numRefs-2))

maxBase = util.GetParamNumber("-maxBase", 300)
maxLevels = util.GetParamNumber("-maxLevels", 20)


RAepsilon = util.GetParamNumber("-RAepsilon", 1)
RAalpha = util.GetParamNumber("-RAalpha", 0)

epsx = util.GetParamNumber("-epsx", 1)
epsy = util.GetParamNumber("-epsy", 1)

bCheck = util.HasParamOption("-bCheck")

bWriteStats = util.HasParamOption("-bWriteStats")
bWriteMat = util.HasParamOption("-bWriteMat")
bRSAMG = util.HasParamOption("-RSAMG") 

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
print("    epsx = "..epsx)
print("    epsy = "..epsy)

function writeln(...)
	write(...)
	write("\n")
end

--------------------------------
-- User Data Functions (begin)
--------------------------------
	function cbDirichletBnd2d(x, y, t)
		return true, 5.9991		
	end
	function cbDirichletBnd3d(x, y, z, t)
		return true, 5.9991		
	end	
	dirchletBnd = LuaBoundaryNumber("cbDirichletBnd"..dim.."d")
	
	
	function cbSinRhs2d(x, y, t)
		local s = 2*math.pi
		return s*s*(math.sin(s*x) + math.sin(s*y))
	end
	sinRhs2d = LuaUserNumber("cbSinRhs2d")
	
	
	function cbSinDirichletBnd2d(x, y, t)
		local s = 2*math.pi
		return true, math.sin(s*x) + math.sin(s*y)
	end
	sinDirchletBnd2d = LuaBoundaryNumber("cbSinDirichletBnd2d")
	
	-- anisotropic diffusion in corners
	function cbAnisoDiffTensor2d(x, y, t)
		local fac1 = 1/1+exp(alpha*(x+y-1)) 
		local fac2 = 1- fac1
		return	fac1*epsx+fac2, 0, 
				0, fac1+fac2*epsy
	end
	anisoDiffTensor2d = LuaUserMatrix("cbAnisoDiffTensor2d")

	
	function CreateRotatedAnisotropyMatrix2d(alpha, epsilon)
		local sinalpha = math.sin(alpha)
		local cosalpha = math.cos(alpha)
		RAmat = ConstUserMatrix()
		-- print(sinalpha.." "..cosalpha.." "..epsilon)
		-- print((sinalpha*sinalpha + epsilon*cosalpha*cosalpha)..", "..(1-epsilon)*sinalpha*cosalpha)
		-- print((1-epsilon)*sinalpha*cosalpha..", "..epsilon*sinalpha*sinalpha + cosalpha*cosalpha)
		RAmat:set_entry(0,0,sinalpha*sinalpha + epsilon*cosalpha*cosalpha)
		RAmat:set_entry(0,1,(1-epsilon)*sinalpha*cosalpha)
		RAmat:set_entry(1,0,(1-epsilon)*sinalpha*cosalpha)
		RAmat:set_entry(1,1, epsilon*sinalpha*sinalpha + cosalpha*cosalpha)
		return RAmat		
	end	


	-- hedgehog diffusion 
	function cbHedgehogDiffTensor2d(x, y, t)
	
		if x<0 then 
			if y<0 then return 1.0, 0, 0, 1.0 
			else return epsx, 0, 0, 1.0 end
		else
			if y<0 then return 1.0, 0, 0, epsy 
			else return 1.0, -1.0, -1.0, 1.0 end			
		end
		
		-- should never happen
		return 0, 0, 0, 0 
	end
	
	
	-- discontinuous coefficient
	function cbJumpDiffTensor2d(x, y, t)
		if (math.abs (x)<0.5 and math.abs(y)<0.5) then
			return	epsx, 0, 0, epsy
		else 
			return	1, 0, 0, 1
		end
	end	
	

if dim == 3 then
	diffusionMatrix =  util.CreateConstUserMatrix3d(1, 0, 0,
													0, 1, 0,
													0, 0, 1)
	velocityField = util.CreateConstUserVector3d(0, 0, 0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd

else
	problem = "rotatedAniso"
	-------------------------------------------
	--  Setup User Functions
	-------------------------------------------
	if problem == "rotatedAniso" then
	diffusionMatrix = CreateRotatedAnisotropyMatrix2d(RAalpha, RAepsilon)
	velocityField = ConstUserVector2d(0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd
	end
	
	if problem == "hedgehog" then
	diffusionMatrix = LuaUserMatrix("cbHedgehogDiffTensor2d")
	velocityField = ConstUserVector2d(0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd
	end
	
	if problem == "jump" then
	diffusionMatrix = LuaUserMatrix("cbJumpDiffTensor2d")
	velocityField = ConstUserVector2d(0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd
	end
end
--------------------------------
-- User Data Functions (end)
--------------------------------
tBefore = os.clock()

-- create Instance of a Domain
dom = Domain()
if LoadDomain(dom, gridName) == false then
print("Loading Domain failed.")
exit()
end

-- create Refiner
assert(numPreRefs < numRefs, "numPreRefs must be smaller than numRefs");

refiner = GlobalDomainRefiner(dom)
for i=1,numPreRefs do
refiner:refine()
end

-- Distribute the domain to all involved processes
assert(DistributeDomain(dom) == true, "Error while Distributing Grid.")

-- Perform post-refine
for i=numPreRefs+1,numRefs do
refiner:refine()
end

tGrid = os.clock()-tBefore

-- get subset handler
sh = dom:subset_handler()
assert(sh:num_subsets() == 2, "Domain must have 2 Subsets for this problem.")
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- write grid to file for test purpose
-- SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)


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
elemDisc:set_upwind(upwind)
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
linOp:set_dof_distribution(approxSpace:surface_dof_distribution())

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
if bMatOutput then
SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")
SaveVectorForConnectionViewer(b, "Rhs.vec")
end

-- create algebraic Preconditioners

print ("create preconditioners... ")
jac = Jacobi()
jac:set_damp(0.66)
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

if bRSAMG == false then
	print ("create FAMG... ")
	-- Testvectors for FAMG ---
	--------------------------	
	function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
		local amgTestvector;
		amgTestvector = GridFunctionVectorWriter()
		amgTestvector:set_reference_grid_function(gridfunction)
		amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
		return amgTestvector	
	end
		
	function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
		local amgDirichlet0 = GridFunctionVectorWriterDirichlet0()
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
		
	if bWriteMat then
		amg:write_testvectors(true)
	end
	
	if true then
		amg:set_external_coarsening(true)
		amg:set_parallel_coarsening(GetColorCoarsening())
		-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
		-- amg:set_parallel_coarsening(GetRS3Coarsening())
	end
	amg:set_testvectorsmoother(jac)
	
	-- amg:set_debug_level_overlap(4, 4)
	-- amg:set_use_precalculate(false)
	-- amg:set_debug_level_get_ratings(4)
	-- amg:set_debug_level_coloring(4)
	-- amg:set_debug_level_communicate_prolongation(4)
	-- amg:set_debug_level_overlap(4,4)
	-- amg:set_debug_level_calculate_parent_pairs(2)
	-- amg:set_debug_level_precalculate_coarsening(4)
	-- amg:set_debug_level_calculate_parent_pairs(4)
else
	print ("create AMG... ")
	amg = RSAMGPreconditioner()
	-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
	-- amg:set_parallel_coarsening(GetColorCoarsening())
	amg:set_parallel_coarsening(GetRS3Coarsening())
	-- amg:set_parallel_coarsening(GetSimpleParallelCoarsening())
	-- amg:enable_aggressive_coarsening_A(2)
end


vectorWriter = GridFunctionPositionProvider()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider(vectorWriter)
if bWriteMat then
amg:set_matrix_write_path("/Users/mrupp/matrices/")
end

amg:set_num_presmooth(2)
amg:set_num_postsmooth(2)
amg:set_cycle_type(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
amg:set_max_levels(maxLevels)

amg:set_min_nodes_on_one_processor(50000)
-- amg:set_preferred_nodes_on_one_processor(1000)
amg:set_max_nodes_for_base(maxBase)
amg:set_max_fill_before_base(0.5)
amg:set_fsmoothing(true)
amg:set_epsilon_truncation(0)
amg:tostring()


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

print("done.")
-- create Linear Solver

if util.HasParamOption("-cg") then
	linSolver = CG()
else
	linSolver = LinearSolver()
end
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
tBefore = os.clock()
linOp:init()
tAssemble = os.clock()-tBefore

-- 2. init solver for linear Operator
print("Init solver for operator.")
linSolver:init(linOp)

-- 3. apply solver


if bCheck then
amg:write_interfaces()
convCheck:set_maximum_steps(4)
linSolver:apply_return_defect(u,b)
SaveVectorForConnectionViewer(b, "b.vec")
SaveVectorForConnectionViewer(u, "u.vec")
amg:check(u,b)
SaveVectorForConnectionViewer(b, "Rhs2.vec")
SaveVectorForConnectionViewer(u, "u2.vec")

else

print("Apply solver.")
tBefore = os.clock()
linSolver:apply_return_defect(u,b)
tSolve = os.clock()-tBefore
print("done")

end


printf = function(s,...)
	print(s:format(...))
end -- function

formatf = function(s, ...)
	return s:format(...)
end

function fsize (file)
	local current = file:seek()      -- get current position
    local size = file:seek("end")    -- get file size
    file:seek("set", current)        -- restore position
    return size
end

if not bCheck then
	if bFileOutput and GetProcessRank() == 0  then
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
			output:write("\ttGrid [s]")
			output:write("\ttAssemble [s]")
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
		output:write("\t"..tGrid)
		output:write("\t"..tAssemble)
		output:write("\n")
		print(s)
		-- else
	end
	
	
	print("NumProcesses: "..GetNumProcesses())
	print("numRefs: "..numRefs)
	print("nr of nodes: "..amg:get_level_information(0):get_nr_of_nodes())
	print("steps: "..convCheck:step())
	print("last reduction: "..convCheck:defect()/convCheck:previous_defect())
	print("tSetupAmg [ms]: "..amg:get_timing_whole_setup_ms())
	print("c_A: "..amg:get_operator_complexity())
	print("c_G: "..amg:get_grid_complexity())
	print("used Levels: "..amg:get_used_levels())
	print("tSolve [s]: "..tSolve)
	print("tGrid [s]: "..tGrid)
	print("tAssemble [s]: "..tAssemble)
	
	
	LI = amg:get_level_information(0)
	printf("Level 0: number of nodes: %d. fill in %.2f%%, nr of interface elements: %d (%.2f%%)", LI:get_nr_of_nodes(), LI:get_fill_in()*100, LI:get_nr_of_interface_elements(), LI:get_nr_of_interface_elements()/LI:get_nr_of_nodes()*100)
	for i = 1, amg:get_used_levels()-1 do
		LI = amg:get_level_information(i)
		printf("Level %d: creation time: %.2f ms. number of nodes: %d. fill in %.2f%%. coarsening rate: %.2f%%, nr of interface elements: %d (%.2f%%)", i, 
			LI:get_creation_time_ms(),  LI:get_nr_of_nodes(), LI:get_fill_in()*100,
			LI:get_nr_of_nodes()*100/amg:get_level_information(i-1):get_nr_of_nodes(), LI:get_nr_of_interface_elements(), LI:get_nr_of_interface_elements()/LI:get_nr_of_nodes()*100)
			
	end
	
	if GetProfilerAvailable() == true then
		create_levelPN = GetProfileNode("c_create_AMG_level")
		
		function PrintParallelProfileNode(name)
			pn = GetProfileNode(name)
			t = pn:get_avg_total_time_ms()/to100 * 100
			tmin = ParallelMin(t)
			tmax = ParallelMax(t)
			printf("%s:\n%.2f %%, min: %.2f %%, max: %.2f %%", name, t, tmin, tmax)
		end
		
		if true then
			print(create_levelPN:call_tree())
			print(create_levelPN:child_self_time_sorted())
		end
		if false then			
			to100 = create_levelPN:get_avg_total_time_ms()
			PrintParallelProfileNode("create_OL2_matrix")
			PrintParallelProfileNode("CalculateTestvector")
			PrintParallelProfileNode("CreateSymmConnectivityGraph")
			PrintParallelProfileNode("calculate_all_possible_parent_pairs")
			PrintParallelProfileNode("color_process_graph")
			PrintParallelProfileNode("FAMG_recv_coarsening_communicate")
			PrintParallelProfileNode("update_rating")
			PrintParallelProfileNode("precalculate_coarsening")
			PrintParallelProfileNode("send_coarsening_data_to_processes_with_higher_color")
			PrintParallelProfileNode("communicate_prolongation")
			PrintParallelProfileNode("create_new_index")
			PrintParallelProfileNode("create_parent_index")
			PrintParallelProfileNode("create_interfaces")
			PrintParallelProfileNode("create_fine_marks")
			PrintParallelProfileNode("create_galerkin_product")
			PrintParallelProfileNode("CalculateNextTestvector")
		end
	
	else
		print("Profiler not available.")
	end	
end
