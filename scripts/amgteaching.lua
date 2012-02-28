----------------------------------------------------------
--
--   Lua - Script performing several test problems
--
--   Author: 	Arne Naegel
--   
-- 				based on laplace.lua by Andreas Vogel 
----------------------------------------------------------

ug_load_script("ug_util.lua")

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
dim = util.GetParamNumber("-dim", 2)

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = "unit_square_tri_four_dirichlet_nodes.ugx"
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 2)
numRefs    = util.GetParamNumber("-numRefs",    5)

-- parallelisation related stuff
-- way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect") -- [grid2d | bisect | metis]

-- amount of output
verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: this sets
						    -- 'fetiSolver:set_debug(dbgWriter)'

-- parameters for diffusion coefficients
epsx = util.GetParamNumber("-epsx", 1.0)
epsy = util.GetParamNumber("-epsy", 1.0)

uvel = util.GetParamNumber("-uvel", 0.0)
vvel = util.GetParamNumber("-vvel", 0.0)


alpha = util.GetParamNumber("-alpha", 0.0)
ilubeta = util.GetParamNumber("-ilu-beta", 0.0)

-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)

print("    epsx = " .. epsx)
print("    epsy = " .. epsy)
print("    alpha = " .. alpha)
print("    ilu-beta = " .. ilubeta)
--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

--------------------------------
-- User Data Functions 2d (begin)
--------------------------------
	
	
	-- A) DIFFUSION
	-- standard (anisotropic) diffusion 
	function stdDiffTensor2d(x, y, t)
	--print (epsx.." "..epsy)
	return	epsx, 0, 0, epsy end
	
	-- anisotropic diffusion in corners
	function anisoDiffTensor2d(x, y, t)
		local fac1 = 1/1+exp(alpha*(x+y-1)) 
		local fac2 = 1- fac1
		return	fac1*epsx+fac2, 0, 
				0, fac1+fac2*epsy
	end
	
	-- hedgehog diffusion 
	function hedgehogDiffTensor2d(x, y, t)
		if x<0.5 then 
			if y<0.5 then return 1.0, 0, 0, 1.0 
			else return epsx, 0, 0, 1.0 end
		else
			if y<0.5 then return 1.0, 0, 0, epsy 
			else return 1.0, -1.0, -1.0, 1.0 end
		end
		
		-- should never happen
		return 0, 0, 0, 0 
	end

	
	
	-- discontinuous coefficient
	function jumpDiffTensor2d(x, y, t)
		if (math.abs (x-0.5)<1/4 and math.abs(y-0.5)<1/4) then
			return	epsx, 0, 0, epsy
		else 
			return	1, 0, 0, 1
		end
	end
	
	-- B) VELOCITIES
	function zeroVelocityField2d(x, y, t)
		return	0, 0
	end
	
	function constVelocityField2d(x, y, t)
		return	uvel, vvel
	end

	function recirculatingVelocityField2d(x, y, t)
		local urecirc = 4.0*x*(1-x)*(1-2*y) 
		local vrecirc = -4.0*y*(1-y)*(1-2*x)
		return	urecirc, vrecirc
	end
	
	-- C) REACTIONS
	function zeroReaction2d(x, y, t)
		return	0
	end
	
	
	-- D) SOURCES 
	function stdSource2d(x, y, t)
		local s = 2*math.pi
		return	s*s*(math.sin(s*x) + math.sin(s*y))
		--return -2*y
--		return 0;
	end
	
	-- E) BOUNDARY
	function stdNeumannBnd2d(x, y, t)
		--local s = 2*math.pi
		--return -s*math.cos(s*x)
		return true, -2*x*y
	end
	
	function stdDirichletBnd2d(x, y, t)
	--	local s = 2*math.pi
	--	return true, math.sin(s*x) + math.sin(s*y)
	return true, 0
		--return true, x*x*y
--		return true, 2.5
	end
--------------------------------
-- User Data Functions 2d (end)
--------------------------------

	
----------------------------------
-- User Data Functions 3d (begin)
----------------------------------
	

-- 

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
dom = Domain()

-- load domain
print("Load Domain from File.")
if LoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end

-- create Refiner
print("Create Refiner")
if numPreRefs > numRefs then
	print("It must be choosen: numPreRefs <= numRefs");
	exit();
end

-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
refiner = GlobalDomainRefiner(dom)

-- Performing pre-refines
print("Perform (non parallel) pre-refinements of grid")
for i=1,numPreRefs do
	print( "PreRefinement step " .. i .. " ...")
	refiner:refine()
	print( "... done!")
end

-- get number of processes
numProcs = GetNumProcesses()

-- Distribute the domain to all involved processes
-- Since only process 0 loaded the grid, it is the only one which has to
-- fill a partitionMap (but every process needs one and has to return his map
-- by calling 'RedistributeDomain()', even if in this case the map is empty
-- for all processes but process 0).
print("Distribute domain with 'distributionType' = " .. distributionType .. "...")
partitionMap = PartitionMap()

if GetProcessRank() == 0 then
	if distributionType == "bisect" then
		util.PartitionMapBisection(partitionMap, dom, numProcs)
		
	elseif distributionType == "grid2d" then
		local numNodesX, numNodesY = util.FactorizeInPowersOfTwo(numProcs / numProcsPerNode)
		util.PartitionMapLexicographic2D(partitionMap, dom, numNodesX,
										 numNodesY, numProcsPerNode)

	elseif distributionType == "metis" then
		util.PartitionMapMetis(partitionMap, dom, numProcs)
										 
	else
	    print( "distributionType not known, aborting!")
	    exit()
	end
	
-- save the partition map for debug purposes
	if verbosity >= 1 then
		SavePartitionMap(partitionMap, dom, "partitionMap_p" .. GetProcessRank() .. ".ugx")
	end
end

if RedistributeDomain(dom, partitionMap, true) == false then
	print("Redistribution failed. Please check your partitionMap.")
	exit()
end
print("... domain distributed!")

--------------------------------------------------------------------------------
-- end of partitioning
--------------------------------------------------------------------------------

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	print( "Refinement step " .. i .. " ...")
	refiner:refine()
	print( "... done!")
end

-- get subset handler
sh = dom:subset_handler()
if sh:num_subsets() ~= 2 then 
print("Domain must have 2 Subsets for this problem.")
exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- write grid to file for test purpose
if verbosity >= 1 then
	refinedGridOutName = "refined_grid_p" .. GetProcessRank() .. ".ugx"
	print("saving domain to " .. refinedGridOutName)
	if SaveDomain(dom, refinedGridOutName) == false then
		print("Saving of domain to " .. refinedGridOutName .. " failed. Aborting.")
		    exit()
	end
	
	hierarchyOutName = "hierachy_p" .. GetProcessRank() .. ".ugx"
	print("saving hierachy to " .. hierarchyOutName)
	if SaveGridHierarchy(dom:grid(), hierarchyOutName) == false then
		print("Saving of domain to " .. hierarchyOutName .. " failed. Aborting.")
		    exit()
	end
end

print("NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--  Order indices using Cuthill-McKee
-- OrderCuthillMcKee(approxSpace, true);


--Order indices lexicographically
OrderLex(approxSpace, "ul");

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.66)
gs = GaussSeidel()
bgs = BackwardGaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
--ilu:set_debug(dbgWriter)
ilu:set_beta(ilubeta)
ilut = ILUT()

------------------------------------------
-- solvers for alebraically smooth error
------------------------------------------
-- create convergence check (solver)
convCheckSmoother = StandardConvergenceCheck()
convCheckSmoother:set_maximum_steps(5)
convCheckSmoother:set_reduction(1e-12)
convCheckSmoother:set_minimum_defect(1e-14)

-- create Jacobi Solver
jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheckSmoother)

-- create Gauss-Seidel Solver
gsSolver = LinearSolver()
gsSolver:set_preconditioner(gs)
gsSolver:set_convergence_check(convCheckSmoother)

-- create SGS Solver
sgsSolver = LinearSolver()
sgsSolver:set_preconditioner(sgs)
sgsSolver:set_convergence_check(convCheckSmoother)

-- create ILU Solver
iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheckSmoother)


-- select smoother/solver

smoother={}
--smoother.solver = jacSolver
--smoother.name = "Jacobi"
--smoother.pre = jac
--smoother.post = jac
smoother.solver = gsSolver
smoother.name = "GS"
smoother.pre = gs
smoother.post = bgs
--smoother.name = "SGS"
--smoother.solver = sgsSolver
--smoother.pre = sgs
--smoother.post = sgs
--smoother.solver = iluSolver
--smoother.name = "ILU"
--smoother.pre = ilu
--smoother.post = ilu

-- create (exact) base solver
-----------------

baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(1000)
baseConvCheck:set_minimum_defect(1e-11)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose_level(false)
--base = LU()
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(ilu)



-- create GMG ---
-----------------

-- Geometric Multi Grid
gmg = GeometricMultiGrid(approxSpace)

gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_parallel_base_solver(true)
gmg:set_smoother(smoother.pre)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(1)
gmg:set_num_postsmooth(1)
if activateDbgWriter >= 1 then
gmg:set_debug(dbgWriter)
end

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




------------------------------------------
-- real solvers
------------------------------------------
-- create convergence check (solver)
convCheckSolver = StandardConvergenceCheck()
convCheckSolver:set_maximum_steps(50)
convCheckSolver:set_minimum_defect(1e-14)
convCheckSolver:set_reduction(1e-12)
--convCheck:set_verbose_level(true)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheckSolver)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheckSolver)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(jac)
bicgstabSolver:set_convergence_check(convCheckSolver)
solver = linSolver
--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------	
function CreateProblemDependentElemDisc(problem, dim)
print ("Setting up Assembling")

-- Select upwind
if dim == 2 then 

--upwind = NoUpwind2d()
upwind = FullUpwind2d()
--upwind = WeightedUpwind2d(); upwind:set_weight(0.0)
--upwind = PartialUpwind2d()
else print("Dim not supported for upwind"); exit() end

elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
if elemDisc:set_upwind(upwind) == false then exit() end
elemDisc:set_diffusion_tensor(problem.diffusion)
elemDisc:set_velocity_field(problem.velocity)
elemDisc:set_reaction(problem.reaction)
elemDisc:set_source(problem.source)
return elemDisc
end

--------------------------------------------------------------------------------
--  Setup problems
--------------------------------------------------------------------------------

problem = {}

-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDiffTensor2d")

-- Problem 1: Poisson's eqn. (w/ anisotropy) 
if false then
	problem["Poisson"] = {}
	problem["Poisson"].desc = "Problem 1: Poisson's eqn. (w/ anisotropy)"
	problem["Poisson"].diffusion = LuaUserMatrix("stdDiffTensor2d")
	problem["Poisson"].velocity = LuaUserVector("zeroVelocityField2d")
	problem["Poisson"].reaction = LuaUserNumber("zeroReaction2d")
	problem["Poisson"].source = LuaUserNumber("stdSource2d")
	problem["Poisson"].neumann = LuaBoundaryNumber("stdNeumannBnd2d")
	problem["Poisson"].dirichlet = LuaBoundaryNumber("stdDirichletBnd2d")
	
	problem["Poisson"].configs = {}
	problem["Poisson"].configs[1] = {["epsx"]= 1.0, ["epsy"]= 1.0}
	problem["Poisson"].configs[2] = {["epsx"]= 1.0, ["epsy"]= 1e-1}
	problem["Poisson"].configs[3] = {["epsx"]= 1.0, ["epsy"]= 0.075}
	problem["Poisson"].configs[4] = {["epsx"]= 1.0, ["epsy"]= 5e-2}
	problem["Poisson"].configs[5] = {["epsx"]= 1.0, ["epsy"]= 2e-2}
	problem["Poisson"].configs[6] = {["epsx"]= 1.0, ["epsy"]= 1e-2}
	problem["Poisson"].configs[7] = {["epsx"]= 1.0, ["epsy"]= 5e-3}
	problem["Poisson"].configs[8] = {["epsx"]= 1.0, ["epsy"]= 1e-3}
	problem["Poisson"].configs[9] = {["epsx"]= 1.0, ["epsy"]= 5e-4}
	problem["Poisson"].configs[10] = {["epsx"]= 1.0, ["epsy"]= 1e-4}
	problem["Poisson"].configs[11] = {["epsx"]= 1.0, ["epsy"]= 1e-5}
	problem["Poisson"].configs[12] = {["epsx"]= 1.0, ["epsy"]= 1e-6}
	--problem["Poisson"].configs[8] = {["epsy"]= 1.0, ["epsx"]= 1e-1}
	--problem["Poisson"].configs[9] = {["epsy"]= 1.0, ["epsx"]= 1e-2}
	--problem["Poisson"].configs[10] = {["epsy"]= 1.0, ["epsx"]= 1e-3}
	--problem["Poisson"].configs[11] = {["epsy"]= 1.0, ["epsx"]= 1e-4}
	--problem["Poisson"].configs[12] = {["epsy"]= 1.0, ["epsx"]= 1e-5}
	--problem["Poisson"].configs[13] = {["epsy"]= 1.0, ["epsx"]= 1e-6}
end

-- Problem 2: Convection-diffusion (w/ anisotropy) 
if false then
	problem["ConvDiff"] = {}
	problem["ConvDiff"].desc = "Problem 2: Convection-diffusion (w/ anisotropy)"
	problem["ConvDiff"].diffusion = LuaUserMatrix("stdDiffTensor2d")
	problem["ConvDiff"].velocity = LuaUserVector("constVelocityField2d")
	problem["ConvDiff"].reaction = LuaUserNumber("zeroReaction2d")
	problem["ConvDiff"].source = LuaUserNumber("stdSource2d")
	problem["ConvDiff"].neumann = LuaBoundaryNumber("stdNeumannBnd2d")
	problem["ConvDiff"].dirichlet = LuaBoundaryNumber("stdDirichletBnd2d")
	
	problem["ConvDiff"].configs = {}
	problem["ConvDiff"].configs[1] = {["epsx"]= 1.0, ["epsy"]= 1.0}
	problem["ConvDiff"].configs[2] = {["epsx"]= 1e-1, ["epsy"]= 1e-1}
	problem["ConvDiff"].configs[3] = {["epsx"]= 1e-2, ["epsy"]= 1e-2}
	problem["ConvDiff"].configs[4] = {["epsx"]= 1e-3, ["epsy"]= 1e-3}
	problem["ConvDiff"].configs[5] = {["epsx"]= 1e-4, ["epsy"]= 1e-4}
	problem["ConvDiff"].configs[6] = {["epsx"]= 1e-5, ["epsy"]= 1e-5}
	problem["ConvDiff"].configs[7] = {["epsx"]= 1e-6, ["epsy"]= 1e-6}
end

-- Problem 3: Recirculating flow (w/ anisotropy) 
if false then
	problem["Recirc"] = {}
	problem["Recirc"].desc = "Problem 2: Recirculating flow (w/ anisotropy)"
	problem["Recirc"].diffusion = LuaUserMatrix("stdDiffTensor2d")
	problem["Recirc"].velocity = LuaUserVector("recirculatingVelocityField2d")
	problem["Recirc"].reaction = LuaUserNumber("zeroReaction2d")
	problem["Recirc"].source = LuaUserNumber("stdSource2d")
	problem["Recirc"].neumann = LuaBoundaryNumber("stdNeumannBnd2d")
	problem["Recirc"].dirichlet = LuaBoundaryNumber("stdDirichletBnd2d")
	
	problem["Recirc"].configs = {}
	problem["Recirc"].configs[1] ={["epsx"]= 1.0, ["epsy"]= 1.0}
	problem["Recirc"].configs[1] = {["epsx"]= 1.0, ["epsy"]= 1.0}
	problem["Recirc"].configs[2] = {["epsx"]= 1e-1, ["epsy"]= 1e-1}
	problem["Recirc"].configs[3] = {["epsx"]= 1e-2, ["epsy"]= 1e-2}
	problem["Recirc"].configs[4] = {["epsx"]= 1e-3, ["epsy"]= 1e-3}
	problem["Recirc"].configs[5] = {["epsx"]= 1e-4, ["epsy"]= 1e-4}
	problem["Recirc"].configs[6] = {["epsx"]= 1e-5, ["epsy"]= 1e-5}
	problem["Recirc"].configs[7] = {["epsx"]= 1e-6, ["epsy"]= 1e-6}
end

-- Problem 4: Hedgehog diffusion (w/ anisotropy, four different subdomains) 
if true then
	problem["Hedgehog"] = {}
	problem["Hedgehog"].desc = "Problem 4: Hedgehog-type Poisson's eqn. (w/ anisotropy)"
	problem["Hedgehog"].diffusion = LuaUserMatrix("hedgehogDiffTensor2d")
	problem["Hedgehog"].velocity = LuaUserVector("zeroVelocityField2d")
	problem["Hedgehog"].reaction = LuaUserNumber("zeroReaction2d")
	problem["Hedgehog"].source = LuaUserNumber("stdSource2d")
	problem["Hedgehog"].neumann = LuaBoundaryNumber("stdNeumannBnd2d")
	problem["Hedgehog"].dirichlet = LuaBoundaryNumber("stdDirichletBnd2d")
	
	problem["Hedgehog"].configs = {}
	problem["Hedgehog"].configs[1] = {["epsx"]= 1.0, ["epsy"]= 1.0}
	problem["Hedgehog"].configs[2] = {["epsx"]= 1e-1, ["epsy"]= 1e-1}
	problem["Hedgehog"].configs[3] = {["epsx"]= 1e-2, ["epsy"]= 1e-2}
	problem["Hedgehog"].configs[4] = {["epsx"]= 1e-3, ["epsy"]= 1e-3}
	problem["Hedgehog"].configs[5] = {["epsx"]= 1e-4, ["epsy"]= 1e-4}
	problem["Hedgehog"].configs[6] = {["epsx"]= 1e-5, ["epsy"]= 1e-5}
	problem["Hedgehog"].configs[7] = {["epsx"]= 1e-6, ["epsy"]= 1e-6}
end

-- Problem 5: Isolator problem (embeded squarent subdomains) 
if true then
	problem["Isolator"] = {}
	problem["Isolator"].desc = "Problem 4: Isolator problem (jumping coefficient)"
	problem["Isolator"].diffusion = LuaUserMatrix("jumpDiffTensor2d")
	problem["Isolator"].velocity = LuaUserVector("zeroVelocityField2d")
	problem["Isolator"].reaction = LuaUserNumber("zeroReaction2d")
	problem["Isolator"].source = LuaUserNumber("stdSource2d")
	problem["Isolator"].neumann = LuaBoundaryNumber("stdNeumannBnd2d")
	problem["Isolator"].dirichlet = LuaBoundaryNumber("stdDirichletBnd2d")
	
	problem["Isolator"].configs = {}
	problem["Isolator"].configs[1] = {["epsx"]= 1e+6, ["epsy"]= 1e+6}
	problem["Isolator"].configs[2] = {["epsx"]= 1e+4, ["epsy"]= 1e+4}
	problem["Isolator"].configs[3] = {["epsx"]= 1e+2, ["epsy"]= 1e+2}
	problem["Isolator"].configs[4] = {["epsx"]= 1e+1, ["epsy"]= 1e+1}
	problem["Isolator"].configs[5] = {["epsx"]= 1.0, ["epsy"]= 1.0}
	problem["Isolator"].configs[6] = {["epsx"]= 1e-1, ["epsy"]= 1e-1}
	problem["Isolator"].configs[7] = {["epsx"]= 1e-2, ["epsy"]= 1e-2}
	problem["Isolator"].configs[8] = {["epsx"]= 1e-4, ["epsy"]= 1e-4}
	problem["Isolator"].configs[9] = {["epsx"]= 1e-6, ["epsy"]= 1e-6}
end



--------------------------------------------------------------------------------
--  Run problems (with several configurations)
--------------------------------------------------------------------------------	

-- for each of the problems
for myname,myproblem in pairs(problem) do
	
	
	
	elemDisc = CreateProblemDependentElemDisc(myproblem, dim)
	print ("Running "..myproblem.desc.." ["..myname.."]")
	--util.PrintTable(myproblem)
	
	
	-- for each configuration
	print(myproblem.configs)
	for key,config in pairs(myproblem.configs) do
		print ("Running parameter set "..key..":")
		print ("  epsx="..config.epsx..", epsy="..config.epsy)
		epsx = config.epsx
		epsy = config.epsy
		
		
		-----------------------------------------------------------------
		--  Setup Neumann Boundary
		-----------------------------------------------------------------
		
		--neumannDisc = FV1NeumannBoundary("Inner")
		--neumannDisc:add(neumann, "c", "NeumannBoundary")
		
		-----------------------------------------------------------------
		--  Setup Dirichlet Boundary
		-----------------------------------------------------------------
		
		dirichletBND = DirichletBoundary()
		dirichletBND:add(myproblem.dirichlet, "c", "DirichletBoundary")
		
		-------------------------------------------
		--  Setup Domain Discretization
		-------------------------------------------
		
		domainDisc = DomainDiscretization(approxSpace)
		domainDisc:add(elemDisc)
		--domainDisc:add(neumannDisc)
		domainDisc:add(dirichletBND)
		
		--------------------------------------------------------------------------------
		--  Algebra
		--------------------------------------------------------------------------------
		print ("Setting up Algebra Solver")
		
		-- create operator from discretization
		linOp = AssembledLinearOperator()
		linOp:set_discretization(domainDisc)
		
		-- get grid function
		u = GridFunction(approxSpace)
		b = GridFunction(approxSpace)
		
		-- debug writer
		dbgWriter = GridFunctionDebugWriter(approxSpace)
		dbgWriter:set_vtk_output(false)
		

		-- update GMG ---
		-----------------
		gmg:set_discretization(domainDisc)	
		
		--------------------------------------------------------------------------------
		--  Apply solver - using method defined in 'operator_util.h',
		--  to get separate profiling for assemble and solve
		--------------------------------------------------------------------------------
		
		-- 0. Reset start solution
		u:set(0.0)
		
		-- 1. init operator
		print("Init operator (i.e. assemble matrix).")
		if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
		print("Could not assemble operator"); exit(); 
		end
		
		-- 1.b write matrix for test purpose
		if verbosity >= 1 then
		SaveMatrixForConnectionViewer(u, linOp, myname.."P"..key.."Stiffness.mat")
		SaveVectorForConnectionViewer(b, myname.."P"..key.."Rhs.vec")
		end
		
		-- 2. apply solver
		print("Apply solver.")
		if ApplyLinearSolver(linOp, u, b, solver) == false then
		print("Could not apply linear solver.");
		end
		
		if verbosity >= 1 then
		WriteGridFunctionToVTK(u, myname.."P"..key.."Solve"..smoother.name)
		end
	
		
		--------------------------------------------------------------------------------
		--  Apply smoother - this generates algebraically smooth errors
		--------------------------------------------------------------------------------
		
		-- 1. Init operator and modify right hand side
		u:set_random(-1.0, 1.0)
		print("Init operator (i.e. assemble matrix).")
		if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
		print("Could not assemble operator"); exit(); 
		end
		b:set(0.0)
		
		-- 2. apply solver
		print("Apply solver.")
		if ApplyLinearSolver(linOp, u, b, smoother.solver) == false then
		print("Could not apply linear solver.");
		end
		
		
		if verbosity >= 1 then
		WriteGridFunctionToVTK(u, myname.."P"..key.."Smooth"..smoother.name)
		end
	
	end -- for each configuration
end -- for each problem



--------------------------------------------------------------------------------
--  Print Profiling
--------------------------------------------------------------------------------

-- check if profiler is available
if GetProfilerAvailable() == true then
-- get node
pn = GetProfileNode("main")
--    pn2 = GetProfileNode("GMG_lmgc")
-- check if node is valid
if pn:is_valid() then
print(pn:call_tree())
--        print(pn2:total_time_sorted())
else
print("main is not known to the profiler.")
end
else
print("Profiler not available.")
end 
