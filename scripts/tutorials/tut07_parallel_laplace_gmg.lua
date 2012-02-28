--------------------------------------------------------------------------------
--	tut07_parallel_laplace_gmg.lua
--
-- In this tutorial a distributed domain is set up that has several vertical
-- interfaces. Those interfaces are created when a part of a grid is 
-- distributed to other processes and further refined. Several interfaces 
-- are e.g. created when a tree-like distribution is performed.
--------------------------------------------------------------------------------
-- include the basic util-methods.
ug_load_script("../ug_util.lua")

-- Get the command line parameters
dim = util.GetParamNumber("-dim", 2)

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2.ugx")
outFileNamePrefix = util.GetParam("-o", "distributed_domain_")
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
numMidRefs = util.GetParamNumber("-numMidRefs", 1)
numPostRefs = util.GetParamNumber("-numPostRefs", 1)
debug = util.HasParamOption("-debug")

-- We will save the created hierarchy to this file (with appended process id)
outHierarchyFilePrefix = util.GetParam("-oh", "hierarchy_on_proc_")

-- Create the domain and load a grid
dom = Domain()

if LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)


-- We'll create a global refiner and prerefine the domain before distribution starts.
refiner = GlobalDomainRefiner(dom)

-- perform pre-refinement
for i = 1, numPreRefs do
	print("refining...")
	refiner:refine()
end

--	This demonstration requires 4 processes
if GetNumProcesses() < 4 then
	print("At least 4 processes are required for this script. Aborting.")
	exit()
end

--	First we'll create a partition map. This is done on all processes, even
--	though the grid is currently only on process 0.
partitionMap = PartitionMap()

--	We'll also store the rank of the current process
procRank = GetProcessRank()

--	first we'll distribute sections of the domain to processes 0, 4, 8 and 12
--firstTargets = {0, 4, 8, 12}
firstTargets = {0, 2}

--	Since only process 0 has a grid at the moment, we'll only fill the
--	partition map on him.
if procRank == 0 then
	-- Iterate over all entries in firstTargets
	for i, target in ipairs(firstTargets) do
		partitionMap:add_target_proc(target)
	end
	
	-- The following call performs the load-balancing and partitions the
	-- domain into several sections, which will then be sent to the different
	-- processes.
	-- Note that several different load-balancers exist. In this example
	-- we'll simply use a regular-grid with 2x2 cells.
	--
	-- Note that the i-th target corresponds to the cell with indices cx, cy,
	-- where i = cy * numCellsY + cx	(numCellsX and numCellsY are specified
	-- during the call to PartitionDomain_RegularGrid).
	--
	-- For this algorithm it is crucial that the number of processes in the
	-- partitionMap match the number of cells specified for the load-balancer.
	-- The last parameter says that only surface-elements shall be partitioned
	PartitionDomain_RegularGrid(dom, partitionMap, 2, 1, true)
	
	-- We'll save the partition map. This should only be done for debugging.
	SaveGrid(dom:grid(), partitionMap:get_partition_handler(),
		"partitionMap_1_p" .. procRank .. ".ugx")
end


--	The partition map is set up. We can now distribute the domain.
--	RedistributeDomain takes three parameters: The domain, the
--	partitionMap, and a boolean indicating whether vertical interfaces
--	shall be created.
print("\nPERFORMING INITIAL DISTRIBUTION")
if RedistributeDomain(dom, partitionMap, true) == false then
	print("First redistribution failed. Please check your partitionMap.")
	exit()
end

if debug == true then
	if TestDomainInterfaces(dom) == false then
		print("Inconsistent grid layouts after first distribution. Aborting.")
		exit()
	end
end

--	The different target processes now all have their part of the grid.
--	we'll now perform the 'mid-refinement' before we'll further
--	distribute the domain.
-- perform mid-refinement
for i = 1, numMidRefs do
	print("refining...")
	refiner:refine()
end

--	We'll further distribute it using regular-grid-partitioning again.
--	Process 0 will distribute its grid among processes 0, 1, 2, 3.
--	Process 4 will distribute its grid among processes 4, 5, 6, 7.
--	...

partitionMap:clear()
if procRank % 2 == 0 and procRank < 4 then
	-- Specify the first process and the number of processes.
	-- This will add processes src, src+1, ..., src+n-1
	partitionMap:add_target_procs(procRank, 2)
	
	-- We again use regular-grid partitioning 
	PartitionDomain_RegularGrid(dom, partitionMap, 1, 2, true)
	
	-- again we'll save the partition map. This should only be done for debugging.
	SaveGrid(dom:grid(), partitionMap:get_partition_handler(),
		"partitionMap_2_p" .. procRank .. ".ugx")
end

--	The partition map is set up. We can now redistribute the domain
print("\nPERFORMING REDISTRIBUTION")
if RedistributeDomain(dom, partitionMap, true) == false then
	print("Second redistribution failed. Please check your partitionMap.")
	exit()
end

--	We can optionally test the interfaces here.
if debug == true then
	if TestDomainInterfaces(dom) == false then
		print("Inconsistent grid layouts after second distribution. Aborting.")
		exit()
	end
end

--	Now perform post-refinement
for i = 1, numPostRefs do
	print("refining...")
	refiner:refine()
end

-- Lets save the domain on each process
outFileName = outFileNamePrefix .. GetProcessRank() .. ".ugx"
if SaveDomain(dom, outFileName) == false then
	print("Saving of domain to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)


-- Now lets save the hierarchy on each process
-- The SaveGridHierarchy routine directly works on the domains grid.
-- We can access the grid of a domain through its grid() member method.
--
-- SaveGridHierarchy outputs a grid, where each level is assigned to a subset.
-- Original subsets are not contained in that file.
outFileName = outHierarchyFilePrefix .. GetProcessRank() .. ".ugx"
if SaveGridHierarchy(dom:grid(), outFileName) == false then
	print("Saving of grid-hierarch to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Again everything went fine.
print("Saved hierarchy to " .. outFileName)

--------------------------------------------------------------------------------
--	New code starts here

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}

if util.CheckSubsets(dom, neededSubsets) == false then
	print("Some subsets are incorrect. Aborting."); exit();
end

print("subsets checked")
ug_load_script("tut04_2_disc_laplace.lua")

-- Create surface functions (vectors) for Au=b (where A=linOp) and initialize them
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- Using the AssembleLaplace method which we created in
-- tut04_2_disc_laplace.lua, we now create the linear operator.
-- Note that we use the default callbacks defined in tut04_2_disc_laplace.lua
-- by passing nil. One could instead simply not specify the arguments
-- (this defaults to nil).
-- You could of course create your callbacks analogous to tutorial 3 and
-- pass their names (as strings) to AssembleLaplace.
-- Make sure that you use the ones with the right dimension!
linOp, approxSpace, domainDisc, dirichletBnd = AssembleLaplace(dom, "Inner", "Boundary", b, nil, nil, nil, nil)

approxSpace:print_layout_statistic()

print("Domain created.")

-- set initial value
u:set(0)
linOp:set_dirichlet_values(u)

-----------------
-- create GMG 
-----------------

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.8)

-- create Base Solver
base = LU()
	
-- create Geometric Multi Grid
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_smoother(jac)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(10000)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

-- create Linear Solver
solver = LinearSolver()
solver:set_preconditioner(gmg)
solver:set_convergence_check(convCheck)

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator
-- ( This is already done in AssembleLaplace, therefore skipped )
--print("Init operator (i.e. assemble matrix).")
--linOp:init()

-- 2. init solver for linear Operator
print("Init solver for operator.")
solver:init(linOp)

-- 3. apply solver
print("Apply solver.")
solver:apply_return_defect(u,b)

-- 4. print solution to vtk file format
WriteGridFunctionToVTK(u, "Solution")

--	New code ends here
--------------------------------------------------------------------------------

-- we're done.
print("")
print("done")
