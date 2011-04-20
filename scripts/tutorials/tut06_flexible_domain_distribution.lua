--------------------------------------------------------------------------------
--	tut06_flexible_domain_distribution.lua
--
--	We already have seen how a grid can be distributed using DistributeDomain(...).
--	In this tutorial however we'll take a look at how distribution can be controled
--	in a more flexible way.
--	We will thus take a look at how to prepare a PartitionMap and how to use
--	it in a call to RedistributeDomain.
--------------------------------------------------------------------------------

-- Since ug supports a bunch of different algebra modules we will choose one here.
-- This should always be the first thing you do in an ug-script.
-- The cpu-algebra is fine for now.
InitAlgebra(CPUAlgebraSelector())

-- include the basic util-methods.
ug_load_script("../ug_util.lua")

-- Get the command line parameters
dim = util.GetParamNumber("-dim", 2)
gridName = util.GetParam("-grid", "unit_square_quads_8x8.ugx")
outFileNamePrefix = util.GetParam("-o", "distributed_domain_")

-- We will save the created hierarchy to this file (with appended process id)
outHierarchyFilePrefix = util.GetParam("-oh", "hierarchy_on_proc_")

-- Create the domain and load a grid
dom = util.CreateDomain(dim)

if util.LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)



--------------------------------------------------------------------------------
--	New code starts here

--	This demonstration requires 16 processes
if GetNumProcesses() < 16 then
	print("At least 16 processes are required for this script. Aborting.")
	exit()
end

--	First we'll create a partition map. This is done on all processes, even
--	though the grid is currently only on process 0.
partitionMap = PartitionMap()

--	We'll also store the rank of the current process
procRank = GetProcessRank()

--	first we'll distribute sections of the domain to processes 0, 4, 8 and 12
firstTargets = {0, 4, 8, 12}
 
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
	PartitionDomain_RegularGrid(dom, partitionMap, 2, 2)
	
	-- We'll save the partition map. This should only be done for debugging.
	SaveGrid(dom:get_grid(), partitionMap:get_partition_handler(),
		"partitionMap_1_p" .. procRank .. ".ugx")
end


--	The partition map is set up. We can now distribute the domain.
--	RedistributeDomain takes three parameters: The domain, the
--	partitionMap, and a boolean indicating whether vertical interfaces
--	shall be created.
if RedistributeDomain(dom, partitionMap, false) == false then
	print("First redistribution failed. Please check your partitionMap.")
	exit()
end

--SaveDomain(dom, "dom_after_first_"..procRank..".ugx")
--exit()

--	The different target processes now all have their part of the grid.
--	We'll further distribute it using regular-grid-partitioning again.
--	Process 0 will distribute its grid among processes 0, 1, 2, 3.
--	Process 4 will distribute its grid among processes 4, 5, 6, 7.
--	...

partitionMap:clear()
for i, src in ipairs(firstTargets) do
	if src == procRank then
		-- Specify the first process and the number of processes.
		-- This will add processes src, src+1, ..., src+n-1
		partitionMap:add_target_procs(src, 4)
		
		-- We again use regular-grid partitioning 
		PartitionDomain_RegularGrid(dom, partitionMap, 2, 2)
		
		-- again we'll save the partition map. This should only be done for debugging.
		SaveGrid(dom:get_grid(), partitionMap:get_partition_handler(),
			"partitionMap_2_p" .. procRank .. ".ugx")
	end
end

--	The partition map is set up. We can now redistribute the domain
if RedistributeDomain(dom, partitionMap, false) == false then
	print("Second redistribution failed. Please check your partitionMap.")
	exit()
end

--	New code ends here
--------------------------------------------------------------------------------




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
-- We can access the grid of a domain through its get_grid() member method.
--
-- SaveGridHierarchy outputs a grid, where each level is assigned to a subset.
-- Original subsets are not contained in that file.
outFileName = outHierarchyFilePrefix .. GetProcessRank() .. ".ugx"
if SaveGridHierarchy(dom:get_grid(), outFileName) == false then
	print("Saving of grid-hierarch to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Again everything went fine.
print("Saved hierarchy to " .. outFileName)

-- we're done.
print("")
print("done")
