----------------------------------------------------------
--
--   Lua - Script to generate lexicographic "right up" partition map,
--         without "second partitioning stage"
--   Author: Ingo Heppner, 05052011
--
----------------------------------------------------------

----------------------------------------------------------
-- ATTENTION - DEPRECIATED!
----------------------------------------------------------
-- This script served as base for "partition_maps.lua". Parts were transfered
-- to ug_util.lua. Please use "partition_maps.lua" instead of this script, since
-- it may no longer be updated.


--------------------------------------------------------------------------------
-- auxiliary functions - maybe to be added to 'ug_util.lua' sometimes
--------------------------------------------------------------------------------
-- function to factorise number which has to be a power of 2 in two factors
-- which differ by a factor of 2 and returns the larger one
function util.factoriseInPowersOfTwo(n)
	if not util.IsPowerOfTwo(n) then
		print("Number to factorise must be a power of 2. Aborting.")
		exit()
	end

	local number firstFactor = n
	local number secFactor = 1

	while (firstFactor > 2*secFactor) do
		firstFactor = firstFactor/2
		secFactor = secFactor*2
	end

	return firstFactor
end

-- Eigentlich in aufrufendem File - oder automatisch - zu setzen:
-- set number of cores per node
numCoresPerCNode = 4 -- on e.g. cekon, jugene; TODO: algorithm works only for this value in the moment!

outFileNamePrefix = util.GetParam("-o", "distributed_domain_") -- lohnt sich nicht wirklich, das als Parameter an die Funktion zu uebergeben ...

--------------------------------------------------------------------------------
-- New lexicographic partitioning in only one step
--------------------------------------------------------------------------------
function util.lexicographicMapping(numProcs, numCoresPerCNode, gridName, verbosity)

	if not util.isPowerOfTwo(numProcs) then
		print("Number of processes must be a power of 2. Aborting.")
		exit()
	end

	local number procRank    = -1
	local number numCNodes   = -1
	local number numCellCols = -1
	local number numCellRows = -1
	local number p0, p1, p2, p3 = 0, 1, 2, 3

	-- get number of compute nodes ("CNodes; this equals number of "cells")
	local numCNodes = math.ceil(numProcs/numCoresPerCNode)

	print("NumProcs is " .. numProcs .. ", NumCoresPerNode is " .. numCoresPerCNode .. ", numCNodes = " .. numCNodes .. " (= number of cells)")
	print("numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)


	-- number of cells in x direction
	local numCellCols = util.factoriseInPowersOfTwo(numCNodes) -- returns larger factor
	-- number of cells in y direction
	local numCellRows = numCNodes/numCellCols

	-- check
	if numProcs ~= (numCellRows*numCellCols*numCoresPerCNode) then
		print("Factorisation in cells failed. Aborting.")
		exit()
	end

	print("   Partitioning divides domain in '" .. numCellCols * numCellRows .. "' cells, organised in '" .. numCellCols .. "' columns and '" .. numCellRows .. "' rows")

	--	First we'll create a partition map. This is done on all processes, even
	--	though the grid is currently only on process 0.
	local partitionMap = PartitionMap()

	--	We'll also store the rank of the current process
	local procRank = GetProcessRank()

	--	Since only process 0 has a grid at the moment, we'll only fill the
	--	partition map on him.
	if procRank == 0 then

		-- The following call performs the load-balancing and partitions the
		-- domain into several sections, which will then be sent to the different
		-- processes.
		-- Note that several different load-balancers exist. In this example
		-- we'll simply use a regular-grid with cx x cy cells.

		-- add procs to map
		--numProcsInX = numCellCols*4

		local p0 = 0;
		local p1 = 1;
		local p2 = 2;
		local p3 = 3;

		for jy = 1, numCellRows, 1 do
			if verbosity >= 1 then
				print("cell row " .. jy .. ": ")
				print("  line 1 ")
			end
			for ix = 1, numCellCols, 1 do
				partitionMap:add_target_proc(p0)
				if verbosity >= 1 then
					print("   target proc '" .. p0 .. "' added!")
				end
				partitionMap:add_target_proc(p1)
				if verbosity >= 1 then
					print("   target proc '" .. p1 .. "' added!")
				end

				p0 = p0 + numCoresPerCNode
				p1 = p1 + numCoresPerCNode
			end

			if verbosity >= 1 then
				print("  line 2 ")
			end
			for ix = 1, numCellCols, 1 do
				partitionMap:add_target_proc(p2)
				if verbosity >= 1 then
					print("   target proc '" .. p2 .. "' added!")
				end
				partitionMap:add_target_proc(p3)
				if verbosity >= 1 then
					print("   target proc '" .. p3 .. "' added!")
				end
				p2 = p2 + numCoresPerCNode
				p3 = p3 + numCoresPerCNode
			end
		end

		-- For this algorithm it is crucial that the number of processes in the
		-- partitionMap match the number of cells specified for the load-balancer.
		-- TODO: Handling of 3. and 4. parameter is not yet general (05052011ih):
		PartitionDomain_RegularGrid(dom, partitionMap, numCellCols*numCoresPerCNode/2, numCellRows*numCoresPerCNode/2, true) -- <== 5. Parameter wieder 'true' seit 08h30, 05052011
	
		-- We'll save the partition map. This should only be done for debugging.
		if verbosity >= 1 then
			SaveGrid(dom:grid(), partitionMap:get_partition_handler(),
				"partitionMap_1_p" .. procRank .. ".ugx")
		end
	end -- end 'if procRank == 0'

	--	The partition map is set up. We can now distribute the domain.
	--	RedistributeDomain takes three parameters: The domain, the
	--	partitionMap, and a boolean indicating whether vertical interfaces
	--	shall be created.
	if RedistributeDomain(dom, partitionMap, true) == false then -- <== 5. Parameter sollte eigentlich wieder 'true' sein, seit 08h30, 05052011
		print("First redistribution failed. Please check your partitionMap.")
		exit()
	end

	if verbosity >= 1 then
		if TestDomainInterfaces(dom) == false then
			print("Inconsistent grid layouts after distribution. Aborting.")
			exit()
		end
	end

	-- Lets save the domain on each process
	if verbosity >= 1 then
		outFileName = outFileNamePrefix .. GetProcessRank() .. ".ugx"
		if SaveDomain(dom, outFileName) == false then
			print("Saving of domain to " .. outFileName .. " failed. Aborting.")
			exit()
		end
		-- Everything seems to went fine.
		print("Saved domain to " .. outFileName)
	end
end
