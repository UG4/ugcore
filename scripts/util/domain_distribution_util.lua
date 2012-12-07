-- created by Sebastian Reiter
-- s.b.reiter@googlemail.com

util = util or {}

--[[!
\file scripts/util/domain_distribution_util.lua
\brief	creates partition maps of different structure and provides some
		easy to use domain distribution methods.
]]--


--!	Distributes the top-level of a domain to the given number of processes.
--! This method has to be called by all processes with the same parameters!
--! @note Some paramters are optional. nil is a valid value for each optional parameter.
--! @return (bool) Returns whether the method was a success
--! @param dom	(Domain) A valid domain instance.
--! @param partitioningMethod	(optional string) Choose the partitioning method.
--!								Valid values are: "bisection", "metis"
--!								default is "bisection". "metis" is only available,
--!								if ug has been built with metis support.
--! @param verticalInterfaces	(optional bool) Trigger creation of vertical interfaces.
--!								vertical interfaces are required for multi-grid solvers.
--!								If you're only going to solve on the surface grid,
--!								vertical interfaces may introduce an unnecessary overhead.
--!								Default is true
--! @param numTargetProcs	(optional integer) The number of target processes to which the domain
--!							will be distributed. This shouldn't be more than
--!							there are elements in the top-level.
--!							Default is GetNumProcesses()
--! @param distributionLevel	(optional integer) Sets the level on which distribution
--!								is performed. Default is the domains top-level.
--!								Currently only supported in partitioningMethod "metis".
function util.DistributeDomain(dom, partitioningMethod, verticalInterfaces,
							   numTargetProcs, distributionLevel)
	local partitionMap = PartitionMap()
		
	if partitioningMethod == nil then
		partitioningMethod = "bisection"
	end
	
	if verticalInterfaces == nil then
		verticalInterfaces = true
	end
	
	if numTargetProcs == nil then
		numTargetProcs = GetNumProcesses()
	end
	
	if distributionLevel == nil then
		distributionLevel = dom:grid():num_levels() - 1
		if distributionLevel < 0 then
			distributionLevel = 0
		end
	end
	
	if partitioningMethod == "bisection" then
		if distributionLevel < dom:grid():num_levels() - 1 then
			print("WARNING in util.DistributeDomain: 'bisection' can currently "
				  .. "only be performed on the top level. Sorry...")
		end 
		util.PartitionMapBisection(dom, partitionMap, numTargetProcs)
	elseif partitioningMethod == "metis" then
		util.PartitionMapMetis(dom, partitionMap, numTargetProcs, distributionLevel)
	else
		print("ERROR in util.DistributeDomain: Unknown partitioning method.")
		print("  Valid partitioning methods are: 'bisection' and 'metis'")
		return
	end
	
	local success = DistributeDomain(dom, partitionMap, verticalInterfaces)
	
	delete(partitionMap)
	
	return success
end


--! create a partition map by performing repeated bisection
--! This method automatically adds target procs 0, ..., numProcs - 1.
--! It is thus not suited for hierarchical redistribution.
--! Use PartitionDomain_Bisection directly in this case.
function util.PartitionMapBisection(dom, partitionMapOut, numProcs)
	partitionMapOut:clear()
	partitionMapOut:add_target_procs(0, numProcs)
	PartitionDomain_Bisection(dom, partitionMapOut, 0)
end

--! create a partition map by using metis graph partitioning.
--! This only works if Metis is available in the current build.
--! This method automatically adds target procs 0, ..., numProcs - 1.
--! It is thus not suited for hierarchical redistribution.
--! Use PartitionDomain_MetisKWay directly in this case.
function util.PartitionMapMetis(dom, partitionMapOut, numProcs, baseLevel)
	partitionMapOut:clear()
	partitionMapOut:add_target_procs(0, numProcs)
	PartitionDomain_MetisKWay(dom, partitionMapOut, numProcs, baseLevel, 1, 1)
end


--! performs lexicographic ordering from the lower left to the upper right.
--! Since each node can consist of multiple processes, we can further
--! perform a 'hierarchical' lexicographic ordering.
--! Set numProcsPerNode to 1 to retrieve a default ordering.
--! Set it to any square of n (n = 2, 3, ...) to retrieve lexicographic ordered
--! subgrids in each node (a node is considered to be a square grid itself).
--! If numProcsPerNode is not a square-number, then it will be defaulted to 1.
--! Since this method adds target procs 0, ..., numNodesX * numNodesY * numProcsPerNode
--! to the given partitionMap, it is not suited for redistribution.
function util.PartitionMapLexicographic2D(dom, partitionMapOut, numNodesX,
										  numNodesY, numProcsPerNode)
	partitionMapOut:clear()
	local strWarning = "WARNING in util.PartitionMap_Lexicographic: "
	
--	get the subgrid width and make sure that numProcsPerNode is a square number.
--	Default them to 1, if numProcsPerNode is not a square number.
	local subGridWidth = math.floor(math.sqrt(numProcsPerNode))
	if subGridWidth * subGridWidth ~= numProcsPerNode then
		print(strWarning .. "numProcsPerNode has to be a square number. Defaulting to 1.")				
		numProcsPerNode = 1
		subGridWidth = 1
	end
	
--	Assign process-ids to the partitionMap. Consider the subgrids...
	local curProcX = 0
	
	for iy = 0, numNodesY - 1 do
		for subY = 0, subGridWidth - 1 do
			local curProcX = iy * numNodesX * numProcsPerNode + subY * subGridWidth
			for ix = 0, numNodesX - 1 do
				for subX = 0, subGridWidth - 1 do
					partitionMapOut:add_target_proc(curProcX + subX)
				end
				curProcX = curProcX + numProcsPerNode
			end
		end
	end
	
--	now perform the actual partitioning.
	local procsX = numNodesX * subGridWidth
	local procsY = numNodesY * subGridWidth
	
--	TODO: the 'surfaceOnly' parameter should be handled in a more flexible way.
	PartitionDomain_RegularGrid(dom, partitionMapOut, procsX, procsY, true)
end