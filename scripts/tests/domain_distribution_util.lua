-- created by Sebastian Reiter
-- s.b.reiter@googlemail.com
-- Jan 2012

-- Create ddu namespace
ddu = ddu or {}

ug_load_script("ug_util.lua")

-- Defines methods which help in grid distribution.
-- This involves pre-refinement, post-refinement and hierarchical redistribution.
-- Methods exist, which automatically parse required parameters. Optionally 
-- parameters can be set from the outside.
--
-- All functions and variables reside in the ddu namespace.

-- This table holds all our parameters. It should not be altered by other scripts.
ddu.numPreRefs = 0
ddu.numRefs = 0
-- the method with which the domain / the grid will be distributed to the processes:
ddu.distributionType = "bisect" -- [grid2d | bisect | metis]
-- number of processes per node (only used if distType == grid2d)
-- should be a square number
ddu.numProcsPerNode = 1

-- number of processes to which the grid will be distributed during the initial
-- distribution step. This value will be recalculated when either
-- ParseAndInitializeParameters or SetAndInitializeParameters
-- is called
ddu.numInitialDistProcs = GetNumProcesses()


-- HIERARCHICAL REDISTRIBUTION
-- Note that hierarchical redistribution is not compatible with
-- distributionType == "grid2d"
-- If this is set to a level > numPreRefs and <= numRefs then hierarchical
-- redistribution is enabled. If it is set to -1, it is disabled (default).
-- Defines the level in which hierarchical redistribution starts.
-- Note that hierarchical redistribution is only performed, if a refinement follows.
ddu.hRedistFirstLevel = -1

-- The number of processes to which a grid is redistributed in a
-- hierarchical redistribution step.
ddu.hRedistNewProcsPerStep = -1

-- defines the number of refinements between hierarchical redistributions
-- if < 1 the hierarchical redistribution is disabled.
-- After s steps, the grid will be redistributed to
-- hRedistNewProcsPerStep processes.
-- Only has effect, if hierarchicalRedistFirstLevel ~= -1.
ddu.hRedistStepSize = 1

-- tells whehter horizontal redistribution is enabled or not (will be set during
-- parameter initialization)
ddu.hRedistEnabled = false

-- only if verbose is set to true, ddu will print something during distribution
ddu.verbosity = false



function ddu.print(msg)
	if ddu.verbosity == true then
		print(msg)
	end
end


function ddu.write(msg)
	if ddu.verbosity == true then
		write(msg)
	end
end


-- Parses default parameters
-- This method takes one parameter: dim (space-dimension of the problem)
-- returns true, if everything went well and false if something went wrong
function ddu.ParseAndInitializeParameters(dim)
	return ddu.SetAndInitializeParameters(
					util.GetParamNumber("-numPreRefs", 1),
					util.GetParamNumber("-numRefs",    3),
					util.GetParam("-distType", "bisect"),
					util.GetParamNumber("-numPPN", 1),
					util.GetParamNumber("-hRedistFirstLevel", -1),
					util.GetParamNumber("-hRedistNewProcsPerStep", math.pow(2, dim)),
					util.GetParamNumber("-hRedistStepSize", 1))
end


-- call this method if you want to specify all parameters manually
-- returns true, if everything went well and false if something went wrong
function ddu.SetAndInitializeParameters(numPreRefs, numRefs, distType,
				numProcsPerNode, hRedistFirstLevel, hRedistNewProcsPerStep,
				hRedistStepSize)
	-- set parameters
	ddu.numPreRefs = numPreRefs
	ddu.numRefs    = numRefs
	
	if ddu.numPreRefs > ddu.numRefs then
		print("It must be choosen: numPreRefs <= numRefs");
		return false
	end

	ddu.distributionType = distType
	ddu.numProcsPerNode = numProcsPerNode
	
	ddu.hRedistFirstLevel = hRedistFirstLevel
	ddu.hRedistNewProcsPerStep = hRedistNewProcsPerStep
	ddu.hRedistStepSize = hRedistStepSize
	
	-- the number of processes to which we will distribute the grid during
	-- initial distribution has to be calculated now.
	-- It depends on whether we perform hierarchical redistribution or not.
	if ddu.hRedistFirstLevel ~= -1 then
	--	make sure that all parameters are valid
		if ddu.hRedistFirstLevel <= ddu.numPreRefs or ddu.hRedistFirstLevel >= ddu.numRefs then
			print("HORIZONTAL REDISTRIBUTION ERROR:")
			print("  Make sure that numPreRefs < hRedistFirstLevel < numRefs")
			return false
		end
	
	--	make sure that distribution type grid2d is not used (the util-script is
	--	not suited for hierarchical redistribution)
		if ddu.distributionType == "grid2d" then
			print("HORIZONTAL REDISTRIBUTION ERROR:")
			print("  distributionType 'grid2d' is not supported for hierarchical"
				  .. " redistribution.")
			return false
		end
		
		ddu.hRedistEnabled = true
		
		local procs = GetNumProcesses()
		local refinements = numRefs - ddu.hRedistFirstLevel
		while refinements > 0 do
			refinements = refinements - ddu.hRedistStepSize
			if procs / ddu.hRedistNewProcsPerStep < 1 then
				break
			else
				procs = math.floor(procs / ddu.hRedistNewProcsPerStep)
			end
		end
		
		ddu.numInitialDistProcs = procs
	elseif hierarchicalRedistFirstLevel ~= -1 then
		print("WARNING: hierarchical redistribution disabled. To enable it, "
			  .. "set hierarchicalRedistFirstLevel to a value > numPreRefs and <= numRefs.")
	end
	
	return true
end


-- Prints all internal parameters
-- Specify a prefix, which is prepended to each line
function ddu.PrintParameters(prefix)
	print(prefix .. "numRefs    = " .. ddu.numRefs)
	print(prefix .. "numPreRefs = " .. ddu.numPreRefs)
	
	print(prefix .. "distType   = " .. ddu.distributionType)
	print(prefix .. "numPPN (numProcsPerNode) = " .. ddu.numProcsPerNode)
	print(prefix .. "hRedistFirstLevel = " .. ddu.hRedistFirstLevel)
	print(prefix .. "hRedistNewProcsPerStep = " .. ddu.hRedistNewProcsPerStep)
	print(prefix .. "hRedistStepSize = " .. ddu.hRedistStepSize)
end


-- Prints information for an analyzer script
function ddu.PrintAnalyzerInfo()
	if ddu.hRedistEnabled == true then
		print("#ANALYZER INFO: hierarchical redistribution: initial distribution at level "
				.. ddu.numPreRefs .. " to " .. ddu.numInitialDistProcs .. " processes.")
	end
	print("#ANALYZER INFO: NumProcs is " .. GetNumProcesses() .. ", numPreRefs = " .. ddu.numPreRefs
		  .. ", numRefs = " .. ddu.numRefs)
end


function ddu.RefineAndDistributeDomain(dom, verbosity)
	ddu.verbosity = verbosity
	-- Create a refiner instance.
	local refiner = GlobalDomainRefiner(dom)
	
	-- Performing pre-refines
	print("Performing (non parallel) pre-refinements")
	for i=1,ddu.numPreRefs do
		write( "PreRefinement step " .. i .. " ...")
		refiner:refine()
		print( " done.")
	end

	local numProcs = GetNumProcesses()
	local numDistProcs = ddu.numInitialDistProcs
	local numProcsWithGrid = 1
	local numCurRefs = ddu.numPreRefs

	-- Distribute the domain to all involved processes
	-- Only processes which already have a grid will fill their partition maps.
	local partitionMap = PartitionMap()

	while numDistProcs > 0 do
		partitionMap:clear()
		if GetProcessRank() < numProcsWithGrid then
		--	add target procs. Make sure to keep a portion on the local process
			partitionMap:add_target_proc(GetProcessRank())
			if numDistProcs > 1 then
				local firstNewProc = numProcsWithGrid + GetProcessRank() * (numDistProcs - 1)
				partitionMap:add_target_procs(firstNewProc, numDistProcs - 1)
			end
		
			if ddu.distributionType == "bisect" then
				PartitionDomain_Bisection(dom, partitionMap, 0)
		
			elseif ddu.distributionType == "grid2d" then
			--	Note that grid2d can not be used if hierarchical redistribution is active.
				local numNodesX, numNodesY = util.FactorizeInPowersOfTwo(numDistProcs / ddu.numProcsPerNode)
				util.PartitionMapLexicographic2D(dom, partitionMap, numNodesX,
												 numNodesY, ddu.numProcsPerNode)
		
			elseif ddu.distributionType == "metis" then
				print("metis is temporarily unsupported.")
				exit()
				--PartitionDomain_MetisKWay(dom, partitionMap, numProcs, baseLevel, 1, 1)
		
			else
			    print( "distributionType not known, aborting!")
			    exit()
			end
	
		-- save the partition map for debug purposes
			if verbosity >= 1 then
				print("saving partition map to 'partitionMap_p" .. GetProcessRank() .. ".ugx'")
				SavePartitionMap(partitionMap, dom, "partitionMap_p" .. GetProcessRank() .. ".ugx")
			end
		end
		
		print("Redistribute domain with 'distributionType' = '" .. ddu.distributionType .. "' ...")
		if RedistributeDomain(dom, partitionMap, true) == false then
			print("Redistribution failed. Please check your partitionMap.")
			exit()
		end
		print("  ... domain distributed.")
	
		numProcsWithGrid = numProcsWithGrid * numDistProcs
		numDistProcs = 0
		
		-- check whether we have to perform another hierarchical distribution.
		-- calculate the number of required refinements in this step.
		local maxRefsInThisStep = ddu.numRefs
		if ddu.hRedistEnabled == true then
			if numCurRefs < ddu.hRedistFirstLevel then
				maxRefsInThisStep = ddu.hRedistFirstLevel
			else
				maxRefsInThisStep = numCurRefs + ddu.hRedistStepSize
			end
			if maxRefsInThisStep >= ddu.numRefs then
				numDistProcs = 0
				maxRefsInThisStep = ddu.numRefs
			else
				numDistProcs = ddu.hRedistNewProcsPerStep
				if numProcsWithGrid * numDistProcs > numProcs then
					numDistProcs = 0
					maxRefsInThisStep = ddu.numRefs
				end
			end
		end
		
		-- Perform post-refine
		print("Refine Parallel Grid")
		for i = numCurRefs + 1, maxRefsInThisStep do
			write( "Refinement step " .. i .. " ...")
			refiner:refine()
			print( " done.")
		end
		numCurRefs = maxRefsInThisStep
	end
	
	-- clean up
	delete(partitionMap)
	delete(refiner)
end
