-- balancerDesc =
-- {
-- 	partitioner = "dynamicBisection",

-- 	hierarchy =
-- 	{
-- 		minElemsPerProcPerLevel		= 32,
-- 		maxRedistProcs				= 256,
-- 		qualityRedistLevelOffset	= 2,

-- 		{-- levels 0 to 1
-- 			upperLvl = 1,
-- 			maxProcs = 1
-- 		},

-- 		{-- levels 2 to 5
-- 			upperLvl = 5,
-- 			maxProcs = 64,
--			maxRedistProcs 8
-- 		},
-- 	}
-- }

-- procH = balancer.CreateProcessHierarchy(dom, balancerDesc.hierarchyDesc)
-- loadBalancer = balancer.CreateLoadBalancer(dom, balancerDesc)


-- IDEAS
--Create a dummy-load-balancer in a serial environment

util = util or {}
util.balancer = util.balancer or {}


util.balancer.defaults =
{
	partitioner = "dynamicBisection",
	hierarchy = "standard",
	qualityThreshold = 0.9,

	partitioners =
	{
		dynamicBisection =
		{
			verbose = false,
			enableXCuts = true,
			enableYCuts = true,
			enableZCuts = true,
			clusteredSiblings = true,
			balanceThreshold = 0.9
		},

		staticBisection =
		{
			verbose = false,
			enableXCuts = true,
			enableYCuts = true,
			enableZCuts = true,
			clusteredSiblings = true,
			balanceThreshold = 0.9
		},

		parmetis =
		{
			balanceWeights = nil,
			communicationWeights = nil,
			childWeight		= 2,
			siblingWeight	= 2,
			itrFactor		= 1000,
			verbose = false,
			clusteredSiblings = true,
			balanceThreshold = 0.9
		}
	},

	partitionPostProcessors =
	{
		smoothPartitionBounds = {
		},

		clusterElementStacks = {
			stackingDir = MakeVec(0, 0, 1)
		}
	},

	hierarchies =
	{
		standard = {
			minElemsPerProcPerLevel		= 32,
			maxRedistProcs				= 64,
			qualityRedistLevelOffset	= 2
		},

		lessRedists = {
			minElemsPerProcPerLevel		= 32,
			maxRedistProcs				= 256,
			qualityRedistLevelOffset	= 3
		}
	}
}


function util.balancer.CondAbort(condition, message)
	if condition == true then
		print("ERROR in util.balancer: " .. message)
		exit()
	end
end

function util.balancer.CreateLoadBalancer(dom, desc)
	if util.tableDesc.IsPreset(desc) then return desc end

	local defaults = util.balancer.defaults
	if desc == nil then desc = defaults end

	--todo:	support a serial load-balancer
	if NumProcs() == 1 then
		return {
			descriptor	= desc,
			rebalance =	function (recordName) end,
			create_quality_record =	function (name) end,
			print_quality_records =	function () print("  ---") end
		}
	end

	local loadBalancer = DomainLoadBalancer(dom)

	print("  util.balancer: creating partitioner...")
	local partitionerDesc = desc.partitioner or defaults.partitioner
	local part = util.balancer.CreatePartitioner(dom, partitionerDesc)

	if desc.partitionPostProcessor ~= nil then
		print("  util.balancer: creating partition post processor...")
		local postProc = util.balancer.CreatePostProcessor(dom, desc.partitionPostProcessor)
		part:set_partition_post_processor(postProc)
	end

	loadBalancer:set_partitioner(part)

	print("  util.balancer: creating process hierarchy...")
	local hierarchyDesc = desc.hierarchy or defaults.hierarchy
	local procH, elemThreshold = util.balancer.CreateProcessHierarchy(dom, hierarchyDesc)
	loadBalancer:set_next_process_hierarchy(procH)

	loadBalancer:set_balance_threshold(desc.qualityThreshold or defaults.qualityThreshold)
	loadBalancer:set_element_threshold(elemThreshold)

	if desc then desc.instance = loadBalancer end

	print("  util.balancer: done")
	return {
		balancer	= loadBalancer,
		descriptor	= desc,
		partitioner	= part,
		hierarchy	= procH,

		rebalance =	function (recordName)
						local hierarchyDesc = desc.hierarchy or util.balancer.defaults.hierarchy
						local ph = util.balancer.CreateProcessHierarchy(dom, hierarchyDesc)
						print("new prochess hierarchy:")
						print(ph:to_string())
						loadBalancer:set_next_process_hierarchy(ph)
						loadBalancer:rebalance()
						if recordName == nil then recordName = "def-rebal" end
						loadBalancer:create_quality_record(recordName)
					end,

		create_quality_record =	function (name)
									if name == nil then name = "def-rebal" end
									loadBalancer:create_quality_record(name)
								end,

		print_quality_records =	function ()
									loadBalancer:print_quality_records()
								end
	}
end


function util.balancer.CreatePartitioner(dom, partitionerDesc)
	if util.tableDesc.IsPreset(partitionerDesc) then return partitionerDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(partitionerDesc)
	local defaults = util.balancer.defaults.partitioners[name]
	if desc == nil then desc = defaults end

	local verbose = true;
	if desc.verbose ~= nil then verbose = desc.verbose
	elseif defaults.verbose ~= nil then verbose = defaults.verbose end

	local partitioner = nil

	if(name == "dynamicBisection") then
		partitioner = Partitioner_DynamicBisection(dom)
		partitioner:set_verbose(verbose)
		if desc.enableXCuts == false then partitioner:enable_split_axis(0, false) end
		if desc.enableYCuts == false then partitioner:enable_split_axis(1, false) end
		if desc.enableZCuts == false then partitioner:enable_split_axis(2, false) end
	elseif(name == "staticBisection") then
		partitioner = Partitioner_DynamicBisection(dom)
		partitioner:set_verbose(verbose)
		partitioner:enable_static_partitioning(true)
		if desc.enableXCuts == false then partitioner:enable_split_axis(0, false) end
		if desc.enableYCuts == false then partitioner:enable_split_axis(1, false) end
		if desc.enableZCuts == false then partitioner:enable_split_axis(2, false) end
	elseif(name == "parmetis") then
		RequiredPlugins({"Parmetis"})

		partitioner = Partitioner_Parmetis(dom)
		partitioner:set_child_weight(desc.childWeight or defaults.childWeight)
		partitioner:set_sibling_weight(desc.siblingWeight or defaults.siblingWeight)
		partitioner:set_itr_factor(desc.itrFactor or defaults.itrFactor)
		partitioner:set_verbose(verbose)

	--todo: create a method that creates balance- and communicationWeights using a descriptor
		if util.tableDesc.IsPreset(desc.balanceWeights) then
			partitioner:set_balance_weights(desc.balanceWeights)
		end
		if util.tableDesc.IsPreset(desc.communicationWeights) then
			partitioner:set_communication_weights(desc.communicationWeights)
		end
	else
		print("ERROR: Unknown partitioner specified in balancer.CreateLoadBalancer: " .. name)
		exit()
	end

	if desc.clusteredSiblings ~= nil then
		partitioner:enable_clustered_siblings(desc.clusteredSiblings)
	elseif defaults.clusteredSiblings ~= nil then
		partitioner:enable_clustered_siblings(defaults.clusteredSiblings)
	else
		partitioner:enable_clustered_siblings(true)
	end

	if desc then desc.instance = partitioner end

	return partitioner
end


function util.balancer.CreatePostProcessor(dom, postProcDesc)
	if util.tableDesc.IsPreset(postProcDesc) then return postProcDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(postProcDesc)
	local defaults = util.balancer.defaults.partitionPostProcessors[name]
	if desc == nil then desc = defaults end

	local postProc = nil
	if name == "smoothPartitionBounds" then
		postProc = SmoothPartitionBounds()
	elseif name == "clusterElementStacks" then
		postProc = ClusterElementStacks(dom:position_attachment(),
										desc.stackingDir or defaults.stackingDir)
	else
		print("ERROR: Unknown partitionPostProcessor specified in " ..
			  " util.balancer. CreatePostProcessor: " .. name)
		exit()
	end

	if desc then desc.instance = postProc end
	return postProc
end


function util.balancer.CreateProcessHierarchy(dom, hierarchyDesc)
	if util.tableDesc.IsPreset(hierarchyDesc) then return hierarchyDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(hierarchyDesc)
	local defaults = util.balancer.defaults.hierarchies[name]
	if desc == nil then desc = defaults end

	local procH = ProcessHierarchy()
	if desc then desc.instance = procH end

	local numComputeProcs = NumProcs()

	local elemThreshold = desc.minElemsPerProcPerLevel or defaults.minElemsPerProcPerLevel
	local qualityRedistLvlOffset = desc.qualityRedistLevelOffset or defaults.qualityRedistLevelOffset
	local defNumRedistProcs = desc.maxRedistProcs or defaults.maxRedistProcs or numComputeProcs
	local defMaxProcs = desc.maxProcs or defaults.maxProcs or numComputeProcs
	
	if numComputeProcs == 1 then
		procH:add_hierarchy_level(0, 1)
		return procH, elemThreshold
	end


	if qualityRedistLvlOffset < 2 then qualityRedistLvlOffset = 2 end


	local domInfo = dom:domain_info()
	if domInfo:num_levels() == 0 then
		procH:add_hierarchy_level(0, 1)
		return procH, elemThreshold
	end

	local largestLvl = 0
	for ilvl = 1, domInfo:num_levels() - 1 do
		if	  domInfo:num_elements_on_level(ilvl)
			> domInfo:num_elements_on_level(largestLvl)
		then
			largestLvl = ilvl
		end
	end

	local maxElemsPerLvl = domInfo:num_elements_on_level(largestLvl)
	if maxElemsPerLvl == 0 then
		procH:add_hierarchy_level(0, 1)
		return procH, elemThreshold
	end

	local defNumProcsInvolved = math.floor(maxElemsPerLvl / elemThreshold)
	if defNumProcsInvolved > numComputeProcs then
		defNumProcsInvolved = numComputeProcs
	elseif defNumProcsInvolved < 1 then
		procH:add_hierarchy_level(0, 1)
		return procH, elemThreshold
	end

	local curProcs = 1
	local lastDistLvl = -2

	for curLvl = 0, domInfo:num_levels() - 1 do
		if lastDistLvl + 1 < curLvl then
			local redistProcs = defNumRedistProcs
			local maxProcs = defMaxProcs
			local hints = nil

			-- check if a subsection exists which contains the current level
			for name, val in pairs(desc) do
				if type(val) == "table" then
					if val.upperLvl == nil or val.upperLvl >= curLvl then
					--	overwrite values from this section
						maxProcs = val.maxProcs or maxProcs
						redistProcs = val.maxRedistProcs or redistProcs
						if type(val.hints) == "table" then
							hints = val.hints
						end
						break
					end
				end
			end

			if maxProcs > numComputeProcs then maxProcs = numComputeProcs end
			if redistProcs > maxProcs then redistProcs = maxProcs end

			-- The following line only ensures that v-interfaces stay on fixed levels
			-- todo: remove the following line as soon as v-interfaces may change levels.
			numProcsInvolved = math.floor(defNumProcsInvolved / redistProcs) * redistProcs

			if maxProcs > numProcsInvolved then maxProcs = numProcsInvolved end

			if curProcs * redistProcs > maxProcs then
				redistProcs = math.floor(maxProcs / curProcs)
				if redistProcs <= 1 then
					if curProcs == numComputeProcs then
						break
					end
				end
			end

			if	redistProcs > 1
				and curProcs * redistProcs * elemThreshold
					<= domInfo:num_elements_on_level(curLvl)
			then
				lastDistLvl = curLvl
				local hlvl = procH:num_hierarchy_levels()
				procH:add_hierarchy_level(curLvl, redistProcs)
				if hints ~= nil then
					for name, val in pairs(hints) do
						procH:add_partition_hint(hlvl, name, Variant(val))
					end
				end
				curProcs = curProcs * redistProcs
			end
		end
	end

	if lastDistLvl < 0 then
		procH:add_hierarchy_level(0, 1)
		return procH, elemThreshold
	end

	if curProcs == numComputeProcs then
		local numQualityRedists = math.floor((domInfo:num_levels() - lastDistLvl) / qualityRedistLvlOffset);
		if numQualityRedists > 0 then
			for i = 1, numQualityRedists do
				local lvl = lastDistLvl + i * qualityRedistLvlOffset
				if lvl < domInfo:num_levels() then
					procH:add_hierarchy_level(lvl, 1);
				end
			end
		end
	end
	return procH, elemThreshold
end
