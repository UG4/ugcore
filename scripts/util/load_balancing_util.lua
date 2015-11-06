-- Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
-- Author: Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


--[[!
\file scripts/util/load_balancing_util.lua
\defgroup scripts_util_loadbalancing Load Balancing Utility
\ingroup scripts_util
\brief Parses parameters related to load-balancing and creates a load-balancer accordingly
\author Sebastian Reiter

All variables and functions in this script are contained in the namespace 'balancer'.
A set of global parameters with initial values is defined together with the methods
balancer.ParseParameters(), which is used to overwrite those initial values with
user-specified parameters, and balancer.CreateLoadBalancer(domain) which creates
a load-balancer for the specified domain according to the global parameters.

If you want to overwrite the default parameter values in your script, make sure to
assign them before calling balancer.ParseParameters but after loading this utility
file.
\{
]]--


ug_load_script("ug_util.lua")

balancer = balancer or {}

balancer.maxLvl 		= 100
balancer.firstDistLvl	= -1
balancer.firstDistProcs	= 256
balancer.redistSteps	= 2
balancer.redistProcs	= 256
balancer.maxDistLvl		= balancer.maxLvl
balancer.maxLvlsWithoutRedist = 4

balancer.parallelElementThreshold	= 32
balancer.qualityThreshold			= 0.8
balancer.balanceWeights				= nil
balancer.communicationWeights		= nil
balancer.childWeight				= 2
balancer.siblingWeight				= 2
balancer.itrFactor					= 1000
balancer.imbalanceFactor			= 1.05
balancer.noSiblingClustering 		= false

balancer.staticProcHierarchy = false

balancer.partitioner		= "bisection"

balancer.qualityRecordName	= "def-rebal"

balancer.parametersParsed	= false
balancer.defaultBalancer	= nil
balancer.defaultPartitioner = nil


--! Parses user-specified parameters related to load-balancing.
function balancer.ParseParameters()
	balancer.maxLvl			= util.GetParamNumber("-maxLvl", balancer.maxLvl,
								"The maximum level allowed. Elements in this level won't be refined")
	balancer.firstDistLvl	= util.GetParamNumber("-firstDistLvl", balancer.firstDistLvl,
								"The first level in which distribution will take place. Set it to -1 if you only want to use 'redistSteps'")
	balancer.firstDistProcs	= util.GetParamNumber("-firstDistProcs", balancer.firstDistProcs,
								"The number of processes onto which the domain will be distributed in firstDistLvl. Ignored if 'firstDistLvl' == -1.")
	balancer.redistSteps	= util.GetParamNumber("-redistSteps", balancer.redistSteps,
								"All levels with index = firstDistLvl + i * redistSteps, i = 0,1,...,n are considered redistribution levels")
	balancer.redistProcs	= util.GetParamNumber("-redistProcs", balancer.redistProcs,
								"The number of processes to which each process distributes grid-parts on each distribution level, i>0")
	balancer.maxDistLvl		= util.GetParamNumber("-maxDistLvl", balancer.maxDistLvl,
								"The maximum distribtion level. elements in levels above won't be partitioned and distributed separately.")

	balancer.maxLvlsWithoutRedist = util.GetParamNumber("-maxLvlsWithoutRedist", balancer.maxLvlsWithoutRedist)
	
	balancer.parallelElementThreshold	= util.GetParamNumber("-parallelElementThreshold", balancer.parallelElementThreshold,
											"No distribution is performed on a given level until each process can potentially receive 'parallelElementThreshold' elements.")
	balancer.qualityThreshold	= util.GetParamNumber("-qualityThreshold", balancer.qualityThreshold,
									"If the distribution quality lies below this threshold, redistribution is triggered")
	balancer.childWeight		= util.GetParamNumber("-childWeight", balancer.childWeight,
									"Children of elements in distribution levels are weighted with this number. integer number >= 1")
	balancer.siblingWeight		= util.GetParamNumber("-siblingWeight", balancer.siblingWeight,
									"Siblings always have to be located on 1 proc. This is considered during redistribution with this weight. integer number >= 1")
	balancer.itrFactor			= util.GetParamNumber("-itrFactor", balancer.itrFactor,
									"Weights the cost of communication versus redistribution. "..
									"Values in the range from 0.000001 to 1000000. A low value means that "..
									"communication time is considered low compared to redistribution time while "..
									"a high value means the contrary. Default is 1000.")
	balancer.imbalanceFactor	= util.GetParamNumber("-imbalanceFactor", balancer.imbalanceFactor, "This factor is set as allowed imbalance for each constraint."..
													  " From METIS manual: "..
													  "'This is an array of size ncon that specifies the allowed load imbalance tolerance "..
	 												  "for each constraint. For the ith partition and jth constraint the allowed weight is "..
	 												  "the ubvec[j]*tpwgts[i*ncon+j] fraction of the jth’s constraint total weight. "..
	 												  "The load imbalances must be greater than 1.0.'")

	balancer.noSiblingClustering = util.HasParamOption("-noSiblingClustering", "siblings should always be clustered if adaptive refinement is performed. Better distribution qualities if disabled.")

	balancer.staticProcHierarchy = balancer.staticProcHierarchy or util.HasParamOption("-staticProcHierarchy")
	
	balancer.partitioner		= util.GetParam("-partitioner", balancer.partitioner,
									"Options: parmetis, bisection, dynBisection. The partitioner which will be used during repartitioning.")
									
	balancer.parametersParsed = true
end

--! Prints the balancing parameters
function balancer.PrintParameters()
	print(" load balancing parameters:")
	print("    maxLvl                   = " .. balancer.maxLvl)
	print("    firstDistLvl             = " .. balancer.firstDistLvl)
	print("    firstDistProcs           = " .. balancer.firstDistProcs)
	print("    redistSteps              = " .. balancer.redistSteps)
	print("    redistProcs              = " .. balancer.redistProcs)
	print("    maxDistLvl               = " .. balancer.maxDistLvl)
	print("    maxLvlsWithoutRedist     = " .. balancer.maxLvlsWithoutRedist)
	print("    parallelElementThreshold = " .. balancer.parallelElementThreshold)
	print("    qualityThreshold         = " .. balancer.qualityThreshold)
	print("    childWeight              = " .. balancer.childWeight)
	print("    siblingWeight            = " .. balancer.siblingWeight)
	print("    itrFactor                = " .. balancer.itrFactor)
	print("    imbalanceFactor          = " .. balancer.imbalanceFactor)
	if balancer.regardAllChildren == true then
		print("    regardAllChildren          active")
	else
		print("    regardAllChildren          inactive")
	end
	if balancer.staticProcHierarchy == true then
		print("    staticProcHierarchy        active")
	else
		print("    staticProcHierarchy        inactive")
	end
	print("    partitioner              = " .. balancer.partitioner)
end


--! Creates a load-balancer for the given domain. The returned class is of the type
--! DomainLoadBalancer. The current set of global parameters in the balancer namespace
--! steers the creation.
--! The returned balancer can then be used to perform dynamic load-balancing and
--! rebalancing of adaptively refined grids.
function balancer.CreateLoadBalancer(domain)
	local loadBalancer = nil
	local numComputeProcs = NumProcs()
	
	if numComputeProcs > 1 then
		loadBalancer = DomainLoadBalancer(domain)
		local partitioner = nil
		if(balancer.partitioner == "parmetis") then
			RequiredPlugins({"Parmetis"})
			print("Creating ParmetisPartitioner")
			partitioner = Partitioner_Parmetis(domain)
			if balancer.balanceWeights ~= nil then
				partitioner:set_balance_weights(balancer.balanceWeights)
			end
			if balancer.communicationWeights ~= nil then
				partitioner:set_communication_weights(balancer.communicationWeights)
			end
			partitioner:set_child_weight(balancer.childWeight)
			partitioner:set_sibling_weight(balancer.siblingWeight)
			partitioner:set_itr_factor(balancer.itrFactor)
			partitioner:set_allowed_imbalance_factor(balancer.imbalanceFactor)
			partitioner:set_verbose(false)
		elseif(balancer.partitioner == "bisection") then
			partitioner = Partitioner_DynamicBisection(domain)
			partitioner:set_verbose(false)
			partitioner:enable_static_partitioning(true)
		elseif(balancer.partitioner == "dynBisection") then
			partitioner = Partitioner_DynamicBisection(domain)
			partitioner:set_verbose(false)
		else
			print("ERROR: Unknown partitioner specified in balancer.CreateLoadBalancer")
			exit()
		end

		partitioner:enable_clustered_siblings(not balancer.noSiblingClustering)
		
		balancer.defaultPartitioner = partitioner
		loadBalancer:set_partitioner(partitioner)
		loadBalancer:set_balance_threshold(balancer.qualityThreshold)
		loadBalancer:set_element_threshold(balancer.parallelElementThreshold)
		
	--	todo: use a more sophisticated algorithm to add distribution levels
		local processHierarchy = ProcessHierarchy()
		local lvl = balancer.redistSteps
		local procsTotal = 1
		if balancer.firstDistLvl ~= -1 then
			local numNew = balancer.firstDistProcs
			if(numNew > numComputeProcs) then
				numNew = numComputeProcs
			end
			processHierarchy:add_hierarchy_level(balancer.firstDistLvl, numNew)
			lvl = balancer.firstDistLvl + balancer.redistSteps
			procsTotal = numNew
		end
		
		if balancer.redistSteps > 0 then
			local numNewProcs = balancer.redistProcs
			--while procsTotal < numComputeProcs do
			for i = 1, balancer.maxLvl do
				if lvl > balancer.maxDistLvl then
					break
				end
				if procsTotal * numNewProcs <= numComputeProcs then
					processHierarchy:add_hierarchy_level(lvl, numNewProcs)
					procsTotal = procsTotal * numNewProcs
				elseif procsTotal < numComputeProcs then
					local numNew = math.floor(numComputeProcs / procsTotal)
					if(numNew > 0) then
						processHierarchy:add_hierarchy_level(lvl, numNew)
						procsTotal = procsTotal * numNew
					end
				else
					--processHierarchy:add_hierarchy_level(lvl, 1)
					break;
				end
				lvl = lvl + balancer.redistSteps
			end
		elseif balancer.firstDistLvl < 0 then
			processHierarchy:add_hierarchy_level(0, numComputeProcs)
			--loadBalancer:rebalance()
		end
		
		loadBalancer:set_next_process_hierarchy(processHierarchy);
	end
	
	if balancer.defaultBalancer == nil then
		balancer.defaultBalancer = loadBalancer
	end
	return loadBalancer
end


--!	If no loadBalancer is specified, this method performs rebalancing using
--! balancer.defaultBalancer. If balancer.defaultBalancer doesn't exist,
--! it will be automatically created using the current balancer parameters.
--! Note that changes to the balancer parameters won't have any effect
--! if balancer.defaultBalancer already exist.
--! However, by setting balancer.defaultBalancer to nil, you trigger
--! the recreation of the defaultBalancer with the new balancer parameters
--! in the next call to balancer.Rebalance.
function balancer.Rebalance(domain, loadBalancer)
	if loadBalancer == nil then
		if balancer.defaultBalancer == nil then
		--	create a default load balancer
			if balancer.parametersParsed == false then
				balancer.ParseParameters()
			end
			balancer.defaultBalancer = balancer.CreateLoadBalancer(domain)
		end
		loadBalancer = balancer.defaultBalancer
	end
	
	if loadBalancer ~= nil then
		if(balancer.staticProcHierarchy == false) then
			local minDistLvl = 0
			if balancer.firstDistLvl > 0 then minDistLvl = balancer.firstDistLvl end
			
			procH = CreateProcessHierarchy(domain, balancer.parallelElementThreshold,
										   balancer.redistProcs, NumProcs(),
										   minDistLvl, balancer.maxLvlsWithoutRedist)
			loadBalancer:set_next_process_hierarchy(procH)
		end
		loadBalancer:rebalance()
		loadBalancer:create_quality_record(balancer.qualityRecordName)
	end
end

--!	If no loadBalancer is specified, this method performs rebalancing using
--! balancer.defaultBalancer. See balancer.Rebalance for more information.
--! The method performs global refinement. After each refinement, rebalancing
--! is triggered.
function balancer.RefineAndRebalanceDomain(domain, numRefs, loadBalancer)
	balancer.Rebalance(domain, loadBalancer)
	if numRefs > 0 then
		refiner = GlobalDomainRefiner(domain)
		
		for i = 1, numRefs do
			TerminateAbortedRun()
			refiner:refine()
			TerminateAbortedRun()
			balancer.Rebalance(domain, loadBalancer)
		end
		
		delete(refiner)
	end
end

--[[!
\}
]]--
