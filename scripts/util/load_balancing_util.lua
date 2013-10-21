-- load balancing util
-- created by Sebastian Reiter
-- August 2013
-- s.b.reiter@gmail.com

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

balancer.parallelElementThreshold = 8
balancer.qualityThreshold	= 0.8
balancer.childWeight		= 2
balancer.siblingWeight		= 2
balancer.itrFactor			= 1000

balancer.partitioner		= "parmetis"


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
								"The number of processes to which each process distributes grid-parts on each distribution level (i>0)")
	balancer.maxDistLvl		= util.GetParamNumber("-maxDistLvl", balancer.maxDistLvl,
								"The maximum distribtion level. elements in levels above won't be partitioned and distributed separately.")

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

	balancer.partitioner		= util.GetParam("-partitioner", balancer.partitioner,
									"(Options: parmetis, bisection) The partitioner which will be used during repartitioning.")
end

--! Prints the balancing parameters
function balancer.PrintParameters()
	print(" load balancing parameters:")
	print("    maxLvl                   = " .. balancer.maxLvl)
	print("    firstDistLvl             = " .. balancer.firstDistLvl)
	print("    firstDistProcs           = " .. balancer.firstDistProcs)
	print("    redistSteps              = " .. balancer.redistSteps)
	print("    redistProcs              = " .. balancer.redistProcs)
	print("    parallelElementThreshold = " .. balancer.parallelElementThreshold)
	print("    qualityThreshold         = " .. balancer.qualityThreshold)
	print("    childWeight              = " .. balancer.childWeight)
	print("    siblingWeight            = " .. balancer.siblingWeight)
	print("    itrFactor                = " .. balancer.itrFactor)
	if balancer.regardAllChildren == true then
		print("    regardAllChildren          active")
	else
		print("    regardAllChildren          inactive")
	end
	print("    partitioner              = " .. balancer.partitioner)
end


--! Creates a load-balancer for the given domain. The returned class is of the type
--! DomainLoadBalancer. The current set of global parameters in the balancer namespace
--! steers the creation.
function balancer.CreateLoadBalancer(domain)
	local loadBalancer = nil
	local numComputeProcs = GetNumProcesses()
	
	if numComputeProcs > 1 then
		loadBalancer = DomainLoadBalancer(domain)
		
		if(balancer.partitioner == "parmetis") then
			if(ParmetisIsAvailable() == true) then
				print("Creating ParmetisPartitioner")
				local partitioner = Partitioner_Parmetis()
				partitioner:set_child_weight(balancer.childWeight)
				partitioner:set_sibling_weight(balancer.siblingWeight)
				partitioner:set_itr_factor(balancer.itrFactor)
				partitioner:set_verbose(false)
				loadBalancer:set_partitioner(partitioner)
			else
				print("ERROR: partitioner 'parmetis' specified in balancer.CreateLoadBalancer but ParMETIS isn't available.")
				exit()
			end
		elseif(balancer.partitioner == "bisection") then
			print("SORRY: partitioner 'bisection' currently isn't available due to lazyness...")
			exit()
		else
			print("ERROR: Unknown partitioner specified in balancer.CreateLoadBalancer")
			exit()
		end
		
		loadBalancer:set_balance_threshold(balancer.qualityThreshold)
		loadBalancer:set_element_threshold(balancer.parallelElementThreshold)
		
	--	todo: use a more sophisticated algorithm to add distribution levels
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

		local processHierarchy = ProcessHierarchy()
		
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
					processHierarchy:add_hierarchy_level(lvl, 1)
				end
				lvl = lvl + balancer.redistSteps
			end
		elseif balancer.firstDistLvl < 0 then
			processHierarchy:add_hierarchy_level(0, numComputeProcs)
			loadBalancer:rebalance()
		end
		
		loadBalancer:set_process_hierarchy(processHierarchy);
	end
	
	return loadBalancer
end

--[[!
\}
]]--
