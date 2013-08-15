-- load balancing util
-- created by Sebastian Reiter
-- s.b.reiter@gmail.com

ug_load_script("ug_util.lua")

balancer = balancer or {}

balancer.maxLvl 		= 100
balancer.firstDistLvl	= -1
balancer.firstDistProcs	= 2
balancer.redistSteps	= 2
balancer.redistProcs	= 2

balancer.parallelElementThreshold = 16
balancer.qualityThreshold	= 0.8
balancer.childWeight		= 2
balancer.siblingWeight		= 2
balancer.regardAllChildren	= false

balancer.partitioner		= "parmetis"


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


	balancer.parallelElementThreshold	= util.GetParamNumber("-parallelElementThreshold", balancer.parallelElementThreshold,
											"No distribution is performed on a given level until each process can potentially receive 'parallelElementThreshold' elements.")
	balancer.qualityThreshold	= util.GetParamNumber("-qualityThreshold", balancer.qualityThreshold,
									"If the distribution quality lies below this threshold, redistribution is triggered")
	balancer.childWeight		= util.GetParamNumber("-childWeight", balancer.childWeight,
									"Children of elements in distribution levels are weighted with this number. integer number >= 1")
	balancer.siblingWeight		= util.GetParamNumber("-siblingWeight", balancer.siblingWeight,
									"Siblings always have to be located on 1 proc. This is considered during redistribution with this weight. integer number >= 1")
	balancer.regardAllChildren	= util.HasParamOption("-regardAllChildren",
									"Regard all children during partitioning to calculate the weight of a parent")

	balancer.partitioner		= util.GetParam("-partitioner", balancer.partitioner,
									"(Options: parmetis, bisection) The partitioner which will be used during repartitioning.")
end


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
				partitioner:set_verbose(false)
				partitioner:set_regard_all_children(balancer.regardAllChildren)
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
		if balancer.firstDistLvl ~= -1 then
			loadBalancer:add_distribution_level(balancer.firstDistLvl,
												balancer.firstDistProcs)
			lvl = balancer.firstDistLvl + balancer.redistSteps
		end
		
		if balancer.redistSteps > 0 then
			local numNewProcs = balancer.redistProcs
			local procsTotal = 1
			--while procsTotal < numComputeProcs do
			for i = 1, balancer.maxLvl do
				procsTotal = procsTotal * numNewProcs
				if procsTotal <= numComputeProcs then
					loadBalancer:add_distribution_level(lvl, numNewProcs)
				else
					loadBalancer:add_distribution_level(lvl, 1)
				end
				lvl = lvl + balancer.redistSteps
			end
		elseif balancer.firstDistLvl < 0 then
			loadBalancer:add_distribution_level(0, numComputeProcs)
			loadBalancer:rebalance()
		end
	end
	
	return loadBalancer
end
