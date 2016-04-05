-- Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

ug_load_script("load_balancing_util_2.lua")

util = util or {}
util.refinement = util.refinement or {}

--! Creates a regular hierarchy through global refinement. Redistribution
--! is handled through 'load_balancing_util_2.lua'. To this end a balancerDesc
--! may optionally be specified, which is used to create the invloved balancer.
--!
--! @param dom			(Domain) the domain on which a hierarchy shall be created
--! @param numRefs		(integer) the number of refinements that shall be performed
--! @param verbose		(boolean)(optional)	if 'true', information on refinements and
--!							distribution-qualities are printed.
--! @param balancerDesc	(string or table)(optional) a descriptor according to
--!							"load_balancing_util_2.lua"
function util.refinement.CreateRegularHierarchy(dom, numRefs, verbose, balancerDesc)
	local ref = GlobalDomainRefiner(dom)
	local bal = util.balancer.CreateLoadBalancer(dom, balancerDesc)

	bal.rebalance()

	for i = 1, numRefs do
		if verbose then
			print("- refining level " .. i-1)
		end
		ref:refine()
		bal.rebalance("adaption-"..i)
	end

	if verbose and NumProcs() > 1 then
		print("\nDistribution quality statistics:")
		bal.print_quality_records()
	end
end