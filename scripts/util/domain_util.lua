-- Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

util = util or {}

--------------------------------------------------------------------------------
-- Subset utils
--------------------------------------------------------------------------------


--! util.CheckSubsets
--! checks if all required subsets are contained in the SubsetHandler
--! @param dom Domain
--! @param neededSubsets List of subsets the SubsetHandler must contain
--! @return true if all subsets are contained, false else
function util.CheckSubsets(dom, neededSubsets)
	sh = dom:subset_handler()
	for i, tval in ipairs(neededSubsets) do
		if sh:get_subset_index(tval) == -1 then
			print("Domain does not contain subset '"..tval.."'.")
			return false
		end
	end
	
	return true
end


--! Creates a new domain and loads the specified grid. The method then performs
--! numRefs global refinements.
--! A list of subset-names can be specified which have to be present in the loaded grid.
--! The method returns the created domain.
--! @note Some parameters are optional. nil is a valid value for each optional parameter.
--! @return	(Domain) the created domain
--! @param gridName	(string) The filename of the grid which shall be loaded.
--!					The grid is searched in a path relative to the current path
--!					first. If it isn't found there, the path is interpreted as
--!					an absolute path. If the grid still can't be found, the method
--!					tries to load it from UG_BASE/data/grids.
--! @param numRefs	(int) The total number of global refinements
--! @param numPreRefs	(int) The number of refinements that are performed before
--!						distribution.
--! @param neededSubsets	(optional, list of strings) The subsets that are required
--!							by the simulation. If not all those subsets are present,
--!							the method aborts. Default is an empty list.
--! @param noIntegrityCheck	(optional, bool) Disables integrity check if 'true'.
function util.CreateDomain(gridName, numRefs, neededSubsets, noIntegrityCheck)

	-- create Instance of a Domain
	local dom = Domain()
	
	-- load domain
	write("Loading Domain "..gridName.." ... ") 
	LoadDomain(dom, gridName)
	write("done.\n")

	if noIntegrityCheck ~= true then
		write("Performing integrity check on domain ... ")
		if CheckForUnconnectedSides(dom:grid()) == true then
			write("WARNING: unconnected sides found (see above).\n")
			local note = "NOTE: You may disable this check by passing 'true' "..
				  		 "to 'noIntegrityCheck' in 'util.CreateDomain'.\n"
			write(note)
			errlog(note)
		end
		write("done.\n")
	end

	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	if numRefs == nil then numRefs = 0 end
	if numRefs > 0 then
		write("Refining("..numRefs.."): ")
		local refiner = GlobalDomainRefiner(dom)
		for i=1,numRefs do
			TerminateAbortedRun()
			refiner:refine()
			write(i .. " ")
		end
		write("done.\n")
		delete(refiner)
	end
	
	-- check whether required subsets are present
	if neededSubsets ~= nil then
		ug_assert(util.CheckSubsets(dom, neededSubsets) == true, 
			"Something wrong with required subsets. Aborting.");
	end
	
	-- return the created domain
	return dom
end


--! Creates a new domain and loads the specified grid. The method then performs
--! numPreRefs refinements before it distributes the domain onto the available
--! processes. The partitioning method can be chosen through distributionMethod.
--! After distribution the domain is refined until a total of numRefs refinement
--! steps has been performed (including numPreRefs).
--! A list of subset-names can be specified. After distribution the methods checks
--! Whether all processes received the required subsets.
--! The method returns the created domain.
--! @note Some parameters are optional. nil is a valid value for each optional parameter.
--! @return	(Domain) the created domain
--! @param gridName	(string) The filename of the grid which shall be loaded.
--!					The grid is searched in a path relative to the current path
--!					first. If it isn't found there, the path is interpreted as
--!					an absolute path. If the grid still can't be found, the method
--!					tries to load it from UG_BASE/data/grids.
--! @param numRefs	(int) The total number of global refinements
--! @param numPreRefs	(int) The number of refinements that are performed before
--!						distribution.
--! @param neededSubsets	(optional, list of strings) The subsets that are required
--!							by the simulation. If not all those subsets are present,
--!							the method aborts. Default is an empty list.
--! @param distributionMethod	(optional, string) The distribution method.
--!								Either "bisection" or "metis". Default is "bisection".
--!								See util.DistributeDomain for more information
--!								(in UG_BASE/scripts/util/domain_distribution.lua)
--! @param verticalInterfaces	(optional, bool) Vertical interfaces are required
--!								by multi-grid solvers. Default is true.
--!								See util.DistributeDomain for more information
--!								(in UG_BASE/scripts/util/domain_distribution.lua)
--! @param numTargetProcs	(optional, int) The number of target processes to which
--!							the domain shall be distributed. Make sure that the
--!							number of target processes is not higher than the
--!							number of elements in the distributionLevel.
--!							Default is NumProcs()
--!							See util.DistributeDomain for more information
--!							(in UG_BASE/scripts/util/domain_distribution.lua)
--! @param distributionLevel	(optional, int) The level on which the distribution
--!								is performed. Default is the domains top-level
--!								after pre-refinement.
--!								See util.DistributeDomain for more information
--!								(in UG_BASE/scripts/util/domain_distribution.lua)
--! @param wFct 			(optional SmartPtr\<EdgeWeighting\>) Sets the weighting function for the
--!							'metisReweigh' partitioning method.
--! @param noIntegrityCheck	(optional, bool) Disables integrity check if 'true'.
function util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs,
										neededSubsets, distributionMethod,
										verticalInterfaces, numTargetProcs,
										distributionLevel, wFct, noIntegrityCheck)

	-- create Instance of a Domain
	local dom = Domain()
	
	-- load domain
	write("Loading Domain "..gridName.." ... ") 
	LoadDomain(dom, gridName)
	write("done.\n")
	
	if noIntegrityCheck ~= true then
		write("Performing integrity check on domain ... ")
		if CheckForUnconnectedSides(dom:grid()) == true then
			write("WARNING: unconnected sides found (see above).\n")
			local note = "NOTE: You may disable this check by passing 'true' "..
				  		 "to 'noIntegrityCheck' in 'util.CreateAndDistributeDomain'.\n"
			write(note)
			errlog(note)
		end
		write("done.\n")
	end

	-- create Refiner
	ug_assert(numPreRefs <= numRefs, "numPreRefs must be smaller than numRefs. Aborting.");
	
	if numPreRefs > numRefs then
		numPreRefs = numRefs
	end
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner = nil
	if numRefs > 0 then
		refiner = GlobalDomainRefiner(dom)
	end
	
	write("Pre-Refining("..numPreRefs.."): ")
	-- Performing pre-refines
	for i=1,numPreRefs do
		TerminateAbortedRun()
		write(i .. " ")
		refiner:refine()
	end
	write("done.\nDistributing...")
	-- Distribute the domain to all involved processes
	if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		ug_error("Error while Distributing Grid. Aborting.")
	end
	write(" done.\nPost-Refining("..(numRefs-numPreRefs).."): ")
	
	if numRefs > 0 then
		-- Perform post-refine
		for i=numPreRefs+1,numRefs do
			TerminateAbortedRun()
			refiner:refine()
			write(i-numPreRefs .. " ")
		end
	end
	write("done.\n")
	
	-- Now we loop all subsets an search for it in the SubsetHandler of the domain
	if neededSubsets ~= nil then
		if util.CheckSubsets(dom, neededSubsets) == false then 
			ug_error("Something wrong with required subsets. Aborting.");
		end
	end
	
	
	--clean up
	if refiner ~= nil then
		delete(refiner)
	end
	
	-- return the created domain
	return dom
end
