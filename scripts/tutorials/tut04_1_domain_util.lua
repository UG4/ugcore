--------------------------------------------------------------------------------
--	tut04_1_domain_util.lua
--
--	This file is part of tut04_modular_programming.lua.
--
--	This file contains the method CreateAndDistributeDomain.
--	The method which will create a domain, load the associated
--	geometry from a file and distriutes it onto all active processes.
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("ug_util.lua")


--------------------------------------------------------------------------------
--	The method loads creates a new domain with the given dimension (dim),
--	loads the grid specidied in gridName, distributes it to all active
--	processes and finally, if savePrefix ~= nil, it will save the received
--	grid on each process to a file with name:
--	savePrefix .. "_" .. GetProcessRank() .. ".ugx"
--	If everything went right, the method returns the created domain object,
--	if an error occured, nil is returned.
--	
--	Params: string gridName, int dim, string savePrefix
--	Returns: Domain
--
function CreateAndDistributeDomain(gridName, dim, savePrefix)
--	Create the domain object
	local dom = Domain()
	
--	Load the domain from file (Note that this is always performed on
--	process 0)
	if LoadDomain(dom, gridName) == false then
		print("Loading of domain " .. gridName .. " failed.")
		return nil
	end
		
--	Distribute the domain to all involved processes
	if DistributeDomain(dom) == false then
		print("Error while Distributing Domain.")
		return nil
	end
	
--	If a savePrefix was specified, we'll save the domain on each process
	if savePrefix ~= nil then
		local outFileName = savePrefix .. "_" .. GetProcessRank() .. ".ugx"
		if SaveDomain(dom, outFileName) == false then
			print("Saving of domain to " .. outFileName .. " failed.")
			return nil
		end
	end
	
--	we're done. Return the domain.
	return dom
end
