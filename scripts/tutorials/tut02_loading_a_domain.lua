--------------------------------------------------------------------------------
--	tut02_loading_a_domain.lua
--
--	This tutorial will give you a first idea on how to use a lua script to
--	steer ug. We will examine how to use command line arguments to specify
--	a filename and a dimension and on how to load and save a domain.
--	Furthermore we will distribute the domain to all involved processes.
--
--	Note that we're using methods and classes in this script file which are
--	not part of lua itself, but which have been registered at lua by ug.
--------------------------------------------------------------------------------

-- We will include a script file which defines some methods often used.
-- Loaded methods are all found in the library util.
ug_load_script("../ug_util.lua")

-- To keep the script flexible we will now define some variables which have
-- a default value that can be overwritten by command line arguments.

-- Depending on the dimension we will choose our domain object
-- (either 1d, 2d or 3d) and associated discretization objects. Note that
-- we're using some methods defined in "ug_util.lua" here.
dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

-- We also need a filename for the grid that shall be loaded.
gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")





-- Now its time to create the domain object. We will use an util method here,
-- which automatically creates the right domain for the given dimension.
dom = Domain()

-- Now that we have a domain, we can load a grid into it. We check the return
-- value whether the loading was successful or not.
-- Note that we use the method LoadDomain instead of LoadDomain. utilLoadDomain
-- has the benefit that grids are automatically searched in the data/grids folder if
-- they were not found at the default locations (execution-path or a path specified
-- in your environments path-variable).
if LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
--	call exit to leave the application right away.
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)


-- Distribute the domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end


-- Lets save the domain on each process
outFileName = "distributedDomainOnProc" .. GetProcessRank() .. ".ugx"
if SaveDomain(dom, outFileName) == false then
	print("Saving of domain to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)
