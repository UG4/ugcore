--------------------------------------------------------------------------------
--	tut02_loading_a_domain.lua
--
--	This tutorial will give you a first idea on how to use a lua script to
--	steer ug. We will examine how to use command line arguments to specify
--	a filename and a dimension and on how to load and save a domain.
--
--	Note that we're using methods and classes in this script file which are
--	not part of lua itself, but which have been registered at lua by ug.
--------------------------------------------------------------------------------

-- First we will include a script file which defines some methods often used.
-- You will recognize those methods by a leading util... in the methods name.
ug_load_script("ug_util.lua")

-- Since ug supports a bunch of different algebra modules we will choose one here.
-- The cpu-algebra is fine for now.
InitAlgebra(CPUAlgebraChooser());


-- To keep the script flexible we will now define some variables which have
-- a default value that can be overwritten by command line arguments.

-- Depending on the dimension we will choose our domain object
-- (either 1d, 2d or 3d) and associated discretization objects. Note that
-- we're using some methods defined in "ug_util.lua" here.
dim = GetParamNumber("-dim", 2) -- default dimension is 2.

-- We also need a filename for the grid that shall be loaded.
gridName = GetParam("-grid", "unit_square_quads_8x8.ugx")

-- Since we want to save the domains grid to a file, we also need an output file.
outFileName = GetParam("-o", "domain.ugx")




-- Now its time to create the domain object. We will use an util method here,
-- which automatically creates the right domain for the given dimension.
dom = utilCreateDomain(dim)

-- Now that we have a domain, we can load a grid into it. We check the return
-- value whether the loading was successful or not.
-- Note that we use the method utilLoadDomain instead of LoadDomain. utilLoadDomain
-- has the benefit that grids are automatically searched in the data/grids folder if
-- they were not found at the default locations (execution-path or a path specified
-- in your environments path-variable).
if utilLoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
--	call exit to leave the application right away.
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)



-- Lets save it to outFileName
if SaveDomain(dom, outFileName) == false then
	print("Saving of domain to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)
