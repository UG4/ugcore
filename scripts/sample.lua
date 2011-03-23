-- print some info about the executed script
print("----------------")
print("- sample.lua - created by Sebastian Reiter.")
print("- execute this script like that:")
print("- ugshell -ex sample.lua gridFileIn gridFileOut")
print("----------------")
print("")

-- print ugargv
for i = 1, ugargc do
	print("ugargv[" .. i .. "]: " .. ugargv[i])
end
print("")

-- use the default Algebra
InitAlgebra(CPUAlgebraSelector());

-- create a new grid and a new subset handler
g = Grid.new()
sh = SubsetHandler.new()

-- assign the grid to the subset handler
-- note that member-methods are called using :
sh:set_grid(g)

-- load a grid from file
-- we'll use the first parameter as filename
loadSucceeded = 0
if ugargv[1] ~= nil then
	print("loading grid from " .. ugargv[1])
	loadSucceeded = ug.load_grid(g, sh, ugargv[1])
	if loadSucceeded == 0 then
		print("load failed.")
	else
		print("load succeeded.")
	end
end

-- save the grid to another file
-- we'll use the second parameter as filename
if ugargv[2] ~= nil and loadSucceeded == 1 then
	print("saving grid to " .. ugargv[2])
	if ug.save_grid(g, sh, ugargv[2]) == 0 then
		print("save failed.")
	else
		print("save succeeded.")
	end
end

-- since we don't need the grid and the subset-handler any more,
-- we'll delete them.
ug.delete(sh)
ug.delete(g)
