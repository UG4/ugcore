-- creates a ugx file
-- author: Martin Scherer

ug_load_script("ug_util.lua")

-- default parameter a, w, h are unit tkd 
a = util.GetParamNumber("-a", 10)
w = util.GetParamNumber("-w", 30)
h = util.GetParamNumber("-h", math.sqrt(6)*a)
d_lipid = util.GetParamNumber("-dl", 3)

if not (w > 2*a) then
	print("Geometric constraint w > 2a not met. Check parameters!")
	exit()
end

-- parameter for stacking the tkd's
rows = util.GetParamNumber("-rows", 5)
cols = util.GetParamNumber("-cols", 4)
layers = util.GetParamNumber("-layers", 3)

-- file name to write grid to
filename = util.GetParam("-filename", "/tmp/tkd/tkdDomain.ugx")

-- subset information
sh_name_corneocyte = util.GetParam("-sh_name_corneo", "corneocytes")
sh_name_lipid = util.GetParam("-sh_name_lipid", "corneocytes")

-- init grid and its associated subset handler
go = GridObject()
grid = go:grid()
sh = go:subset_handler()

-- create the grid
CreateTKDDomain(grid, sh, a, w, h, d_lipid, rows, cols, layers)

-- optionally override subset information
green = MakeVec(0, 1, 0, 0)
blue = MakeVec(0, 0, 1, 0)
SetTKDSubsetInfos(sh, sh_name_corneocyte, sh_name_lipid, blue, green)

if SaveGrid(grid, sh, filename) then
	print("TKD Domain successfully written to " .. filename)
else
	print("Problem writing to given file: " .. filename)
end