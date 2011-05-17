-- created by Sebastian Reiter
-- s.b.reiter@googlemail.com

ug_load_script("ug_util.lua")

gridInName = util.GetGridPath() .. "unit_square/unit_square_quads_8x8.ugx"

go = GridObject()

if LoadGridObject(go, gridInName) == false then
	print("grid " .. gridInName .. " not found.")
end

hr = HangingNodeRefiner_Grid()
hr:assign_grid(go:get_grid())
CreateFractal(go:get_grid(), hr, 0.6, 4)
SaveGridObject(go, "tmp.ugx")
