-- created by Sebastian Reiter
-- s.b.reiter@googlemail.com

gridInName = "unit_square_quads_16x16.obj"

go = GridObject()
if LoadGridObject(go, gridInName) == false then
	printf("grid " .. gridInName .. " not found.")
end

hr = HangingNodeRefiner()
hr:assign_grid(go:get_grid())
CreateFractal(go:get_grid(), hr, 0.6, 4)
SaveGridObject(go, "tmp.ugx")
