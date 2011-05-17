----------------------------------------------------------
--
--   Interface Test - Testing for parallel conversions using interfaces
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

-- make sure that ug_util is in the right path.
-- currently only the path in which you start your application is valid.
ug_load_script("ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraSelector());

-- constants
dim = 3

if dim == 2 then
	gridName = "unit_square/unit_square_tri_2x2.ugx"
end
if dim == 3 then
	gridName = "unit_square/unit_cube_hex.ugx"
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

numPreRefs = 0
numRefs = 1

-- create Instance of a Domain
dom = util.CreateDomain(dim)
sh = dom:get_subset_handler()

-- load domain
if util.LoadDomain(dom, gridName) == false then
   print("Loading Domain failed.")
   exit()
end

-- create Refiner
print("Create Refiner")
if numPreRefs > numRefs then
	print("numPreRefs must be smaller/equal than numRefs");
	exit();
end

refiner = GlobalMultiGridRefiner()
refiner:assign_grid(dom:get_grid())
for i=1,numPreRefs do
	refiner:refine()
end

if util.DistributeDomain(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

for i=numPreRefs+1,numRefs do
	util.GlobalRefineParallelDomain(dom)
end

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = util.CreateApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init()

-- get grid function
u = approxSpace:create_surface_function()

--------------------------
-- some parallel testing
--------------------------

-- set to consistent and u == 1.0
u:set(1.0)
WriteGridFunctionToVTK(u, "Sol_Constistent_1");
SaveVectorForConnectionViewer(u, "Sol_Constistent_1.mat")

-- change to unique and write (should be zero in slaves now)
u:change_storage_type_by_string("unique");
WriteGridFunctionToVTK(u, "Sol_Unique_1");
SaveVectorForConnectionViewer(u, "Sol_Unique_1.mat")

-- reset vector to consistent and u == 2.0
u:set(2.0)
WriteGridFunctionToVTK(u, "Sol_Constistent_2");
SaveVectorForConnectionViewer(u, "Sol_Constistent_2.mat")

-- manually set type to additive __without__ changing values in vector
-- this is a somehow forbidden action for real simulations but good for testing
u:set_storage_type_by_string("additive");

-- change type to unique (should be 2 * (num_slave_copies+1) now in interface nodes)
u:change_storage_type_by_string("unique");
WriteGridFunctionToVTK(u, "Sol_Unique_2");
SaveVectorForConnectionViewer(u, "Sol_Unique_2.mat")
