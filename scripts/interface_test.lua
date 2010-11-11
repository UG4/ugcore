----------------------------------------------------------
--
--   Interface Test - Testing for parallel conversions using interfaces
--
--   Author: Andreas Vogel
--
----------------------------------------------------------

-- make sure that ug_util is in the right path.
-- currently only the path in which you start your application is valid.
dofile("../scripts/ug_util.lua")

-- choose algebra
InitAlgebra(CPUAlgebraChooser());

-- constants
dim = 3

if dim == 2 then
	gridName = "unit_square_tri.ugx"
end
if dim == 3 then
	gridName = "unit_cube_hex.ugx"
	--gridName = "unit_cube_tets_regular.ugx"
end

numPreRefs = 0
numRefs = 1

-- create Instance of a Domain
print("Create Domain.")
dom = utilCreateDomain(dim)
sh = dom:get_subset_handler()

-- load domain
print("Load Domain from File.")
if utilLoadDomain(dom, gridName) == false then
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

if utilDistributeDomain(dom) == false then
	print("Error while Distributing Grid.")
	exit()
end

print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	utilGlobalRefineParallelDomain(dom)
end

-- create function pattern
print("Create Function Pattern")
pattern = P1ConformFunctionPattern()
pattern:set_subset_handler(sh)
AddP1Function(pattern, "c", dim)
pattern:lock()

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = utilCreateApproximationSpace(dom, pattern)

-- get grid function
u = approxSpace:create_surface_function("u", true)
b = approxSpace:create_surface_function("b", true)

-- some parallel testing
u:set(1.0)
WriteGridFunctionToVTK(u, "Sol_Constistent_1");
SaveVectorForConnectionViewer(u, "Sol_Constistent_1.mat")

u:change_storage_type_by_string("unique");
WriteGridFunctionToVTK(u, "Sol_Unique_1");
SaveVectorForConnectionViewer(u, "Sol_Unique_1.mat")

u:set(2.0)
WriteGridFunctionToVTK(u, "Sol_Constistent_2");
SaveVectorForConnectionViewer(u, "Sol_Constistent_2.mat")

u:set_storage_type_by_string("additive");
u:change_storage_type_by_string("unique");

WriteGridFunctionToVTK(u, "Sol_Unique_2");
SaveVectorForConnectionViewer(u, "Sol_Unique_2.mat")
