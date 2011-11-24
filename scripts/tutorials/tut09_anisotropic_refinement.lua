--------------------------------------------------------------------------------
--	tut09_anisotropic_refinement.lua
--
--	The purpose of this tutorial is to demonstrate how one can use adaptive
--	refinement to refine geometries in a more sophisticated way.
--	The tutorial is very similar to tut05_global_refinement.lua
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("../ug_util.lua")

-- Get the command line parameters
dim = util.GetParamNumber("-dim", 2)

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

if dim == 2 then
	gridName = util.GetParam("-grid", "anisotropic_rect.ugx")
elseif dim == 3 then
	gridName = util.GetParam("-grid", "anisotropic_hexa.ugx")
else
	print("Only dim == 2 and dim == 3 are supported in the moment")
	exit()
end

outFileNamePrefix = util.GetParam("-o", "distributed_domain_")

-- We will save the created hierarchy to this file (with appended process id)
outHierarchyFilePrefix = util.GetParam("-oh", "hierarchy_on_proc_")

-- We additionally use parameters which allow to specify the number of
-- pre- and total-refinement steps (wrt domain distribution).
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
numTotalRefs = util.GetParamNumber("-numTotalRefs", 3)

-- Calculate the number of post-refs and make sure that the result makes sense.
numPostRefs = numTotalRefs - numPreRefs
if numPostRefs < 0 then
	print("WARNING:\tnumPreRefs exceeds the number of total refinements.")
	print("\t\t\tNo refinement will be preformed after distribution.")
	numPostRefs = 0
end

--	The edgeRatio determines which faces are considered to be anisotropic.
--	It should always be in the interval [0, 1].
edgeRatio = util.GetParamNumber("-edgeRatio", 0.5)
if edgeRatio < 0 or edgeRatio > 1 then
	print("edgeRatio should be in the interval [0, 1]")
	exit()
end

-- Create the domain and load a grid
dom = Domain()

if LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)


-- We will create a refiner now. Since we want to use adaptive refinement
-- we'll use a hanging node refiner
refiner = HangingNodeDomainRefiner(dom)

-- perform pre-refinement
for i = 1, numPreRefs do
--	we'll mark anisotropic elements for anisotropic refinement here.
--	Note that many different MarkForRefinement_... algorthms exist.
--	e.g. MarkForRefinement_VerticesInSphere or MarkForRefinement_VerticesInCube.
	MarkForRefinement_AnisotropicElements(dom, refiner, edgeRatio)
	refiner:refine()
end


-- Distribute the refined domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end


-- perform post-refinement
for i = 1, numPostRefs do
	MarkForRefinement_AnisotropicElements(dom, refiner, edgeRatio)
	refiner:refine()
	PrintGridElementNumbers(dom:grid())
end


-- Lets save the domain on each process
outFileName = outFileNamePrefix .. GetProcessRank() .. ".ugx"
if SaveDomain(dom, outFileName) == false then
	print("Saving of domain to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Everything seems to went fine.
print("Saved domain to " .. outFileName)


-- Now lets save the hierarchy on each process
-- The SaveGridHierarchy routine directly works on the domains grid.
-- We can access the grid of a domain through its grid() member method.
--
-- SaveGridHierarchy outputs a grid, where each level is assigned to a subset.
-- Original subsets are not contained in that file.
outFileName = outHierarchyFilePrefix .. GetProcessRank() .. ".ugx"
if SaveGridHierarchy(dom:grid(), outFileName) == false then
	print("Saving of grid-hierarch to " .. outFileName .. " failed. Aborting.")
	exit()
end

-- Again everything went fine.
print("Saved hierarchy to " .. outFileName)

-- we're done.
print("")
print("done")

