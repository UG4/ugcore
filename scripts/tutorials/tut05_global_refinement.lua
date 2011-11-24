--------------------------------------------------------------------------------
--	tut05_global_refinement.lua
--
--	In this tutorial we would like to show you how a domain can be refined
--	globally. This script again works in a parallel environment.
--	We will examine how the underlying grid can be refined prior to distribution
--	and after distribution. Note that ug4 also supports adaptive refinement.
--	This however is subject to a later tutorial.
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("../ug_util.lua")

-- Get the command line parameters
dim = util.GetParamNumber("-dim", 2)

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")
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



-- Create the domain and load a grid
dom = Domain()

if LoadDomain(dom, gridName) == false then
	print("Loading of domain " .. gridName .. " failed. Aborting.")
	exit() 
end

-- Now that we're here, the domain was successfully loaded
print("Loaded domain from " .. gridName)


-- We will create a refiner now. This is a tool associated with a domain.
-- UG defines factory methods for refiners, which automatically choose
-- the right refiner for the given context, i.e. different refiners are
-- created depending on whether we are in a parallel or in a serial environment.
-- Note that another factory method is HangingNodeDomainRefiner, which is
-- subject to a later tutorial.
refiner = GlobalDomainRefiner(dom)

-- perform pre-refinement
for i = 1, numPreRefs do
	refiner:refine()
end


-- Distribute the refined domain to all involved processes
if DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end


-- perform post-refinement
for i = 1, numPostRefs do
	refiner:refine()
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
