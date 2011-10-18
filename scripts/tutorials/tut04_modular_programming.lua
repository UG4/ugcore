--------------------------------------------------------------------------------
--	tut04_modular_programming.lua
--
--	This tutorial will not introduce new numerical routines. Instead it will
--	show you how you can organize your code to make it easy for others to
--	reuse it. Note that similar but slightly extended methods as the ones we
--	create in this tutorial can also be found in scripts/prefabs and scripts/util,
--	where they are located in the "discs" and "util" namespaces.
--
--	Note that this tutorial spans over multiple files:
--	* tut04_modular_programming.lua
--	* tut04_1_domain_util.lua
--	* tut04_2_disc_laplace.lua
--
--	Note that the code runs in a parallel evironment too. However the exact
--	LU decomposition can only be used for serial environments. One would
--	thus have to replace the solver before the code can be executed in parallel.
--	See the next tutorials for an overview of alternate solvers.
--------------------------------------------------------------------------------

-- include the basic util-methods.
ug_load_script("../ug_util.lua")

-- include the files in which we defined our helper methods
ug_load_script("tut04_1_domain_util.lua")
ug_load_script("tut04_2_disc_laplace.lua")


-- Get the command line parameters
dim = util.GetParamNumber("-dim", 2)

-- Since ug supports a bunch of different dimensions and algebra modules 
-- we will choose a combination here. This should always be the first thing 
-- you do in an ug-script. The cpu-algebra is fine for now.
InitUG(dim, AlgebraType("CPU", 1))

gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")
outFileNamePrefix = util.GetParam("-o", "distributed_domain_")


-- Create the domain through the CreateAndDistributeDomain method,
-- which we created in tut04_1_domain_util.lua
dom = CreateAndDistributeDomain(gridName, dim, outFileNamePrefix) 

-- Create surface functions (vectors) for Au=b (where A=linOp) and initialize them
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)

-- Using the AssembleLaplace method which we created in
-- tut04_2_disc_laplace.lua, we now create the linear operator.
-- Note that we use the default callbacks defined in tut04_2_disc_laplace.lua
-- by passing nil. One could instead simply not specify the arguments
-- (this defaults to nil).
-- You could of course create your callbacks analogous to tutorial 3 and
-- pass their names (as strings) to AssembleLaplace.
-- Make sure that you use the ones with the right dimension!
linOp, approxSpace = AssembleLaplace(dom, "Inner", "Boundary", b, nil, nil, nil, nil)

u:set(0)

-- We again choose the LU solver and make sure that it is only used in a
-- serial environment
if GetNumProcesses() > 1 then
	print("Can't apply LU decomposition in parallel environment. Aborting.")
	exit()
end

-- Create the solver and initialize it with the linear operator
solver = LU()
solver:init(linOp)

-- The solver is ready. We can now apply it.
solver:apply(u, b)

-- Write the result to a file.
WriteGridFunctionToVTK(u, "Solution")

-- We're done.
print("")
print("done.")
