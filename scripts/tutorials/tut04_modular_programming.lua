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

-- Since ug supports a bunch of different algebra modules we will choose one here.
-- This should always be the first thing you do in an ug-script.
-- The cpu-algebra is fine for now.
InitAlgebra(CPUAlgebraChooser())

-- include the basic util-methods.
ug_load_script("../ug_util.lua")

-- include the files in which we defined our helper methods
ug_load_script("tut04_1_domain_util.lua")
