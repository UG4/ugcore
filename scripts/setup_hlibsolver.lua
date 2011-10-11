----------------------------------------------------------
--
--   Lua - Script which creates and configures a HLIB solver (linear solver).
--
--   Author: Ingo Heppner (using HLIB solver specific code from 'hlibtest.lua'.
--   Begin: 29082011.
--
--   Provides: Function which set up a HLIB solver object.
--          Input parameters:
--          'linMaxIterations', 'activateDbgWriter', 'verbosity'.
--          Return value:
--          Reference of the fully configured HLIB solver object.
--
--   Requires: ug4 must be configured with cmake parameter '-DHLIBPRO=ON'.
--          See 'ug4/README_hlibpro.txt' for details.
--
--   Usage example:
--	1. In the calling LUA script:
--
--	ug_load_script("setup_hlibsolver.lua")
--	solver = SetupHLIBSolver(lsMaxIter,
--				 numProcs,
--				 activateDbgWriter,
--				 verbosity)
--
--      where 'solver' is the solver object which is called by
--      'ApplyLinearSolver(., ., ., <solver>)', i.e.:
--
--	print("Apply solver.")
--	if ApplyLinearSolver(linOp, u, b, solver) == false then
--		print("Could not apply linear solver.");
--	end
--
--	2. Execution of the calling LUA script:
--[[
# Local:
########
setenv UGARGS "-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 0 -lsType hlib"
openmpirun -np 4 ugshell $UGARGS -numRefs 5

# cekon (if configured with HLIB, see requirements above):
########
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 0 -lsType hlib"

salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 5
]]

----------------------------------------------------------
-- function 'SetupHLIBSolver()':
----------------------------------------------------------
function SetupHLIBSolver(linMaxIterations,
			 activateDbgWriter,
			 dbgWriter,
			 verbosity)

	print("    'setup_hlibsolver.lua': Setting up HLIB solver...")

	-- Stuff for the HLIB solver:
        local hlibSolver

        local hlibConvCheck
	local linAbsLimit  = 1e-7
	local linReduction = 1e-16

	--------------------------------------------------------------------------------
	-- Checking for HLIB solver related parameters (begin)
	--------------------------------------------------------------------------------
	if util.HasParamOption("-geom") == true then
		clustering = "geom"
	else
		clustering = "alg"
	end

	print(" HLIB parameters chosen:")
	print("    clustering =        " .. clustering)

	----------------------------------------------------------
	-- create and configure HLIB Solver
	----------------------------------------------------------
	hlibSolver = HLIBSolver()
	hlibSolver:set_hlib_accuracy_H(1.e-4)  -- default: 1.e-4
	hlibSolver:set_hlib_accuracy_LU(1.e-4) -- default: 1.e-4

	-- define construction of cluster tree
	--   first  arg: "clustering type" \in [algebraic | geometric (not yet implemented)]; algebraic is default
	--   second arg: "clustering mode" \in [nested dissection | empty/everything else]; nested dissection is default 
	hlibSolver:set_clustering_method("algebraic", "nested dissection")

	hlibSolver:set_hlib_verbosity(4) -- '>= 2': create HLIB related postscripts; '>=3' also create plots of matrix entries
	hlibSolver:set_ps_basename("hlib")

	-- define convergence criteria for the HLIB solver
	hlibConvCheck = StandardConvergenceCheck()

	print("    'setup_hlibsolver.lua': linMaxIterations = " .. linMaxIterations .. ", linAbsLimit = " .. linAbsLimit .. ", linReduction = " .. linReduction)

	hlibConvCheck:set_maximum_steps(linMaxIterations)
	hlibConvCheck:set_minimum_defect(linAbsLimit)
	hlibConvCheck:set_reduction(linReduction)
	--hlibConvCheck:set_verbose_level(false)

	hlibSolver:set_convergence_check(hlibConvCheck)

	if activateDbgWriter >= 1 then
		print( "activating debug writer for hlibSolver")
		hlibSolver:set_debug(dbgWriter)
	end

	----------------------------------------------------------
	print("    'setup_hlibsolver.lua': returning HLIB solver 'hlibSolver', ready for application!")
	return hlibSolver
end
