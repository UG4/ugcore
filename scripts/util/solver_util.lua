--[[!
\file solver_util.lua
\defgroup scripts_util_solver Solver Utility
\ingroup scripts_util
\brief Utility functions to easily create solver toolchains
\note Use these functions to reduce copy-and-paste code in your Lua scripts.
\{

<hr>
\b GetPreconditioner

Creates a preconditioner based on an abbreviated name and additional, optional
parameters for that preconditioner.

\em Parameters
- \c name -
  Required. Abbreviation of the preconditioner. Possible options:
  - \c jac - Jacobi
  - \c gs - Gauß-Seidel
  - \c sgs - symmetric Gauß-Seidel
  - \c bgs - backward Gauß-Seidel
  - \c ilu - ILU
  - \c ilut - ILUT
  - \c gmg - Geometric Multi-Grid
- \c jac_damp -
  damping factor for Jacobi solver (default: 0.6)
- \c gmg_approxSpace -
  approximation space for GMG solver (required)
- \c gmg_disc -
  discretisation for GMG solver
- \c gmg_base -
  base solver for GMG solver (default: "lu"). 
- \c gmg_baseLevel -
  base level for GMG solver (default: 0)
- \c gmg_parallelBase -
  whether the base solver of GMG should run in parallel (default: \c false)
- \c gmg_smoother -
  name of the smoother for GMG (default: "ilu")
- \c gmg_cycleType -
  cycle type of the GMG (default: 1)
- \c gmg_numPreSmooth -
  number of pre smoothing steps of GMG (default: 2)
- \c gmg_numPostSmooth -
  number of post smoothing steps of GMG (default: 2)

\note The values for the smoother and base solver can be given as well as tables
  with the same options as for this function (util.GetPreconditioner is getting
  called recursively) or as previously created (ug4)-objects.

\em Examples
\code{.lua}
// standard ILU preconditioner
iluPrecond = util.GetPreconditioner( "ilu" ) // or util.GetPreconditioner( { name="ilu" } )

// Jacobi with damping factor of 0.75
jacPrecond = util.GetPreconditioner( {
  name = "jac",
  jac_damp = 0.75
} )

// GMG with LU base solver (default), ILU as smoother (default), V-cycle (default) 
// and 1 soomthing step each
gmgPrecond = util.GetPreconditioner( {
  name = "gmg",
  gmg_approxSpace = approxSpace, // obviously this has to be created beforehand
  gmg_numPreSmooth = 1,
  gmg_numPostSmooth = 1
} )

// GMG with LU as base solver (default), Jacobi (damp=0.8) as smoother
gmgPrecond = util.GetPreconditioner( {
  name = "gmg",
  gmg_approxSpace = approxSpace,
  gmg_smoother = {
    name = "jac",
    jac_damp = 0.8 }
} )
\endcode

<hr>
\b GetConvCheck

Prepares the convergence check

\em Parameters

- \c default -
  \c true for using defaults for all parameters (\c false will have no effect)
- \c maxSteps -
  maximum number of steps (default: 50)
- \c reduction -
  desired reduction (default: 1e-10)
- \c minDef -
  desired minimal defect (default: 1e-7)

<hr>
\b GetSolver

Creates a solver based on a name with given preconditioner name and convergence check

\em Parameters
- \c name -
  Required. Abbreviated name of the solver. Possible options:
  - \c linear - linear solver
  - \c bicgstab - BiCGStab
  - \c lu - LU
- \c precond -
  preconditioner, which can be defined by parameters for util.GetPreconditioner
  or as a pre-created (ug)-object.
  Default: ILU
- \c convCheck -
  convergence check, which can be defined by parameters for util.GetConvCheck
  or as a pre-created (ug)-object
  Default: default ConvCheck (see util.GetConvCheck)

\em Examples
\code{.lua}
// LU solver
luSolver = util.GetSolver( "lu" ) // or util.GetSolver( { name = "lu" } )

// Linear MG Solver with GMG preconditioner (with Jacobi and ILU)
mgSolver = util.GetSolver( {
  name = "linear",
  precond = {
    name = "gmg",
    gmg_approxSpace = approxSpace,
    gmg_smoother = {
      name = "jac",
      jac_damp = 0.8 } },
  convCheck = {
    maxSteps = 100,
    minDef = 1e-9,
    reduction = 1e-12 }
} )
\endcode
]]--

util = util or {}

--[[
Creates a preconditioner based on an abbreviated name
]]--
util.GetPreconditioner = util.CreateFancyFunction(
	{
		{"name"},
		{"jac_damp", 0.6},
		{"gmg_approxSpace"},
		{"gmg_disc"},
		{"gmg_base", "lu"},
		{"gmg_baseLevel", 0},
		{"gmg_parallelBase", false},
		{"gmg_smoother", "ilu"},
		{"gmg_cycleType", 1},
		{"gmg_numPreSmooth", 2},
		{"gmg_numPostSmooth", 2},
		{"gmg_restriction_post_process"},
		{"gmg_prolongation_post_process"}
	},
	function( name, jac_damp, gmg_approxSpace, gmg_disc, gmg_base, gmg_baseLevel, gmg_parallelBase, gmg_smoother, gmg_cycleType, gmg_numPreSmooth, gmg_numPostSmooth )
		if not name then
			print( "Specifying the name of the preconditioner is mandatory. Aborting." )
			exit()
		end
		
-- 		print( "DBG >> Calling util.GetPreconditioner with name="..name )
		
		if name == "jac" then
			jac = Jacobi()
			jac:set_damp( jac_damp )
			return jac
		elseif name == "gs" then
			return GaussSeidel()
		elseif name == "sgs" then
			return SymmetricGaussSeidel()
		elseif name == "bgs" then
			return BackwardGaussSeidel()
		elseif name == "ilu" then
			ilu = ILU()
			return ilu
		elseif name == "gmg" then
			if not gmg_approxSpace then
				print( "Approximation space must be defined. Aborting." )
				exit()
			end
			if not gmg_disc then
				print( "WARNING: Not specifying the discretization may cause errors. Use 'gmg_disc'." )
			end
			
			-- get base solver
			if type(gmg_base) == "table" then
-- 				print( "DBG >>   GMG Base Solver from table:" )
-- 				util.PrintTable(gmg_base)
				baseSolver = util.GetSolver( gmg_base )
			elseif type(gmg_base) == "string" then
				baseSolver = util.GetSolver( { name = gmg_base } )
			else
				baseSolver = gmg_base
			end
			
			-- get Smoother
			if type(gmg_smoother) == "table" then
-- 				print( "DBG >>   GMG Smoother from table:" )
-- 				util.PrintTable(gmg_smoother)
				smoother = util.GetPreconditioner( gmg_smoother )
			elseif type(gmg_smoother) == "string" then
				smoother = util.GetPreconditioner( { name = gmg_smoother } )
			else
				smoother = gmg_smoother
			end
			
			-- Geometric Multi-Grid
			gmg = GeometricMultiGrid( gmg_approxSpace )
			if gmg_disc then gmg:set_discretization( gmg_disc ) end
			gmg:set_base_level( gmg_baseLevel )
			gmg:set_parallel_base_solver( gmg_parallelBase )
			gmg:set_base_solver( baseSolver )
			gmg:set_smoother( smoother )
			gmg:set_cycle_type( gmg_cycleType )
			gmg:set_num_presmooth( gmg_numPreSmooth )
			gmg:set_num_postsmooth( gmg_numPostSmooth )
			if gmg_restriction_post_process then
				gmg:add_restriction_post_process( gmg_restriction_post_process )
			end
			if gmg_prolongation_post_process then
				gmg:add_prolongation_post_process( gmg_prolongation_post_process )
			end
			return gmg
			
		elseif name == "ilut" then
			return ILUT()
		end
	end
)

--[[!
Prepares the convergence check
]]--
util.GetConvCheck = util.CreateFancyFunction(
	{
		{"default"},
		{"maxSteps", 50},
		{"reduction", 1e-10},
		{"minDef", 1e-7}
	},
	function( default, maxSteps, reduction, minDef )
-- 		print( "DBG >> Calling util.GetConvCheck with maxSteps="..maxSteps..", reduction="..reduction..", minDef="..minDef )
		
		convCheck = ConvCheck()
		convCheck:set_maximum_steps( maxSteps )
		convCheck:set_reduction( reduction )
		convCheck:set_minimum_defect( minDef )
		return convCheck
	end
)

--[[
Creates a solver based on a name with given preconditioner and convergence check
]]--
util.GetSolver = util.CreateFancyFunction(
	{
		{"name"},
		{"precond", "ilu" },
		{"convCheck", "default" }
	},
	function( name, precond, convCheck )
		if not name then
			print( "Specifying the name of the solver is mandatory. Aborting." )
			exit()
		end
		
-- 		print( "DBG >> Calling util.GetSolver with name="..name )
		
		-- get preconditioner
		if type(precond) == "table" then
-- 			print( "DBG >>   Preconditioner from table:" )
-- 			util.PrintTable(precond)
			precondObj = util.GetPreconditioner( precond )
		elseif type(precond) == "string" then
			precondObj = util.GetPreconditioner( { name = precond } )
		else
			precondObj = precond
		end
		
		-- get convergence check
		if convCheck == "default" then
			convCheckObj = util.GetConvCheck( { default = true } )
		elseif type(convCheck) == "table" then
-- 			print( "DBG >>   ConvergenceCheck from table:" )
-- 			util.PrintTable(convCheck)
			convCheckObj = util.GetConvCheck( convCheck )
		else
			convCheckObj = convCheck
		end
		
		if name == "linear" then
-- 			print( "DBG >>   Creating linear solver" )
			linSolver = LinearSolver()
			linSolver:set_preconditioner( precondObj )
			linSolver:set_convergence_check( convCheckObj )
			return linSolver

		elseif name == "bicgstab" then
-- 			print( "DBG >>   Creating BiCGStab solver" )
			bicgstab = BiCGStab()
			bicgstab:set_preconditioner( precondObj )
			bicgstab:set_convergence_check( convCheckObj )
			return bicgstab

		elseif name == "lu" then
-- 			print( "DBG >>   Creating LU solver" )
			return LU()
		end
	end
)

-- end group scripts_util_solver
--[[! \} ]]--
