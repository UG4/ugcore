-- Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
-- Authors: Sebastian Reiter, Andreas Vogel
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.



--[[!
\defgroup solver_util Solver Util
\ingroup scripts_util
\brief Table based solver creation

Given a (possibly nested) solver description, e.g.,
\code
	solver = {
		type = "bicgstab",
		precond = "ilu",
		convCheck = "standard"
	}
\endcode
a call e.g. to util.solver.CreateSolver(solver)
creates and returns the requested solver.

All components may either be specified through a string-id or through a
descriptor. You only have to specify those components which deviate from
the defaults (solver.util.defaults), the rest will be ammended.

To use a geometric multigrid preconditioner, you either have
to specify an approxSpace entry in its descriptor or set
'util.solver.defaults.approxSpace' before calling any 'util.solver.Create...'
method.


<br>
<h2>Examples</h2>
Several examples follow:

<h3>Example: Creating an exact solver</h3>
\code
	solver = util.solver.CreateSolver("lu")
\endcode


<br>
<h3>Example: Creating a simple iterative jacobi solver</h3>
\code
	solverDesc = {
		type = "linear",
		precond = "jac"
	}

	solver = util.solver.CreateSolver(solverDesc)
\endcode


<br>
<h3>Example: Creating a simple iterative jacobi solver with custom parameters</h3>
\code
	solverDesc = {
		type = "linear",
		precond = {
			type = "jac",
			damping = 0.5
		}
	}

	solver = util.solver.CreateSolver(solverDesc)
\endcode


<br>
<h3>Example: Creation of a bicgstab solver with geometric multigrid preconditioner</h3>
here someApproxSpace points to a previously created approximation space
\code
	solverDesc = 
	{
		type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
		precond = 
		{	
			type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
			smoother 	= "gs",		-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
			cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
			preSmooth	= 3,		-- number presmoothing steps
			postSmooth 	= 3,		-- number postsmoothing steps
			rap			= false,	-- comutes RAP-product instead of assembling if true 
			baseLevel	= 0,		-- gmg - baselevel
			baseSolver	= {			-- better options are most likely "lu" or "superlu"
				type	  = "bicgstab",
				precond	  = "gs",
				convCheck = {
					type		= "standard",
					iterations	= 1000,
					absolute	= 1e-12,
					reduction	= 1e-10,	-- higher values may suffice and may be much
											-- more efficent (e.g. 1e-2, 1e-4, ...).
					verbose		= false
				}
			},
			approxSpace	= someApproxSpace,
		},

		convCheck = "standard"
	}
	solver = util.solver.CreateSolver(solverDesc)
\endcode


<br>
<h3>Example: Providing custom solver and preconditioner objects (advanced)</h3>
You may also provide instances of preconditioners instead of descriptor tables:
\code
	solverDesc =
	{
		type = "linear",
		precond = ILU()		-- for each component one can optionally supply a custom instance
	}
	solver = util.solver.CreateSolver(solverDesc)
\endcode


<br>
<h3>Example: Accessing the created solvers and preconditioners (advanced)</h3>
For conveniance, all the methods will add an 'instance' entry to their respective
descriptor tables, assigning the created solver component.
If you thus setup a solver, e.g.
\code
	solver = {
		type = "bicgstab",
		precond = {
			type = "ilu"
		}
	}
\endcode
you may later access 'solver.instance' and 'solver.precond.instance'.
Note that this is only possible if you specified a component through a descriptor.
If you use the short form, e.g.
\code
	solver = {
		type = "bicgstab",
		precond = "ilu"
	}
\endcode
you will still be able to access 'solver.instance', but you can't access
'solver.precond.instance'.


<br>
<h2>Non-Linear Solvers</h2>
The following listing gives an overview over available non-linear solvers and their default
parameters. Pass the descriptor of a non-linear solver to
\code
	util.solver.CreateSolver(...)
\endcode
to obtain the described solved.

\b linSolver can be any linear solver listed in the <b>Linear Solvers</b> section.

\b convCheck can be any convergence check listed in the <b>Convergence Checks</b> section.

\b lineSearch can be any line search method listed in the <b>Line Search</b> section.

Currently only the Newton method is available as non-linear solver.

<h3>Newton Method</h3>
\code
{
	type		= "newton",
	convCheck	= "standard",
	linSolver	= "bicgstab",
	lineSearch	= nil
}
\endcode


<br>
<h2>Linear Solvers</h2>
The following listing gives an overview over available linear-solvers and their default
parameters. Pass the descriptor of a linear solver to
\code
	util.solver.CreateSolver(...)
\endcode
to obtain the described solved.

\b precond can be any preconditioner listed in the \b Preconditioners section.

\b convCheck can be any convergence check listed in the <b>Convergence Checks</b> section.

<h3>Linear Solver</h3>
\code
{
	type		= "linear",
	precond		= "ilu",
	convCheck	= "standard"
}
\endcode

<h3>Conjugate Gradients</h3>
\code
{
	type		= "cg",
	precond		= "ilu",
	convCheck	= "standard"
}
\endcode

<h3>BiCGStab</h3>
\code
{
	type		= "bicgstab",
	precond		= "ilu",
	convCheck	= "standard"
}
\endcode


<br>
<h2>Preconditioners</h2>
The following listing gives an overview over available preconditioners and their default
parameters.

<h3>Geometric Multigrid</h3>
\code
{
	type 							= "gmg",
	adaptive 						= false,
	approxSpace 					= nil,
	baseLevel 						= 0,
	baseSolver 						= "lu",		-- any solver listed in the 'Solvers' section
	cycle 							= "V",
	discretization 					= nil,		-- only necessary if the underlying matrix is not of type AssembledLinearOperator
	gatheredBaseSolverIfAmbiguous	= false,
	preSmooth 						= 3,
	postSmooth 						= 3,
	rap 							= false,
	rim 							= false,
	emulateFullRefined 				= false,
	smoother 						= "gs",		-- any preconditioner listed in the 'Preconditioners' section
	transfer 						= "std",	-- any transfer listed in the 'Transfers' section
	debug 							= false,
	mgStats							= nil		-- any mgStats listed in the 'MGStats' section
}
\endcode

<br>
<h3>ILU</h3>
\code
{
	type 			= "ilu",
	beta 			= 0,
	damping 		= 1,
	sortEps 		= 1.e-50,
	inversionEps 	= 1.e-8
}
\endcode

<br>
<h3>ILUT</h3>
\code
{
	type 		= "ilut",
	threshold 	= 1e-6
}
\endcode

<br>
<h3>Jacobi</h3>
\code
{
	type 	= "jac",
	damping = 0.66
}
\endcode


<br>
<h3>Gauss Seidel</h3>
\code
{
	type = "gs"
}
\endcode

<br>
<h3>Block Gauss Seidel</h3>
\code
{
	type = "bgs"
}
\endcode

<br>
<h3>Symmetric Gauss Seidel</h3>
\code
{
	type = "sgs"
}
\endcode

<br>
<h3>Element Gauss Seidel</h3>
\code
{
	type = "egs"
}
\endcode

<br>
<h3>Schur</h3>
\code
{
	type 			= "schur",
	dirichletSolver	= "lu",
	skeletonSolver	= "lu"
}
\endcode


<br>
<h2>Convergence Checks</h2>
The following listing gives an overview over available convergence checks and their
default parameters.

<h3>Standard convergence check</h3>
\code
{
	type		= "standard",
	iterations	= 100,
	absolute	= 1e-12,
	reduction	= 1e-6,
	verbose		= true
}
\endcode


<br>
<h2>Transfer Operators</h2>
The following listing gives an overview over available transfer operators and their
default parameters.
<h3>Standard Transfer Operator</h3>
\code
{
	type				= "std",
	restrictionDamp 	= 1.0,
	prolongationDamp 	= 1.0,
	enableP1LagrangeOptimization = true,
	debug 				= false
}
\endcode


<br>
<h2>Line Search Methods</h2>
The following listing gives an overview over available line search methods and their
default parameters.

<h3>Standard Line Search Method</h3>
\code
{
	type			= "standard",
	maxSteps		= 10,
	lambdaStart		= 1,
	lambdaReduce	= 0.5,
	acceptBest 		= true,
	checkAll		= false
}
\endcode


<br>
<h2>MGStats</h2>
The following listing gives an overview over available MGStats objects.
MGStats objects are used to record statistics on individual multigrid cycles.
They add some overhead, so one should only use them for debugging.

<h3>MGStats</h3>
\code
{
	type			= "standard",
	filenamePrefix	= "mgstats",
	exitOnError		= false,
	writeErrVecs	= false,
	writeErrDiffs	= false,
	activeStages	= nil
}
\endcode

\{
]]--	


ug_load_script("util/table_desc_util.lua")
ug_load_script("util/table_util.lua")

util = util or {}
util.solver = util.solver or {}
util.solver.util = util.solver.util or {}

util.solver.defaults =
{
	approxSpace = nil,
	discretization = nil,

	nonlinearSolver =
	{
		newton = {
			convCheck	= "standard",
			linSolver	= "bicgstab",
			lineSearch	= nil
		}
	},
	
	linearSolver = 
	{
		linear = {
			precond		= "ilu",
			convCheck	= "standard"
		},

		cg = {
			precond		= "ilu",
			convCheck	= "standard"
		},

		bicgstab = {
			precond		= "ilu",
			convCheck	= "standard"
		},
	},

	preconditioner =
	{
		gmg = {
			adaptive = false,
			approxSpace = nil,
			baseLevel = 0,
			baseSolver = "lu",
			cycle = "V",
			discretization = nil,	-- only necessary if the underlying matrix is not of type AssembledLinearOperator
			gatheredBaseSolverIfAmbiguous = false,
			preSmooth = 3,
			postSmooth = 3,
			rap = false,
			rim = false,
			emulateFullRefined = false,
			smoother = "gs",
			transfer = "std",
			debug = false,
			mgStats = nil
		},

		ilu = {
			beta 			= 0,
			damping 		= 1,
			sort			= false,
			sortEps 		= 1.e-50,
			inversionEps 	= 1.e-8,
			consistentInterfaces = false,
			overlap 		= false
		},

		ilut = {
			threshold = 1e-6
		},
		
		gs = {
			consistentInterfaces = false
		},

		sgs = {
			consistentInterfaces = false
		},

		jac = {
			damping = 0.66
		},

		schur = {
			dirichletSolver	= "lu",
			skeletonSolver	= "lu"
		}
	},
	
	transfer = 
	{
		std = {
			restrictionDamp = 1.0,
			prolongationDamp = 1.0,
			enableP1LagrangeOptimization = true,
			debug = false
		}
	},

	convCheck =
	{
		standard = {
			iterations	= 100,		-- number of iterations
			absolute	= 1e-12,	-- absolute value of defect to be reached;
			reduction	= 1e-6,		-- reduction factor of defect to be reached;
			verbose		= true
		}
	},
	
	lineSearch =
	{
		standard = {
			maxSteps		= 10,
			lambdaStart		= 1,
			lambdaReduce	= 0.5,
			acceptBest 		= true,
			checkAll		= false
		}
	},

	mgStats =
	{
		standard = {
			filenamePrefix	= "mgstats",
			exitOnError		= false,
			writeErrVecs	= false,
    		writeErrDiffs	= false,
			activeStages	= nil
		}
	}
}


function util.solver.CondAbort(condition, message)
	if condition == true then
		print("ERROR in util.solver: " .. message)
		exit()
	end
end

--! if solverDesc is a table, solverDesc.solverutil will be set to solverutil or
--! initialized as an empty table. The method will then return solverDesc.solverutil.
--! If solverDesc is not a table, the method returns solverutil or {}, if
--! solverutil == nil.
function util.solver.PrepareSolverUtil(solverDesc, solverutil)
	if type(solverDesc) == "table" then
		solverDesc.solverutil = solverutil or solverDesc.solverutil or {}
		return solverDesc.solverutil
	else
		return solverutil or {}
	end
end

--! Creates a solver.
--! @param solverutil	You may OPTIONALLY pass a table solverutil in which
--! 					solver related information will be stored. If solverDesc
--!						is a table, you may also acacess this information through
--!						solverDesc.solverutil (even if the parameter solverutil == nil).
function util.solver.CreateSolver(solverDesc, solverutil)
	solverutil = util.solver.PrepareSolverUtil(solverDesc, solverutil)

	if util.tableDesc.IsPreset(solverDesc) then return solverDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(solverDesc)
	local defaults   = util.solver.defaults.nonlinearSolver[name]
	if desc == nil then desc = defaults end
	
	if name == "newton" then
		local newtonSolver = NewtonSolver()
		newtonSolver:set_linear_solver(
			util.solver.CreateLinearSolver(
				desc.linSolver or defaults.linSolver, solverutil))

		newtonSolver:set_convergence_check(
			util.solver.CreateConvCheck(
				desc.convCheck or defaults.convCheck, solverutil))
		
		local lineSearch = desc.lineSearch or defaults.lineSearch
		if lineSearch and lineSearch ~= "none" then
			newtonSolver:set_line_search(
				util.solver.CreateLineSearch(lineSearch))
		end
		
		util.solver.SetDebugWriter(newtonSolver, solverDesc, defaults)

		return newtonSolver
	else
		return util.solver.CreateLinearSolver(solverDesc, solverutil)
	end
end

function util.solver.CreateLinearSolver(solverDesc, solverutil)
	solverutil = util.solver.PrepareSolverUtil(solverDesc, solverutil)

	if util.tableDesc.IsPreset(solverDesc) then return solverDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(solverDesc)
	local defaults   = util.solver.defaults.linearSolver[name]
	if desc == nil then desc = defaults end


	-- create solver
	local linSolver = nil
	local createPrecond = false
	local createConvCheck = false

	if name == "linear" then
		linSolver = LinearSolver()
		createPrecond = true
		createConvCheck = true

	elseif name == "cg" then
		linSolver = CG()
		createPrecond = true
		createConvCheck = true

	elseif name == "bicgstab" then
		linSolver = BiCGStab()
		createPrecond = true
		createConvCheck = true

	elseif name == "lu"	then
		if HasClassGroup("SuperLU") then
			linSolver = AgglomeratingSolver(SuperLU());
		else
			linSolver = AgglomeratingSolver(LU())
		end

	elseif name == "uglu" then
		linSolver = AgglomeratingSolver(LU());
		
	elseif name == "superlu" then
		linSolver = AgglomeratingSolver(SuperLU());
	end

	util.solver.CondAbort(linSolver == nil, "Invalid linear solver specified: " .. name)
	
	if createPrecond == true then
		linSolver:set_preconditioner(
			util.solver.CreatePreconditioner(desc.precond or defaults.precond, solverutil))
	end

	if createConvCheck == true then
		linSolver:set_convergence_check(
			util.solver.CreateConvCheck(desc.convCheck or defaults.convCheck, solverutil))
	end
	
	util.solver.SetDebugWriter(linSolver, solverDesc, defaults)

	if desc then desc.instance = linSolver end
	return linSolver
end

function util.solver.CreateTransfer(transferDesc)

	local name, desc = util.tableDesc.ToNameAndDesc(transferDesc)
	local defaults   = util.solver.defaults.transfer[name]
	if desc == nil then desc = defaults end

	local transfer = nil
	if name == "std" then
		transfer = StdTransfer()
	end

	util.solver.CondAbort(transfer == nil, "Invalid transfer specified: " .. name)

	transfer:enable_p1_lagrange_optimization(desc.enableP1LagrangeOptimization or defaults.enableP1LagrangeOptimization)

	if type(desc.restrictionDamp) == "number" and desc.restrictionDamp ~= 1.0 then	
		transfer:set_restriction_damping(desc.restrictionDamp)
	end
	if type(desc.prolongationDamp) == "number" and desc.restrictionDamp ~= 1.0 then	
		transfer:set_prolongation_damping(desc.prolongationDamp)
	end
	
	util.solver.SetDebugWriter(transfer, transferDesc, defaults)
	
	return transfer
end

-- create a preconditioner
-- solverutil may be nil
function util.solver.CreatePreconditioner(precondDesc, solverutil)
	solverutil = util.solver.PrepareSolverUtil(solverDesc, solverutil)
	if util.tableDesc.IsPreset(precondDesc) then return precondDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(precondDesc)
	local defaults   = util.solver.defaults.preconditioner[name]
	if desc == nil then desc = defaults end

	local precond = nil

	local approxSpace = nil
	if desc then
		approxSpace = desc.approxSpace or util.solver.defaults.approxSpace
	end

	if name == "ilu"  then
		precond = ILU ()
		precond:set_beta (desc.beta or defaults.beta)
		precond:set_damp(desc.damping or defaults.damping)
		precond:set_sort(desc.sort or defaults.sort)
		precond:set_sort_eps(desc.sortEps or defaults.sortEps)
		precond:set_inversion_eps(desc.inversionEps or defaults.inversionEps)
		precond:enable_consistent_interfaces(desc.consistentInterfaces or defaults.consistentInterfaces)
		precond:enable_overlap(desc.overlap or defaults.overlap)
	elseif name == "ilut" then precond = ILUT (desc.threshold or defaults.threshold);
	elseif name == "jac"  then precond = Jacobi (desc.damping or defaults.damping);
	elseif name == "bgs"  then precond = BlockGaussSeidel ();
	elseif name == "gs"   then
		precond = GaussSeidel ()
		precond:enable_consistent_interfaces(desc.consistentInterfaces or defaults.consistentInterfaces)
	elseif name == "sgs"  then
		precond = SymmetricGaussSeidel ()
		precond:enable_consistent_interfaces(desc.consistentInterfaces or defaults.consistentInterfaces)
	elseif name == "egs"  then precond = ElementGaussSeidel ();

	elseif name == "gmg"  then 
		local smoother =
				util.solver.CreatePreconditioner(
					desc.smoother or defaults.smoother, solverutil)
		
		local baseSolver = 
				util.solver.CreateLinearSolver(
					desc.baseSolver or defaults.baseSolver, solverutil)

		if approxSpace == nil then
			print("An ApproximationSpace is required to create a 'gmg' solver.")
			print("Please specify one through the gmg-descriptor (member 'approxSpace')")
			print("or specify it through util.solver.defaults.approxSpace")
			exit()
		end

		local gmg = GeometricMultiGrid(approxSpace)
		gmg:set_base_solver					(baseSolver)
		gmg:set_smoother 					(smoother)
		gmg:set_base_level					(desc.baseLevel or defaults.baseLevel)
		gmg:set_cycle_type					(desc.cycle or defaults.cycle)
		gmg:set_num_presmooth				(desc.preSmooth or defaults.preSmooth)
		gmg:set_num_postsmooth				(desc.postSmooth or defaults.postSmooth)
		gmg:set_rap 						(desc.rap or defaults.rap)
		gmg:set_smooth_on_surface_rim		(desc.rim or defaults.rim)
		gmg:set_emulate_full_refined_grid	(desc.emulateFullRefined or defaults.emulateFullRefined)
		gmg:set_gathered_base_solver_if_ambiguous (
				desc.gatheredBaseSolverIfAmbiguous or
				defaults.gatheredBaseSolverIfAmbiguous)

		local transfer = util.solver.CreateTransfer(desc.transfer or defaults.transfer)
		if desc.adaptive == true then
		--	next three lines are obsolete, since gmg:set_transfer overwrites gmg:set_projection, anyways!?
			local project = StdTransfer()
			project:enable_p1_lagrange_optimization(false)
			gmg:set_projection(project)

			transfer:enable_p1_lagrange_optimization(false)
		end
		gmg:set_transfer(transfer)

		if (desc.debug or defaults.debug) then
			gmg:set_debug(GridFunctionDebugWriter(approxSpace))
		end

		local discretization = desc.discretization or util.solver.defaults.discretization
		if discretization then
			gmg:set_discretization(discretization)
		end

		local mgStats = util.solver.CreateMGStats(desc.mgStats or defaults.mgStats)
		if(mgStats) then
			gmg:set_mg_stats(mgStats)
		end

		precond = gmg

	elseif name == "schur" 	then
		local dirichletSolver =
				util.solver.CreateLinearSolver(
					desc.dirichletSolver or defaults.dirichletSolver, solverutil)

		local skeletonSolver =
				util.solver.CreateLinearSolver(
					desc.skeletonSolver or defaults.skeletonSolver, solverutil)

		local schur = SchurComplement()
		schur:set_dirichlet_solver(dirichletSolver)
		schur:set_skeleton_solver(SchurInverseWithFullMatrix(skeletonSolver))
		precond = schur
	end

	util.solver.CondAbort(precond == nil, "Invalid preconditioner specified: " .. name)
  
  

  -- check:   
  local rprecond  = precond
  if desc and desc.debugSolver then
    print(solver)
    local dsolver = 
      util.solver.CreateLinearSolver(
          desc.debugSolver, solverutil)
          
    local err = GridFunction(approxSpace)
    rprecond = DebugIterator(precond, dsolver)
    rprecond:set_solution(err)
  else  
  
  
   
  end
  
  -- create debug writer (optional)
  if desc and ((desc.debug == true) or (desc.debugSolver and desc.debugSolver.debug==true)) then
    if approxSpace == nil then
      print("An ApproximationSpace is required to create a DebugWriter for the '" .. name .. "'' preconditioner.")
      print("Consider setting the 'approxSpace' property of your preconditioner,")
      print("or alternatively the util.solver.defaults.approxSpace property or")
      print("alternatively set the option 'debug=false' in your preconditioner.")
      exit()
    end

    local dbgWriter = GridFunctionDebugWriter(approxSpace)
    dbgWriter:set_vtk_output(true)
    rprecond:set_debug(dbgWriter)
  end
  
  if desc then desc.instance = rprecond end
  
	return rprecond
end

-- solverutil may be nil
function util.solver.CreateConvCheck(convCheckDesc, solverutil)
	solverutil = util.solver.PrepareSolverUtil(convCheckDesc, solverutil)
	if util.tableDesc.IsPreset(convCheckDesc) then return convCheckDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(convCheckDesc)
	local defaults	 = util.solver.defaults.convCheck[name]
	if desc == nil then desc = defaults end
	
	local cc = nil
	
	local verbose = true;
	if desc.verbose ~= nil then verbose = desc.verbose
	elseif defaults.verbose ~= nil then verbose = defaults.verbose end

	if name == "standard" then
		cc = ConvCheck()
		cc:set_maximum_steps	(desc.iterations	or defaults.iterations)
		cc:set_minimum_defect	(desc.absolute		or defaults.absolute)
		cc:set_reduction		(desc.reduction		or defaults.reduction)
		cc:set_verbose			(verbose)
	elseif name == "composite" then
		local approxSpace = desc.approxSpace or util.solver.defaults.approxSpace
		if approxSpace == nil then 
			print("For a composite convergence check the approximation space must be passed.")
			exit()
		end
		cc = CompositeConvCheck(approxSpace, 
			desc.iterations	or defaults.iterations, 
			desc.absolute		or defaults.absolute,
			desc.reduction		or defaults.reduction)
		cc:set_verbose			(verbose)
		if desc.sub and type(desc.sub) == "table" then
			for _, v in pairs(desc.sub) do
				cc:set_component_check(v.cmp, 
					(v.absolute or desc.absolute or defaults.absolute), 
					(v.relative or desc.reduction or defaults.reduction))
			end
		end
	end

	util.solver.CondAbort(cc == nil, "Invalid conv-check specified: " .. name)
	if desc then
		desc.instance = cc
		solverutil.convCheckDescs = solverutil.convCheckDescs or {}
		table.insert(solverutil.convCheckDescs, desc)
	end

	return cc
end

function util.solver.CreateLineSearch(lineSearchDesc)
	if util.tableDesc.IsPreset(lineSearchDesc) then return lineSearchDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(lineSearchDesc)
	local defaults	 = util.solver.defaults.lineSearch[name]
	if desc == nil then desc = defaults end
	
	local ls = nil
	if name == "standard" then
	--	battle booleans
		if desc.acceptBest == nil then desc.acceptBest = defaults.acceptBest end
		if desc.checkAll == nil then desc.checkAll = defaults.checkAll end
		
		ls = StandardLineSearch(desc.maxSteps or defaults.maxSteps,
								desc.lambdaStart or defaults.lambdaStart,
								desc.lambdaReduce or defaults.lambdaReduce,
								desc.acceptBest,
								desc.checkAll)
		if desc.verbose ~= nil then ls:set_verbose(desc.verbose) end
		if desc.suffDesc ~= nil then ls:set_suff_descent_factor(desc.suffDesc) end
		if desc.maxDefect ~= nil then ls:set_maximum_defect(desc.maxDefect) end
	end
	
	util.solver.CondAbort(ls == nil, "Invalid line-search specified: " .. name)
	if desc then
		desc.instance = ls
	end

	return ls
end


function util.solver.CreateMGStats(mgStatsDesc)
	if mgStatsDesc == nil then return nil end
	
	if util.tableDesc.IsPreset(mgStatsDesc) then return mgStatsDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(mgStatsDesc)
	local defaults	 = util.solver.defaults.mgStats[name]
	if desc == nil then desc = defaults end
	
	local mgStats = nil
	if name == "standard" then
		if desc.exitOnError == nil then desc.exitOnError = defaults.exitOnError end
		if desc.writeErrVecs == nil then desc.writeErrVecs = defaults.writeErrVecs end
		if desc.writeErrDiffs == nil then desc.writeErrDiffs = defaults.writeErrDiffs end

		mgStats = MGStats()
		mgStats:set_filename_prefix(desc.filenamePrefix or defaults.filenamePrefix)
		mgStats:set_exit_on_error(desc.exitOnError)
		mgStats:set_write_err_vecs(desc.writeErrVecs)
		mgStats:set_write_err_diffs(desc.writeErrDiffs)
		if desc.activeStages or defaults.activeStages then
			mgStats:set_active_stages(desc.activeStages or defaults.activeStages)
		end
	end
	
	util.solver.CondAbort(mgStats == nil, "Invalid mgStats specified: " .. name)
	if desc then
		desc.instance = mgStats
	end

	return mgStats
end

function util.solver.SetDebugWriter(obj, desc, defaults)
	
	local dbgDesc = desc.debug
	if dbgDesc == nil then
		if defaults then dbgDesc = defaults.debug end
	end
	
	if dbgDesc then
		local debug = false
		local vtk = true
		local conn_viewer = false
	
		if type (dbgDesc) == "boolean" then
			debug = dbgDesc
		elseif type (dbgDesc) == "table" then
			debug = true
			if dbgDesc.vtk ~= nil then vtk = dbgDesc.vtk end
			if dbgDesc.conn_viewer ~= nil then conn_viewer = dbgDesc.conn_viewer end
		else
			print("An unrecognized type of the debug writer specification.")
			exit()
		end
	
		if debug then
	
			local approxSpace = nil
			if desc then
				approxSpace = desc.approxSpace or util.solver.defaults.approxSpace
			end
	
			if approxSpace == nil then
			  print("An ApproximationSpace is required to create a DebugWriter.")
			  print("Consider setting the 'approxSpace' property of of the object to debug,")
			  print("or alternatively the util.solver.defaults.approxSpace property.")
			  print("Otherwiese set the option 'debug=false' for the object.")
			  exit()
			end
	
			local dbgWriter = GridFunctionDebugWriter(approxSpace)
			dbgWriter:set_vtk_output(vtk)
			dbgWriter:set_conn_viewer_output(conn_viewer)
			obj:set_debug(dbgWriter)
		
		end
	end
end


--! Prepares a solution step, e.g. during nested iterations.
--! @param desc			a solver-desc which was used during solver-creation in one of
--! 					the solver_util methods or a solverutil table, which is was 
--! 					also created in one of the solver_util methods or simply
--!						a table which holds a bunch of Convergence-Check descriptor tables.
--!
--! @param nestedStep	(optional) the step of the nested iteration (default: 1)
--!
--! @param timeStep		(optional) the current time step (default: 1)
function util.solver.PrepareStep(desc, nestedStep, timeStep)
	if nestedStep == nil then nestedStep = 1 end
	if timeStep == nil then timeStep = 1 end

	local convCheckDescs = (desc.solverutil and desc.solverutil.convCheckDescs)
						 or desc.convCheckDescs
						 or desc

	if convCheckDescs then
		for _, ccDesc in pairs(convCheckDescs) do
			local cc = ccDesc.instance
			for _, d in pairs(ccDesc) do
				if type(d) == "table" then
					if	(d.nestedStep == nil or d.nestedStep == nestedStep) and
						(d.timeStep == nil or d.timeStep == timeStep)
					then
						if d.iterations		then cc:set_maximum_steps(d.iterations) end
						if d.absolute		then cc:set_minimum_defect(d.absolute) end
						if d.reduction		then cc:set_reduction(d.reduction) end
						if d.verbose ~= nil	then cc:set_verbose(d.verbose) end
					end
				end
			end
		end
	end
end
