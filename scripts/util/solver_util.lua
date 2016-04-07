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


-- todo: enhance 'CreateSolver', which can also create a NewtonSolver and similar
--		 if requested.

-- Given a (possibly nested) solver description, e.g.,
-- 
-- 	solver = {
--		type = "bicgstab",
--		precond = "ilu"
--		convCheck = "standard"
-- 	}
--
-- a call e.g. to util.solver.CreateSolver(solver)
-- creates and returns the requested solver.
--
-- For conveniance, all the methods will add an 'instance' entry to their respective
-- descriptor tables, assigning the created solver component.
-- If you thus setup a solver, e.g.
--
-- 	solver = {
--		type = "bicgstab",
--		precond = {
--			type = "ilu"
--		}
-- 	}
--
-- you may later access 'solver.instance' and 'solver.precond.instance'.
-- Note that this is only possible if you specified a component through a descriptor.
-- If you use the short form, e.g.
--
-- solver = {
--		type = "bicgstab",
--		precond = "ilu"
-- 	}
--
-- you will still be able to access 'solver.instance', but you can't access
-- 'solver.precond.instance'.
--
--
-- All components may either be specified through a string-id or through a
-- descriptor. You only have to specify those components which deviate from
-- the defaults (solver.util.defaults), the rest will be ammended.
--
--
-- If you want to use a geometric multigrid preconditioner, you either have
-- to specify an approxSpace entry in its descriptor or set
-- 'util.solver.defaults.approxSpace' before calling any 'util.solver.Create...'
-- method.
--
--
-- Two not necessarily meaningful examples follow:
-- (here someApproxSpace points to a previously created approximation space)
--
--  -- EXAMPLE 1
--	solverDesc = 
--	{
--		type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
--		precond = 
--		{	
--			type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
--			smoother 	= "gs",		-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
--			cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
--			preSmooth	= 3,		-- number presmoothing steps
--			postSmooth 	= 3,		-- number postsmoothing steps
--			rap			= false,	-- comutes RAP-product instead of assembling if true 
--			baseLevel	= 0,		-- gmg - baselevel
--			baseSolver	= {			-- better options are most likely "lu" or "superlu"
--				type	  = "bicgstab",
--				precond	  = "gs",
--				convCheck = {
--					type		= "standard"
--					iterations	= 1000,
--					absolute	= 1e-12,
--					reduction	= 1e-10,
--					verbose		= false
--				}
--			}
--			approxSpace	= someApproxSpace,
--		},
--	
--		convCheck = "standard"
--	}
--
--	solver = util.solver.CreateSolver(solverDesc)
--
--
--  -- EXAMPLE 2
--	solverDesc =
--	{
--		type = "linear",
--		precond = ILU()		-- for each component one can optionally supply a custom instance
--	}
--
--	solver = util.solver.CreateSolver(solverDesc)

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
			smoother = "gs"
		},

		ilu = {
			beta 			= 0,
			damping 		= 1,
			sortEps 		= 1.e-50,
			inversionEps 	= 1.e-8
		},

		ilut = {
			threshold = 1e-6
		},

		jac = {
			damping = 0.66
		},

		schur = {
			dirichletSolver	= "lu",
			skeletonSolver	= "lu"
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
				util.solver.CreateLineSearch(lineSearch, solverutil))
		end

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

	if desc then desc.instance = linSolver end
	return linSolver
end


function util.solver.CreatePreconditioner(precondDesc, solverutil)
	solverutil = util.solver.PrepareSolverUtil(solverDesc, solverutil)
	if util.tableDesc.IsPreset(precondDesc) then return precondDesc end

	local name, desc = util.tableDesc.ToNameAndDesc(precondDesc)
	local defaults   = util.solver.defaults.preconditioner[name]
	if desc == nil then desc = defaults end

	local precond = nil

	if name == "ilu"  then
		precond = ILU ()
		precond:set_beta (desc.beta or defaults.beta)
		precond:set_damp(desc.damping or defaults.damping)
		precond:set_sort_eps(desc.sortEps or defaults.sortEps)
		precond:set_inversion_eps(desc.inversionEps or defaults.inversionEps)
	elseif name == "ilut" then precond = ILUT (desc.threshold or defaults.threshold);
	elseif name == "jac"  then precond = Jacobi (desc.damping or defaults.damping);
	elseif name == "bgs"  then precond = BlockGaussSeidel ();
	elseif name == "gs"   then precond = GaussSeidel ();
	elseif name == "sgs"  then precond = SymmetricGaussSeidel ();
	elseif name == "egs"  then precond = ElementGaussSeidel ();

	elseif name == "gmg"  then 
		local smoother =
				util.solver.CreatePreconditioner(
					desc.smoother or defaults.smoother, solverutil)
		
		local baseSolver = 
				util.solver.CreateLinearSolver(
					desc.baseSolver or defaults.baseSolver, solverutil)

		local approxSpace = desc.approxSpace or util.solver.defaults.approxSpace
		if approxSpace == nil then
			print("An ApproximationSpace is required to create a 'gmg' solver.")
			print("Please specify one through the gmg-descriptor (member 'approxSpace')")
			print("or specify it through util.solver.defaults.approxSpace")
			exit()
		end

		local gmg = GeometricMultiGrid(approxSpace)
		gmg:set_base_solver		(baseSolver)
		gmg:set_smoother 		(smoother)
		gmg:set_base_level		(desc.baseLevel or defaults.baseLevel)
		gmg:set_cycle_type		(desc.cycle or defaults.cycle)
		gmg:set_num_presmooth	(desc.preSmooth or defaults.preSmooth)
		gmg:set_num_postsmooth	(desc.postSmooth or defaults.postSmooth)
		gmg:set_rap 			(desc.rap or defaults.rap)
		gmg:set_gathered_base_solver_if_ambiguous (
				desc.gatheredBaseSolverIfAmbiguous or
				defaults.gatheredBaseSolverIfAmbiguous)

		if desc.adaptive == true then
			local transfer = StdTransfer()
			transfer:enable_p1_lagrange_optimization(false)
			gmg:set_transfer(transfer)
		end

		local discretization = desc.discretization or util.solver.defaults.discretization
		if discretization then
			gmg:set_discretization(discretization)
		end

		precond = gmg

	elseif name == "schur" 	then
		local dirichletSolver =
				solver.util.CreateLinearSolver(
					desc.dirichletSolver or defaults.dirichletSolver, solverutil)

		local skeletonSolver =
				solver.util.CreateLinearSolver(
					desc.skeletonSolver or defaults.skeletonSolver, solverutil)

		local schur = SchurComplement()
		schur:set_dirichlet_solver(dirichletSolver)
		schur:set_skeleton_solver(SchurInverseWithFullMatrix(skeletonSolver))
		precond = schur
	end

	util.solver.CondAbort(precond == nil, "Invalid preconditioner specified: " .. name)
	if desc then desc.instance = precond end

	return precond
end


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


function util.solver.CreateLineSearch(lineSearchDesc, solverutil)
	solverutil = util.solver.PrepareSolverUtil(lineSearchDesc, solverutil)
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
	end
	
	util.solver.CondAbort(ls == nil, "Invalid line-search specified: " .. name)
	if desc then
		desc.instance = ls
	end

	return ls
end

--! Prepares a solution step, e.g. during nested iterations.
--! @param desc			a solver-desc which was used during solver-creation in one of
--! 					the solver_util_2 methods or a solverutil table, which is was 
--! 					also created in one of the solver_util_2 methods or simply
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
