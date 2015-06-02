-- created by Sebastian Reiter and Andreas Vogel

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
--			basesolver	= {			-- better options are most likely "lu" or "superlu"
--				type	  = "bicgstab",
--				precond	  = "gs",
--				convCheck = {
--					type		= "standard"
--					iteration	= 1000,
--					absolute		= 1e-12,
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

ug_load_script("table_desc_util.lua")

util = util or {}
util.solver = util.solver or {}
util.solver.util = util.solver.util or {}

util.solver.defaults =
{
	approxSpace = nil,
	discretization = nil,

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
			beta = 0,
			damping = 1
		},

		ilut = {
			threshold = 1e-6
		},

		jac = {
			damping = 0.66
		},

		schur = {
			dirichletSolver	= "superlu",
			skeletonSolver	= "superlu"
		}
	},

	convCheck =
	{
		standard = {
			iterations	= 100,		-- number of iterations
			absolute		= 1e-12,	-- absolute value of defect to be reached;
			reduction	= 1e-6,		-- reduction factor of defect to be reached;
			verbose		= true
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
--	todo: add non-linear solvers
	return util.solver.CreateLinearSolver(solverDesc, solverutil)
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
		linSover = CG()
		createPrecond = true
		createConvCheck = true

	elseif name == "bicgstab" then
		linSolver = BiCGStab()
		createPrecond = true
		createConvCheck = true

	elseif name == "lu"		then
		linSolver = AgglomeratingSolver(LU())

	elseif name == "superlu"	then
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
	solverutil = util.solver.PrepareSolverUtil(solverDesc, solverutil)
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
	end

	util.solver.CondAbort(cc == nil, "Invalid conv-check specified: " .. name)
	if desc then
		desc.instance = cc
		solverutil.convCheckDescs = solverutil.convCheckDescs or {}
		table.insert(solverutil.convCheckDescs, desc)
	end

	return cc
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
