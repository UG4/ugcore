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
--					absolut		= 1e-12,
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

util = util or {}
util.solver = util.solver or {}
util.solver.util = util.solver.util or {}

util.solver.defaults =
{
	approxSpace = nil,

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
			beta = 0
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
			absolut		= 1e-12,	-- absolut value of defect to be reached;
			reduction	= 1e-6,		-- reduction factor of defect to be reached;
			verbose		= true
		}
	}
}


function util.solver.util.CondAbort(condition, message)
	if condition == true then
		print(message)
		exit()
	end
end

function util.solver.util.IsPreset(desc)
	if type(desc) == "userdata" then
		return true
	else
		return false
	end
end

function util.solver.util.ToNameAndDesc(descOrName)
	if type(descOrName) == "string" then
		return descOrName, nil
	elseif type(descOrName) == "table" then
		return descOrName.type, descOrName
	end
	print("Invalid name or descriptor specified!")
	exit()
end


function util.solver.CreateSolver(solverDesc)
--	todo: add non-linear solvers
	return util.solver.CreateLinearSolver(solverDesc)
end


function util.solver.CreateLinearSolver(solverDesc)
	if util.solver.util.IsPreset(solverDesc) then return solverDesc end

	local name, desc = util.solver.util.ToNameAndDesc(solverDesc)
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

	util.solver.util.CondAbort(linSolver == nil, "Invalid linear solver specified: " .. name)
	
	if createPrecond == true then
		linSolver:set_preconditioner(
			util.solver.CreatePreconditioner(desc.precond or defaults.precond))
	end

	if createConvCheck == true then
		linSolver:set_convergence_check(
			util.solver.CreateConvCheck(desc.convCheck or defaults.convCheck))
	end

	if desc then desc.instance = linSolver end
	return linSolver
end


function util.solver.CreatePreconditioner(precondDesc)
	if util.solver.util.IsPreset(precondDesc) then return precondDesc end

	local name, desc = util.solver.util.ToNameAndDesc(precondDesc)
	local defaults   = util.solver.defaults.preconditioner[name]
	if desc == nil then desc = defaults end

	local precond = nil

	if name == "ilu"  then
		precond = ILU ()
		precond:set_beta (desc.beta or defaults.beta)
	
	elseif name == "ilut" then precond = ILUT (desc.threshold or defaults.threshold);
	elseif name == "jac"  then precond = Jacobi (desc.damping or defaults.damping);
	elseif name == "bgs"  then precond = BlockGaussSeidel ();
	elseif name == "gs"   then precond = GaussSeidel ();
	elseif name == "sgs"  then precond = SymmetricGaussSeidel ();
	elseif name == "egs"  then precond = ElementGaussSeidel ();

	elseif name == "gmg"  then 
		local smoother =
				util.solver.CreatePreconditioner(
					desc.smoother or defaults.smoother)
		
		local baseSolver = 
				util.solver.CreateLinearSolver(
					desc.baseSolver or defaults.baseSolver)

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

		if desc.discretization then
			gmg:set_discretization(desc.discretization)
		end

		precond = gmg

	elseif name == "schur" 	then
		local dirichletSolver =
				solver.util.CreateLinearSolver(
					desc.dirichletSolver or defaults.dirichletSolver)

		local skeletonSolver =
				solver.util.CreateLinearSolver(
					desc.skeletonSolver or defaults.skeletonSolver)

		local schur = SchurComplement()
		schur:set_dirichlet_solver(dirichletSolver)
		schur:set_skeleton_solver(SchurInverseWithFullMatrix(skeletonSolver))
		precond = schur
	end

	util.solver.util.CondAbort(precond == nil, "Invalid preconditioner specified: " .. name)
	if desc then desc.instance = precond end

	return precond
end


function util.solver.CreateConvCheck(convCheckDesc)
	if util.solver.util.IsPreset(convCheckDesc) then return convCheckDesc end

	local name, desc = util.solver.util.ToNameAndDesc(convCheckDesc)
	local defaults	 = util.solver.defaults.convCheck[name]
	if desc == nil then desc = defaults end
	
	local cc = nil
	
	local verbose = true;
	if desc.verbose ~= nil then verbose = desc.verbose
	elseif defaults.verbose ~= nil then verbose = defaults.verbose end

	if name == "standard" then
		cc = ConvCheck()
		cc:set_maximum_steps	(desc.iterations	or defaults.iterations)
		cc:set_minimum_defect	(desc.absolut		or defaults.absolut)
		cc:set_reduction		(desc.reduction		or defaults.reduction)
		cc:set_verbose			(verbose)
	end

	util.solver.util.CondAbort(cc == nil, "Invalid conv-check specified: " .. name)
	if desc then desc.instance = cc end
	return cc
end
