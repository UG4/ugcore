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
-- A not necessarily meaningful example follows:
-- (here someApproxSpace points to a previously created approximation space)
--
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

	iterativeSolver =
	{
		ilut = {
			threshold = 1e-6
		},

		jacobi = {
			damping = 0.66
		}
	},

	preconditioner =
	{
		gmg = {
			adaptive = true,
			baseLevel = 0,
			baseSolver = "lu",
			cycle = "V",
			smoother = "gs",
			preSmooth = 3,
			postSmooth = 3,
			rap = false,
			gatheredBaseSolverIfAmbiguous = false,
			approxSpace = nil,
			discretization = nil
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
	return util.solver.CreateLinearSolver(solverDesc, true)
end


function util.solver.CreateLinearSolver(solverDesc, agglomerating)
	if agglomerating == nil then agglomerating = true end

	local name, desc = util.solver.util.ToNameAndDesc(solverDesc)
	
	local defaults = util.solver.defaults.linearSolver[name]
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
		if agglomerating then
			linSolver = AgglomeratingSolver(LU())
		else
			linSolver = LU()
		end
	elseif name == "superlu"	then
		if agglomerating then
			linSolver = AgglomeratingSolver(SuperLU());
		else
			linSolver = SuperLU()
		end

	elseif name == "schur" 	then
		-- todo: Move "schur" to preconditioners?
		local skeletonSolver = AgglomeratingSolver(SuperLU())
		local schur = SchurComplement()
		schur:set_dirichlet_solver(SuperLU())
		schur:set_skeleton_solver(SchurInverseWithFullMatrix(skeletonSolver))
		linSolver = LinearSolver()
		linSolver:set_preconditioner(schur)

	else 
		print("Linear solver '"..name.."' not found."); 
		exit(); 
	end
	
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
	local name, desc = util.solver.util.ToNameAndDesc(precondDesc)

	local precond = nil
	if name == "gmg" then 
		local defaults = util.solver.defaults.preconditioner.gmg
		if desc == nil then desc = defaults end

		local smoother =
				util.solver.CreateIterativeSolver(
					desc.smoother or defaults.smoother)
		
		local baseSolver = 
				util.solver.CreateLinearSolver(
					desc.baseSolver or defaults.baseSolver, false)

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
	else
		precond = util.solver.CreateIterativeSolver(precondDesc)
	end

	util.solver.util.CondAbort(precond == nil, "Invalid preconditioner specified: " .. name)
	if desc then desc.instance = precond end

	return precond
end


function util.solver.CreateIterativeSolver(iterSolverDesc)
	local name, desc = util.solver.util.ToNameAndDesc(iterSolverDesc)
	local defaults = util.solver.defaults.iterativeSolver[name]
	if desc == nil then desc = defaults end

	local solver = nil
	if 	    name == "ilu"  then solver = ILU ();
	elseif 	name == "ilut" then solver = ILUT (desc.threshold or defaults.threshold);
	elseif 	name == "jac"  then solver = Jacobi (desc.damping or defaults.damping);
	elseif 	name == "bgs"  then solver = BlockGaussSeidel ();
	elseif 	name == "gs"   then solver = GaussSeidel ();
	elseif 	name == "sgs"  then solver = SymmetricGaussSeidel ();
	elseif  name == "egs"  then solver = ElementGaussSeidel ();
	else
		print("Iterative solver " .. " '" .. name .."' not found.")
		exit()
	end

	if desc then desc.instance = solver end
	return solver
end


function util.solver.CreateConvCheck(convCheckDesc)
	
	local name, desc = util.solver.util.ToNameAndDesc(convCheckDesc)
	local cc = nil
	
	defaults = util.solver.defaults.convCheck[name]
	if desc == nil then desc = defaults end
	
	if name == "standard" then
		cc = ConvCheck()
		cc:set_maximum_steps	(desc.iterations	or defaults.iterations)
		cc:set_minimum_defect	(desc.absolut		or defaults.absolut)
		cc:set_reduction		(desc.reduction		or defaults.reduction)
		cc:set_verbose			(desc.verbose		or defaults.verbose)
	end

	util.solver.util.CondAbort(cc == nil, "Invalid conv-check specified: " .. name)
	if desc then desc.instance = cc end
	return cc
end
