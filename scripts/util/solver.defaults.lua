return {
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
			precond	= "ilu",
			convCheck	= "standard"
		},

		cg = {
			precond	= "ilu",
			convCheck	= "standard"
		},

		bicgstab = {
			precond	= "ilu",
			convCheck	= "standard"
		},
		
		lu = {
			info = false,
			showProgress = true
		},
		
		uglu = {
			info = false,
			showProgress = true
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
			inversionEps 		= 1.e-8,
			consistentInterfaces   = false,
			overlap 		= false,
			ordering 		= nil
		},

		ilut = {
			threshold = 1e-6,
			ordering = "NativeCuthillMcKee"
		},
		
		gs = {
			consistentInterfaces 	= false,
			overlap 	 	= false
		},

		sgs = {
			consistentInterfaces	= false,
			overlap		= false
		},

		jac = {
			damping = 0.66
		},

		schur = {
			dirichletSolver	= "lu",
			skeletonSolver		= "lu"
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
			reduction	= 1e-6,	-- reduction factor of defect to be reached;
			verbose	= true
		}
	},
	
	lineSearch =
	{
		standard = {
			maxSteps		= 10,
			lambdaStart		= 1,
			lambdaReduce		= 0.5,
			acceptBest 		= true,
			checkAll		= false
		}
	},

	mgStats =
	{
		standard = {
			filenamePrefix	= "mgstats",
			exitOnError	= false,
			writeErrVecs	= false,
			writeErrDiffs	= false,
			activeStages	= nil
		}
	},

	timeSolver = 
	{
		fixed = {
			scheme = "ImplEuler",
			dt = 0.1,
			start = 0,
			stop = 1
		}
	}
}