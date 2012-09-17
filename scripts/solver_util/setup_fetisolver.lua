----------------------------------------------------------
--[[!
--   \file setup_fetisolver.lua
--   \brief Lua - Script which creates and configures a FETI solver (linear solver).
-- 	 \author Ingo Heppner (using FETI solver specific code from 'fetitest.lua').--
--
--   Begin: 06072011.
--   Provides: Function which set up a FETI solver object.
--
--          Input parameters:
--          - 'str_problem' -- string describing problem (not yet used)
--          - 'dim',
--          - 'linMaxIterations'
--          - 'numProcs'
--          - 'u'
--          - 'dirichletBND', 'approxSpace',
--          - 'activateDbgWriter', 'verbosity', logfileName.
--          Return value:
--          Reference of the fully configured FETI solver object.
--
--   Note: To use FAMG as sub problem solvers configure ug4 as follows (since FAMG
--   is now contained in the 'amg' plugin):
--   cmake <other cmake options> -Damg=ON ..
--
--   Description of some parameters / options:
--      -AMGwriteMat:   write testvectors (only if Neumann or Dirichlet problem solver is of type "famg")

--   Usage example:
--	1. In the calling LUA script:
--
--  \code
--	ug_load_script("setup_fetisolver.lua")
--	solver = SetupFETISolver(str_problem,
--				 dim
--				 linMaxIterations,
--				 numProcs,
--				 u,
--				 dirichletBND, approxSpace,
--				 activateDbgWriter,
--				 verbosity, logfileName)
--  \endcode
--      where 'solver' is the solver object which is called by
--      'ApplyLinearSolver(., ., ., <solver>)', i.e.:
--	
--  \code
--	print("Apply solver.")
--	if ApplyLinearSolver(linOp, u, b, solver) == false then
--		print("Could not apply linear solver.");
--	end
--  \endcode
--
--	2. Execution of the calling LUA script:

\verbatim
UGARGS="-ex ../apps/scaling_tests/modular_scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti"

# Local:
########
openmpirun -np 4 ugshell $UGARGS -numPPSD 1 -numRefs 5 # etc.
openmpirun -np 4 ugshell $UGARGS -numPPSD 1 -numRefs 5 -dbgw 1 # with debug writer
openmpirun -np 4 ugshell $UGARGS -numPPSD 1 -numRefs 5 -logtofile bla -rlf # with logfile and automatic naming
	
# cekon:
########
salloc -n  4 mpirun ./ugshell $UGARGS -numPPSD 1 -numRefs 5
\endverbatim
!]]--
--[[
# JuGene:
#########
# Interactive - Job lief jedoch nicht mehr (01092011):
llrun -v  -np 16 -exe ./ugshell -mode VN -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ -args $UGARGS -numPPSD 1 -outproc 0

# Batch:
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -numPPSD 1"
mpirun -np  64 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 8 -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs8_pe64.txt"

mpirun -np 256 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 8 -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs8_pe256.txt -verb 1"   - geht!

mpirun -np 256 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 9 -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs9_pe256.txt" - hier geht "step 3.b" auf folgenden Procs nicht:
95,253,1,117,215,213,127,245,247,223,221,64,65,69,20,21,16,17,4,5
Umgeordnet:
echo "95,253,1,117,215,213,127,245,247,223,221,64,65,69,20,21,16,17,4,5" | tr "," "\012" | sort -n | tr "\012" ","
1,4,5,16,17,20,21,64,65,69,95,117,127,213,215,221,223,245,247,253 - 68, das ebenfalls am linken Rand liegt, taucht nicht auf; 65 liegt am oberen Rand
Automatisch, direkt aus dem Output-File auf JuGene:
grep "Could not solve Dirichlet problem " ug4_laplace_feti.204348.out_feti-sd1_8x8-quad_prerefs3-refs9_bisect_pe256 | awk ' { printf( "%d\n", $2) }' | sort -n | tr "\012" ","
1,4,5,16,17,20,21,64,65,69,95,117,127,213,215,221,223,245,247,253

# Mit '-distType grid2d':
mpirun -np 256 -exe ./ugshell -mode VN -mapfile TXYZ -verbose 2 -env LD_LIBRARY_PATH=/bgsys/drivers/ppcfloor/comm/lib/ \
                       -args "$UGARGS -numRefs 9 -distType grid2d -logtofile jugene_feti-sd1_8x8-quad_prerefs3-refs9_grid2d_pe256.txt -verb 1"
79,176,239,16,111,159,144,96,127,128,191,112,64,143,48,207,32,223,175,80,12,243,11,8,254,253,252,1,3,2,247,244,5,4,7,6,249,250,251,248
Umgeordnet:
echo "79,176,239,16,111,159,144,96,127,128,191,112,64,143,48,207,32,223,175,80,12,243,11,8,254,253,252,1,3,2,247,244,5,4,7,6,249,250,251,248" | tr "," "\012" | sort -n | tr "\012" ","
1,2,3,4,5,6,7,8,11,12,16,32,48,64,79,80,96,111,112,127,128,143,144,159,175,176,191,207,223,239,243,244,247,248,249,250,251,252,253,254
Automatisch, direkt aus dem Output-File auf JuGene:
grep "Could not solve Dirichlet problem " ug4_laplace_feti.204720.out_feti-sd1_8x8-quad_prerefs3-refs9_grid2d_pe256 | awk ' { printf( "%d\n", $2) }' | sort -n | tr "\012" ","
1,2,3,4,5,6,7,8,11,12,16,32,48,64,79,80,96,111,112,127,128,143,144,159,175,176,191,207,223,239,243,244,247,248,249,250,251,252,253,254

  Number of DoF's, starting from coarse grid 8 x 8 quadrilaterals:
  Level  0 ==> (2^3 * 2^{ 0} + 1)^2 nodes = (2^{ 3} + 1)^2 nodes =              81 nodes
  Level  1 ==> (2^3 * 2^{ 1} + 1)^2 nodes = (2^{ 4} + 1)^2 nodes =             289 nodes
  Level  2 ==> (2^3 * 2^{ 2} + 1)^2 nodes = (2^{ 5} + 1)^2 nodes =           1'089 nodes
  Level  3 ==> (2^3 * 2^{ 3} + 1)^2 nodes = (2^{ 6} + 1)^2 nodes =           4'225 nodes
  Level  4 ==> (2^3 * 2^{ 4} + 1)^2 nodes = (2^{ 7} + 1)^2 nodes =          16'641 nodes
  Level  5 ==> (2^3 * 2^{ 5} + 1)^2 nodes = (2^{ 8} + 1)^2 nodes =          66'049 nodes
  Level  6 ==> (2^3 * 2^{ 6} + 1)^2 nodes = (2^{ 9} + 1)^2 nodes =         263'169 nodes
  Level  7 ==> (2^3 * 2^{ 7} + 1)^2 nodes = (2^{10} + 1)^2 nodes =       1'050'625 nodes - das war bei Klawonn & Rheinbach (2010) die Problemgroesse fuer 16 PE ...
  Level  8 ==> (2^3 * 2^{ 8} + 1)^2 nodes = (2^{11} + 1)^2 nodes =       4'198'401 nodes
  Level  9 ==> (2^3 * 2^{ 9} + 1)^2 nodes = (2^{12} + 1)^2 nodes =      16'785'409 nodes
  Level 10 ==> (2^3 * 2^{10} + 1)^2 nodes = (2^{13} + 1)^2 nodes =      67'125'249 nodes
  Level 11 ==> (2^3 * 2^{11} + 1)^2 nodes = (2^{14} + 1)^2 nodes =     268'468'225 nodes
  Level 12 ==> (2^3 * 2^{12} + 1)^2 nodes = (2^{15} + 1)^2 nodes =   1'073'807'361 nodes
  Level 13 ==> (2^3 * 2^{13} + 1)^2 nodes = (2^{16} + 1)^2 nodes =   4'295'098'369 nodes
  Level 14 ==> (2^3 * 2^{14} + 1)^2 nodes = (2^{17} + 1)^2 nodes =  17'180'131'329 nodes
  Level 15 ==> (2^3 * 2^{15} + 1)^2 nodes = (2^{18} + 1)^2 nodes =  68'720'001'025 nodes
!]]--


----------------------------------------------------------
-- auxiliary functions for FAMG
-- Testvectors for FAMG ---
----------------------------------------------------------
function ourTestvector2d_0_0(x, y, t)
	return 0
end

function ourTestvector2d_1_1(x, y, t)
	return math.sin(math.pi*x)*math.sin(math.pi*y)
end

function ourTestvector2d_2_1(x, y, t)
	return math.sin(2*math.pi*x)*math.sin(math.pi*y)
end


function ourTestvector2d_1_2(x, y, t)
	return math.sin(math.pi*x)*math.sin(2*math.pi*y)
end


function ourTestvector2d_2_2(x, y, t)
	return math.sin(2*math.pi*x)*math.sin(2*math.pi*y)
end

function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
	local amgTestvector;
	print("          Create writer for testvector via grid function for FAMG ...")
	amgTestvector = GridFunctionVectorWriter()
	amgTestvector:set_reference_grid_function(gridfunction)
	amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
	return amgTestvector	
end

function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	local amgDirichlet0
	print("          Create vector writer for 'Dirichlet 0, constant 1 else' testvector for FAMG ...")
	amgDirichlet0 = GridFunctionVectorWriterDirichlet0()
	amgDirichlet0:init(dirichletBND, approxSpace)
	return amgDirichlet0
end

----------------------------------------------------------
-- function 'SetupFETISolver()' (first parameters only for (re)naming of logfile):
----------------------------------------------------------
function SetupFETISolver(str_problem,
			 dim,
			 linMaxIterations,
			 numProcs,
			 u,                         -- for testvector writer for FAMG (created by 'CreateAMGTestvector()')
			 dirichletBND, approxSpace, -- for testvector writer for FAMG (created by 'CreateAMGTestvectorDirichlet0()'); 'approxSpace' also for 'GridFunctionDebugWriter()'
			 activateDbgWriter,
			 verbosity, logfileName)

	print("    Setting up FETI solver (begin of 'SetupFETISolver()') ...")

	--------------------------------------------------------------------------------
	-- preconditioners, convergence checks, sub solvers and FETI solver:
	--------------------------------------------------------------------------------

	-- coarse problem solver stuff:
        local cpJac = Jacobi()
        local cpGS  = GaussSeidel()
        local cpILU = ILU()

        local cpConvCheck
        local coarseproblemSolver

	-- Neumann problem solver stuff:
        local npJac = Jacobi()
        local npGS  = GaussSeidel()
        local npILU = ILU()

        local neumannConvCheck
        local neumannSolver

	-- Dirichlet problem solver stuff:
        local dpJac = Jacobi()
        local dpGS  = GaussSeidel()
        local dpILU = ILU()

        local dirichletConvCheck
        local dirichletSolver

	-- Stuff for the FETI-DP solver itself:
        local fetiConvCheck
	local linAbsLimit  = 1e-16 -- before 25102011: 1e-7
	local linReduction = 1e-7  -- acc. to Klawonn&Rheinbach, 2010; before 25102011: 1e-16

        local fetidbgWriter

        local fetiSolver
        
     -- options for AMG
        local bWriteMat = util.HasParamOption("-AMGwriteMat")
        local bAggressiveCoarsening = false


	--------------------------------------------------------------------------------
	-- Checking for FETI solver related parameters (begin)
	-- (preconditioners are fix - edit script!)
	--------------------------------------------------------------------------------
	-- get number of processes per subdomain:
	local numProcsPerSubdomain       = util.GetParamNumber("-numPPSD", 1)

	-- types of sub solvers:
	local coarseProblemSolverType    = util.GetParam("-cps", "exact") -- choose one in ["exact" | "cg" | "hlib" ]
	local neumannProblemSolverType   = util.GetParam("-ns",    "cg") -- choose one in ["exact" | "ls" | "cg" | "bicg" | "rsamg" | "famg" ]
	local dirichletProblemSolverType = util.GetParam("-ds",    "cg") -- choose one in ["exact" | "ls" | "cg" | "bicg" | "rsamg" | "famg" ]


	-- Display parameters (or defaults):
	print("       FETI solver related parameters chosen (or defaults):")
	print("          numPPSD (numProcsPerSubdomain)   = " .. numProcsPerSubdomain)
	
	print("          cps (coarseProblemSolverType)    = " .. coarseProblemSolverType)
	print("          ns  (neumannProblemSolverType)   = " .. neumannProblemSolverType)
	print("          ds  (dirichletProblemSolverType) = " .. dirichletProblemSolverType)
	
	--------------------------------------------------------------------------------
	-- Checking for FETI solver related parameters (end)
	--------------------------------------------------------------------------------
	
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	-- Gather info for domain decomposition
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------

	-- check number of processes
	if numProcs < 2 then
		print("WARNING: number of processes is smaller than 2 - huh??")
	end
	
	-- check number of processes per subdomain
	write("       Check if numPPSD = '" .. numProcsPerSubdomain .. "' process(es) per subdomain makes sense ... " )
	
	if not util.IsPowerOfTwo(numProcsPerSubdomain) then
		write("WARNING: numPPSD = '" .. numProcsPerSubdomain .. "' is not a power of 2! " )
	--	exit()
	end
	
	-- compute number of subdomains
	numSubdomains = numProcs / numProcsPerSubdomain
	
	-- check that numSubdomains is greater 1 && \in \N && a power of 2.
	if numSubdomains < 2 then
		print("ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is smaller than 2!? Aborting!" )
		exit()
	end

	if not util.IsNaturalNumber(numSubdomains) then
		print("ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is NOT a natural number!? Aborting!" )
		exit()
	end
	if not util.IsPowerOfTwo(numSubdomains) then
		write("WARNING: numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is not a power of 2! " )
	-- TODO: Maybe switch to a default value then
	--	exit() -- in this case the partition can be quite erratic (at least on small (triangular) grids)..
	end
	print(" ok!")
	
	print("          numSubDomains is " .. numSubdomains .. " (numProcs = " .. numProcs .. ", numProcsPerSubdomain = " .. numProcsPerSubdomain .. ")")
	--------------------------------------------------------------------------------
	
	-- create subdomain info
	print("       Create domainDecompInfo")
	domainDecompInfo = StandardDomainDecompositionInfo()
	domainDecompInfo:set_num_subdomains(numSubdomains)

	-- test one to many interface creation
	if verbosity >= 1 then
		for i=0,numProcs-1 do
			print("       subdom of proc " .. i .. ": " .. domainDecompInfo:map_proc_id_to_subdomain_id(i))
		end
	end
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	-- Setup of FETI solver:
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	
	--if verbosity >= 1 then
		print("       Setup of FETI coarse and sub problem solvers ...")
	--end
	----------------------------------------------------------
	-- preconditioners used by FETI sub solvers
	----------------------------------------------------------
	cpJac = Jacobi()
	cpJac:set_damp(0.8)
	cpGS  = GaussSeidel()
	cpILU = ILU()
	cpSGS = SymmetricGaussSeidel()
	
	npJac = Jacobi()
	npJac:set_damp(0.8)
	npGS  = GaussSeidel()
	npILU = ILU()
	npSGS = SymmetricGaussSeidel()
	
	dpJac = Jacobi()
	dpJac:set_damp(0.8)
	dpGS  = GaussSeidel()
	dpILU = ILU()
	dpSGS = SymmetricGaussSeidel()
	
	----------------------------------------------------------
	-- create and configure coarse problem solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("          Create and configure coarse problem solver of type '" .. coarseProblemSolverType .. "' ...")
	--end
	-- choose solver for coarse problem
	if coarseProblemSolverType == "exact" then
	
		coarseproblemSolver = LU()
	
	elseif coarseProblemSolverType == "cg" then
	
		coarseproblemSolver = CG()
		coarseproblemSolver:set_preconditioner(cpILU) -- Absturz: "Signal: Bus error (10)", "Signal code:  (5583)", etc ...
		--coarseproblemSolver:set_preconditioner(cpJac) -- "ERROR in 'JacobiPreconditioner::apply': Cannot change parallel status of correction to consistent."
		--coarseproblemSolver:set_preconditioner(cpGS) -- "ERROR in 'Gauss-Seidel::apply': Cannot change parallel storage type of correction to consistent."
		-- Ganz ohne preconditioner: "Cannot convert z to consistent vector."
	
	elseif coarseProblemSolverType == "hlib" then
	
		coarseproblemSolver = HLIB() -- create HLIB Solver
		coarseproblemSolver:set_hlib_accuracy_H(1.e-4)
	
	else
		print ("ERROR: coarse problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the coarse problem solver
	cpConvCheck = ConvCheck()
	cpConvCheck:set_maximum_steps(20)
	cpConvCheck:set_minimum_defect(1e-10)
	cpConvCheck:set_reduction(1e-16)
	cpConvCheck:set_verbose(false)
	
	coarseproblemSolver:set_convergence_check(cpConvCheck)
	
	----------------------------------------------------------
	-- create and configure sub problem solvers
	----------------------------------------------------------

	----------------------------------------------------------
	-- create and configure Neumann problem solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("          Create and configure Neumann problem solver of type '" .. neumannProblemSolverType .. "' ...")
	--end
	-- choose solver for Neumann problems
	if neumannProblemSolverType == "exact" then
	
		neumannSolver = LU()
	
	elseif neumannProblemSolverType == "ls" then
	
		neumannSolver = LinearSolver()
		neumannSolver:set_preconditioner(npILU) -- npJac
	
	elseif neumannProblemSolverType == "cg" then
	
		neumannSolver = CG()
		neumannSolver:set_preconditioner(npILU) -- npJac
		--neumannSolver:set_preconditioner(npSGS)
	
	elseif neumannProblemSolverType == "bicg" then
	
		neumannSolver = BiCGStab() -- ERROR in 'PrimalSubassembledMatrixInverse::init': Could not solve local problem to compute Schur complement w.r.t. primal unknowns. - war das bei 1x1?
		-- BiCGStab statt CG hatte z.T. eine Verdopplung bis Verdreifachung der Iterationszahlen
		-- (apply_F), z.T. eine Verringerung von einem Viertel bis zur Haelfte (backsolve; compute_d) bewirkt!?
		neumannSolver:set_preconditioner(npILU) -- npJac
	
	elseif neumannProblemSolverType == "rsamg" or neumannProblemSolverType == "famg" then

		if neumannProblemSolverType == "famg" then

			npAMG = FAMGPreconditioner()	
			npAMG:set_delta(0.5)
			npAMG:set_theta(0.95)
			npAMG:set_aggressive_coarsening(bAggressiveCoarsening)
			npAMG:set_write_f_values(false)

			if bWriteMat then
				npAMG:write_testvectors(true)
			end
			local FAMGtestvectorSmoother = Jacobi()
			FAMGtestvectorSmoother:set_damp(0.66)
				
			npAMG:set_testvector_smooths(1)
			npAMG:set_damping_for_smoother_in_interpolation_calculation(0.66)
			npAMG:set_testvector_smoother(FAMGtestvectorSmoother)
--[[
			-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
			--testvectorwriter = CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
			testvectorwriter = CreateAMGTestvector(u, "ourTestvector2d_1_1", dim)
			testvector = GridFunction(approxSpace)
			testvectorwriter:update(testvector)	
			npAMG:add_testvector(testvectorwriter, 1.0)
   ]]				
			npAMG:set_testvector_from_matrix_rows()
				
			if bExternalCoarsening then
				npAMG:set_external_coarsening(true)
				npAMG:set_parallel_coarsening(GetColorCoarsening())
				-- npAMG:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
				-- npAMG:set_parallel_coarsening(GetRS3Coarsening())
				if dim == 2 then
					npAMG:set_strong_connection_external(0.2)
				else
					npAMG:set_strong_connection_external(0.1)
				end				
			end
			
			--if bWriteMat then
			--	npAMG:write_testvectors(true)
			--end
			npAMG:set_galerkin_truncation(1e-9)
			npAMG:set_H_reduce_interpolation_nodes_parameter(0.1)
			npAMG:set_prereduce_A_parameter(0.01)
		else
			print ("       Create RSMG as Neumann problem solver ... ")
			npAMG = RSAMGPreconditioner()
			-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
			-- amg:set_parallel_coarsening(GetColorCoarsening()) --
			npAMG:set_parallel_coarsening(GetRS3Coarsening()) --
			-- amg:set_parallel_coarsening(GetSimpleParallelCoarsening())
			if bAggressiveCoarsening then
				npAMG:enable_aggressive_coarsening_A(2)
			end
		end
		
		maxBase = util.GetParamNumber("-maxBase", 1000) -- TODO (maybe): make 'maxBase' different for "np solver" and "dp solver"!
		
		if bWriteMat then
			vectorWriter = GridFunctionPositionProvider()
			vectorWriter:set_reference_grid_function(u)
			npAMG:set_position_provider(vectorWriter)
			npAMG:set_matrix_write_path("/Users/mrupp/matrices/")
		end
		
		
		npAMGGS = GaussSeidel()
		npAMGBase = LU()
		
		npAMG:set_num_presmooth(3)
		npAMG:set_num_postsmooth(3)
		npAMG:set_cycle_type(1)
		npAMG:set_presmoother(npAMGGS)
		npAMG:set_postsmoother(npAMGGS)
		npAMG:set_base_solver(npAMGBase)
		npAMG:set_max_levels(20)
		npAMG:set_max_nodes_for_base(maxBase)
		npAMG:set_max_fill_before_base(0.7)
		npAMG:set_fsmoothing(true)
		
		print("\nThe AMG solver:")
		npAMG:tostring()
		print("")
		
		neumannSolver = LinearSolver()
		neumannSolver:set_preconditioner(npAMG)
	
	else
		print ("ERROR: Neumann problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the Neumann problem solver
	npsMaxIterations = 2000
	npsAbsLimit      = 1e-10
	npsReduction     = 1e-16
	print("             npsMaxIterations = " .. npsMaxIterations .. ", npsAbsLimit = " .. npsAbsLimit .. ", npsReduction = " .. npsReduction)
	neumannConvCheck = ConvCheck()
	neumannConvCheck:set_maximum_steps(npsMaxIterations)
	neumannConvCheck:set_minimum_defect(npsAbsLimit)
	neumannConvCheck:set_reduction(npsReduction)
	neumannConvCheck:set_verbose(true)
	
	neumannSolver:set_convergence_check(neumannConvCheck)
	
	----------------------------------------------------------
	-- create and configure Dirichlet problem solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("          Create and configure Dirichlet problem solver of type '" .. dirichletProblemSolverType .. "' ...")
	--end
	-- choose solver for Dirichlet problems
	if dirichletProblemSolverType == "exact" then
	
		dirichletSolver = LU()
	
	elseif dirichletProblemSolverType == "ls" then
	
		dirichletSolver = LinearSolver()
		dirichletSolver:set_preconditioner(dpILU) -- dpJac
	
	elseif dirichletProblemSolverType == "cg" then
	
		dirichletSolver = CG()
		dirichletSolver:set_preconditioner(dpILU) -- dpJac
		--dirichletSolver:set_preconditioner(npSGS)
	
	elseif dirichletProblemSolverType == "bicg" then
	
		dirichletSolver = BiCGStab() -- not yet tested for Dirichlet problem!
		dirichletSolver:set_preconditioner(dpILU) -- dpJac
	
	elseif dirichletProblemSolverType == "rsamg" or dirichletProblemSolverType == "famg" then

		if dirichletProblemSolverType == "famg" then

			dpAMG = FAMGPreconditioner()	
			dpAMG:set_delta(0.5)
			dpAMG:set_theta(0.95)
			dpAMG:set_aggressive_coarsening(bAggressiveCoarsening)
			dpAMG:set_write_f_values(false)
				
			if bWriteMat then
				dpAMG:write_testvectors(true)
			end
			local FAMGtestvectorSmoother = Jacobi()
			FAMGtestvectorSmoother:set_damp(0.66)
				
			dpAMG:set_testvector_smooths(1)
			dpAMG:set_damping_for_smoother_in_interpolation_calculation(0.66)
			dpAMG:set_testvector_smoother(FAMGtestvectorSmoother)
	--[[
			-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
			--testvectorwriter = CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
			testvectorwriter = CreateAMGTestvector(u, "ourTestvector2d_1_1", dim)
			testvector = GridFunction(approxSpace)
			testvectorwriter:update(testvector)	
			dpAMG:add_testvector(testvectorwriter, 1.0)
   ]]	
			dpAMG:set_testvector_from_matrix_rows()
				
			if bExternalCoarsening then
				dpAMG:set_external_coarsening(true)
				dpAMG:set_parallel_coarsening(GetColorCoarsening())
				-- dpAMG:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
				-- dpAMG:set_parallel_coarsening(GetRS3Coarsening())
				if dim == 2 then
					dpAMG:set_strong_connection_external(0.2)
				else
					dpAMG:set_strong_connection_external(0.1)
				end				
			end
			
			--if bWriteMat then
			--	dpAMG:write_testvectors(true)
			--end
			dpAMG:set_galerkin_truncation(1e-9)
			dpAMG:set_H_reduce_interpolation_nodes_parameter(0.1)
			dpAMG:set_prereduce_A_parameter(0.01)
		else
			print ("       Create RSMG as Dirichlet problem solver ... ")
			dpAMG = RSAMGPreconditioner()
			-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
			-- amg:set_parallel_coarsening(GetColorCoarsening()) --
			dpAMG:set_parallel_coarsening(GetRS3Coarsening()) --
			-- amg:set_parallel_coarsening(GetSimpleParallelCoarsening())
			if bAggressiveCoarsening then
				dpAMG:enable_aggressive_coarsening_A(2)
			end
		end
		
		maxBase = util.GetParamNumber("-maxBase", 1000) -- TODO (maybe): make 'maxBase' different for "np solver" and "dp solver"!
		
		if bWriteMat then
			vectorWriter = GridFunctionPositionProvider()
			vectorWriter:set_reference_grid_function(u)
			dpAMG:set_position_provider(vectorWriter)
			dpAMG:set_matrix_write_path("/Users/mrupp/matricesNew/")
		end
		
		dpAMGGS = GaussSeidel()
		dpAMGBase = LU()
		
		dpAMG:set_num_presmooth(3)
		dpAMG:set_num_postsmooth(3)
		dpAMG:set_cycle_type(1)
		dpAMG:set_presmoother(dpAMGGS)
		dpAMG:set_postsmoother(dpAMGGS)
		dpAMG:set_base_solver(dpAMGBase)
		dpAMG:set_max_levels(20)
		dpAMG:set_max_nodes_for_base(maxBase)
		dpAMG:set_max_fill_before_base(0.7)
		dpAMG:set_fsmoothing(true)
		
		print("\nThe AMG solver:")
		dpAMG:tostring()
		print("")
		
		dirichletSolver = LinearSolver()
		dirichletSolver:set_preconditioner(dpAMG)
	
	else
		print ("ERROR: Dirichlet problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the Dirichlet problem solver
	dpsMaxIterations = 2000
	dpsAbsLimit      = 1e-10
	dpsReduction     = 1e-16
	print("             dpsMaxIterations = " .. dpsMaxIterations .. ", dpsAbsLimit = " .. dpsAbsLimit .. ", dpsReduction = " .. dpsReduction)
--	print("#ANALYZER INFO")
	dirichletConvCheck = ConvCheck()
	dirichletConvCheck:set_maximum_steps(dpsMaxIterations)
	dirichletConvCheck:set_minimum_defect(dpsAbsLimit)
	dirichletConvCheck:set_reduction(dpsReduction)
if dirichletProblemSolverType == "famg" then -- TMP
	dirichletConvCheck:set_verbose(true) -- TMP
else -- TMP
	dirichletConvCheck:set_verbose(false)
end -- TMP
	dirichletSolver:set_convergence_check(dirichletConvCheck)
	
	----------------------------------------------------------
	-- create and configure FETI Solver
	----------------------------------------------------------
	--if verbosity >= 1 then
		print("          Create and configure FETI main solver ...")
	--end
	fetiSolver = FETI()

	fetiSolver:set_domain_decomp_info(domainDecompInfo)
	
	fetiSolver:set_neumann_solver(neumannSolver)
	fetiSolver:set_dirichlet_solver(dirichletSolver)
	fetiSolver:set_coarse_problem_solver(coarseproblemSolver)
	
	-- define convergence criteria for the FETI solver
	fetiConvCheck = ConvCheck()

	print("             linMaxIterations = " .. linMaxIterations .. ", linAbsLimit = " .. linAbsLimit .. ", linReduction = " .. linReduction)

	fetiConvCheck:set_maximum_steps(linMaxIterations)
	fetiConvCheck:set_minimum_defect(linAbsLimit)
	fetiConvCheck:set_reduction(linReduction)
	
	fetiSolver:set_convergence_check(fetiConvCheck)

	if activateDbgWriter >= 1 then
		print("       Setting debug writer for 'fetiSolver' (no make consistent (former: 'raw data')")
		-- debug writer
		fetidbgWriter = GridFunctionDebugWriter(approxSpace)
		fetidbgWriter:set_vtk_output(true)
		-- alt:fetidbgWriter:set_print_raw_data(true) -- if 'true' print "raw" data (no "make consistent" before printing): for checking temporary results in FETI
		fetidbgWriter:set_print_consistent(false) -- if 'false' print "raw" data (no "make consistent" before printing): for checking temporary results in FETI

		fetiSolver:set_debug(fetidbgWriter)		
		fetiSolver:set_test_one_to_many_layouts(true)
	end
	
	

--[[ Layout tests sind irgendwann rausgeflogen ...
	if verbosity >= 2 then
		print("       Performing layout tests:")
		--BuildDomainDecompositionLayoutsTest2d(u, domainDecompInfo);
		--OneToManyTests2d(u)
	end
]]

	----------------------------------------------------------
-- add feti specific name parts to new logfile name
	str_spsolvers = "spss" -- "sub problem solvers"
	str_nppsd     = "nppsd" -- "number of processes per Subdomain"

	str_nppsd     = str_nppsd     .. "-" .. numProcsPerSubdomain
	logfileName   = logfileName   .. "_" .. str_nppsd

	str_spsolvers = str_spsolvers .. "-" .. dirichletProblemSolverType
                                      .. "-" .. neumannProblemSolverType
                                      .. "-" .. coarseProblemSolverType
	logfileName   = logfileName   .. "_" .. str_spsolvers

	print("    Returning FETI solver 'fetiSolver', ready for application (end of 'SetupFETISolver()')!")
	return fetiSolver, fetiConvCheck, logfileName
end
