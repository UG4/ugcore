----------------------------------------------------------
--
--   Lua - Script which creates and configures a FETI solver (linear solver).
--
--   Author: Ingo Heppner (using FETI solver specific code from 'fetitest.lua'.
--   Begin: 06072011.
--
--   Provides: Function which set up a FETI solver object.
--          Input parameters:
--          'linMaxIterations', 'activateDbgWriter', 'verbosity'.
--          Return value:
--          Reference of the fully configured FETI solver object.
--
--   Usage example:
--	1. In the calling LUA script:
--
--	ug_load_script("setup_fetisolver.lua")
--	solver = SetupFETISolver(lsMaxIter,
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
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 0 -lsType feti -nPPSD 1"
openmpirun -np 4 ugshell $UGARGS -numRefs 5 # etc.
	
# cekon:
########
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -nPPSD 1"

salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 5 -logtofile feti-sd1_8x8-quad_prerefs3-refs5_pe04.txt - geht!
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 6 -logtofile feti-sd1_8x8-quad_prerefs3-refs6_pe04.txt - geht, zumindest mit 2000 Schritten Dirichlet-Solver!
salloc -n  4 mpirun ./ugshell $UGARGS -numRefs 7 -logtofile feti-sd1_8x8-quad_prerefs3-refs7_pe04.txt - "Could not solve Dirichlet problem (step 3.b)", auf allen 4 Procs!

salloc -n 16 mpirun ./ugshell $UGARGS -numRefs 6 -logtofile feti-sd1_8x8-quad_prerefs3-refs6_pe16.txt - geht!
salloc -n 16 mpirun ./ugshell $UGARGS -numRefs 7 -logtofile feti-sd1_8x8-quad_prerefs3-refs7_pe16.txt - "Could not solve Dirichlet problem (step 3.b)", nur auf Proc 8!?

salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 7 -logtofile feti-sd1_8x8-quad_prerefs3-refs7_pe64.txt - geht!
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 8 -logtofile feti-sd1_8x8-quad_prerefs3-refs8_pe64.txt - "Could not solve Dirichlet problem (step 3.b)", auf Proc 53, 55, 61, 4, 5, 1 und 16
(Das sind in der Tat alle, bei denen es schief geht - keine weiteren durch Auskommentieren des "return(false)" in der entsprechenden Bedingung!)

#Mit '-distType grid2d':
salloc -n 64 mpirun ./ugshell $UGARGS -numRefs 8 -distType grid2d -logtofile feti-sd1_8x8-quad_prerefs3-refs8_grid2d_pe64.txt - jetzt geht "step 3.b" auf folgenden Procs nicht:
39, 47, 31, 55, 60, 61, 62, 3, 1, 2

# JuGene:
#########
UGARGS="-ex ../scripts/tests/scalability_test.lua -dim 2 -grid ../data/grids/unit_square_01/unit_square_01_quads_8x8.ugx -lsMaxIter 100 -numPreRefs 3 -lsType feti -nPPSD 1"
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

]]


----------------------------------------------------------
-- function 'SetupFETISolver()':
----------------------------------------------------------
function SetupFETISolver(domain,
			 linMaxIterations,
			 numProcs,
			 activateDbgWriter,
			 verbosity)

	print("    'setup_fetisolver.lua': Setting up FETI solver...")

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
	local linAbsLimit = 1e-7
	local linReduction = 1e-16

        local fetiSolver


	--------------------------------------------------------------------------------
	-- Checking for FETI solver related parameters (begin)
	-- (preconditioners are fix - edit script!)
	--------------------------------------------------------------------------------
	-- get number of processes per subdomain:
	local numProcsPerSubdomain       = util.GetParamNumber("-nPPSD", 1)

	-- types of sub solvers:
	local coarseProblemSolverType    = util.GetParam("-cps", "exact") -- choose one in ["exact" | "cg" | "hlib" ]
	local neumannProblemSolverType   = util.GetParam("-nps",    "cg") -- choose one in ["exact" | "ls" | "cg" | "bicg" ]
	local dirichletProblemSolverType = util.GetParam("-dps",    "cg") -- choose one in ["exact" | "ls" | "cg" | "bicg" ]


	-- Display parameters (or defaults):
	print("    FETI solver related parameters chosen (or defaults):")
	print("        nPPSD (numProcsPerSubdomain)     = " .. numProcsPerSubdomain)
	
	print("        cps (coarseProblemSolverType)    = " .. coarseProblemSolverType)
	print("        nps (neumannProblemSolverType)   = " .. neumannProblemSolverType)
	print("        dps (dirichletProblemSolverType) = " .. dirichletProblemSolverType)
	
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
	
	if not util.IsPowerOfTwo(numProcsPerSubdomain) then
		print("WARNING: nPPSD = '" .. numProcsPerSubdomain .. "' is not a power of 2!" )
	--	return
	end
	
	print("    Check if nPPSD = '" .. numProcsPerSubdomain .. "' process(es) per subdomain makes sense ..." )
	
	-- compute number of subdomains
	numSubdomains = numProcs / numProcsPerSubdomain
	
	-- check that numSubdomains is greater 1 && \in \N && a power of 2.
	if numSubdomains < 2 then
		print("ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is smaller than 2!? Aborting!" )
		return
	end

	if not util.IsNaturalNumber(numSubdomains) then
		print("ERROR:   numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is NOT a natural number!? Aborting!" )
		return
	end
	
	if not util.IsPowerOfTwo(numSubdomains) then
		print("WARNING: numSubdomains = numProcs / numProcsPerSubdomain = '" .. numSubdomains .. "' is not a power of 2! Continuing ..." )
	-- TODO: Maybe switch to a default value then
	--	return -- in this case the partition can be quite erratic (at least on small (triangular) grids)..
	end
	
	print("    NumProcs is " .. numProcs .. ", NumSubDomains is " .. numSubdomains )
	--------------------------------------------------------------------------------
	
	-- create subdomain info
	print("    Create domainDecompInfo")
	domainDecompInfo = StandardDomainDecompositionInfo()
	domainDecompInfo:set_num_subdomains(numSubdomains)

	-- test one to many interface creation
	if verbosity >= 1 then
		for i=0,numProcs-1 do
			print("    subdom of proc " .. i .. ": " .. domainDecompInfo:map_proc_id_to_subdomain_id(i))
		end
	end
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	-- Setup of FETI solver:
	--------------------------------------------------------------------------------
	--------------------------------------------------------------------------------
	
	----------------------------------------------------------
	-- preconditioners used by FETI sub solvers
	----------------------------------------------------------
	cpJac = Jacobi()
	cpJac:set_damp(0.8)
	cpGS  = GaussSeidel()
	cpILU = ILU()
	
	npJac = Jacobi()
	npJac:set_damp(0.8)
	npGS  = GaussSeidel()
	npILU = ILU()
	
	dpJac = Jacobi()
	dpJac:set_damp(0.8)
	dpGS  = GaussSeidel()
	dpILU = ILU()
	
	----------------------------------------------------------
	-- create and configure coarse problem solver
	----------------------------------------------------------
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
	cpConvCheck = StandardConvergenceCheck()
	cpConvCheck:set_maximum_steps(2000)
	cpConvCheck:set_minimum_defect(1e-10)
	cpConvCheck:set_reduction(1e-16)
	cpConvCheck:set_verbose_level(false)
	
	coarseproblemSolver:set_convergence_check(cpConvCheck)
	
	----------------------------------------------------------
	-- create and configure Neumann solver
	----------------------------------------------------------
	-- choose solver for Neumann problems
	if neumannProblemSolverType == "exact" then
	
		neumannSolver = LU()
	
	elseif neumannProblemSolverType == "ls" then
	
		neumannSolver = LinearSolver()
		--neumannSolver:set_preconditioner(npJac)
		neumannSolver:set_preconditioner(npILU)
	
	elseif neumannProblemSolverType == "cg" then
	
		neumannSolver = CG()
		--neumannSolver:set_preconditioner(npJac)
		neumannSolver:set_preconditioner(npILU)
	
	elseif neumannProblemSolverType == "bicg" then
	
		neumannSolver = BiCGStab() -- ERROR in 'PrimalSubassembledMatrixInverse::init': Could not solve local problem to compute Schur complement w.r.t. primal unknowns. - war das bei 1x1?
		-- BiCGStab statt CG hatte z.T. eine Verdopplung bis Verdreifachung der Iterationszahlen
		-- (apply_F), z.T. eine Verringerung von einem Viertel bis zur Haelfte (backsolve; compute_d) bewirkt!?
	
	else
		print ("ERROR: Neumann problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the coarse problem solver
	neumannConvCheck = StandardConvergenceCheck()
	neumannConvCheck:set_maximum_steps(2000)
	neumannConvCheck:set_minimum_defect(1e-10)
	neumannConvCheck:set_reduction(1e-16)
	neumannConvCheck:set_verbose_level(false)
	
	neumannSolver:set_convergence_check(neumannConvCheck)
	
	----------------------------------------------------------
	-- create and configure Dirichlet solver
	----------------------------------------------------------
	-- choose solver for Dirichlet problems
	if dirichletProblemSolverType == "exact" then
	
		dirichletSolver = LU()
	
	elseif dirichletProblemSolverType == "ls" then
	
		dirichletSolver = LinearSolver()
		dirichletSolver:set_preconditioner(dpILU)
	
	elseif dirichletProblemSolverType == "cg" then
	
		dirichletSolver = CG()
		dirichletSolver:set_preconditioner(dpILU)
	
	elseif dirichletProblemSolverType == "bicg" then
	
		dirichletSolver = BiCGStab() -- not yet tested for Dirichlet problem!
		dirichletSolver:set_preconditioner(dpILU)
	
	elseif dirichletProblemSolverType == "rsamg" then

		maxBase = util.GetParamNumber("-maxBase", 1000)
			
		dpRSAMG = RSAMGPreconditioner()
		
		dpRSAMGGS = GaussSeidel()
		dpRSAMGBase = LU()
		
		dpRSAMG:set_num_presmooth(3)
		dpRSAMG:set_num_postsmooth(3)
		dpRSAMG:set_cycle_type(1)
		dpRSAMG:set_presmoother(dpRSAMGGS)
		dpRSAMG:set_postsmoother(dpRSAMGGS)
		dpRSAMG:set_base_solver(dpRSAMGBase)
		dpRSAMG:set_max_levels(2)
		dpRSAMG:set_max_nodes_for_base(maxBase)
		dpRSAMG:set_max_fill_before_base(0.7)
		dpRSAMG:set_fsmoothing(true)
		dpRSAMG:set_epsilon_truncation(0)		
		dpRSAMG:tostring()	
		
		dirichletSolver = LinearSolver()
		dirichletSolver:set_preconditioner(dpRSAMG)
		
	
	else
		print ("ERROR: Dirichlet problem solver not specified ==> exit")
		exit()
	end
	
	-- define convergence criteria for the coarse problem solver
	dirichletConvCheck = StandardConvergenceCheck()
	dirichletConvCheck:set_maximum_steps(100)
	dirichletConvCheck:set_minimum_defect(1e-10)
	dirichletConvCheck:set_reduction(1e-16)
	dirichletConvCheck:set_verbose_level(false)
	
	dirichletSolver:set_convergence_check(dirichletConvCheck)
	
	----------------------------------------------------------
	-- create and configure FETI Solver
	----------------------------------------------------------
	fetiSolver = FETI()

	fetiSolver:set_domain_decomp_info(domainDecompInfo)
	
	fetiSolver:set_neumann_solver(neumannSolver)
	fetiSolver:set_dirichlet_solver(dirichletSolver)
	fetiSolver:set_coarse_problem_solver(coarseproblemSolver)
	
	-- define convergence criteria for the coarse problem solver
	fetiConvCheck = StandardConvergenceCheck()

	print("    'setup_fetisolver.lua': linMaxIterations = " .. linMaxIterations .. ", linAbsLimit = " .. linAbsLimit .. ", linReduction = " .. linReduction)

	fetiConvCheck:set_maximum_steps(linMaxIterations)
	fetiConvCheck:set_minimum_defect(linAbsLimit)
	fetiConvCheck:set_reduction(linReduction)
	
	fetiSolver:set_convergence_check(fetiConvCheck)

	if activateDbgWriter >= 1 then
		fetiSolver:set_debug(dbgWriter)
	end
	
--[[ Layout tests sind irgendwann rausgeflogen ...
	if verbosity >= 2 then
		print("    Performing layout tests:")
		--BuildDomainDecompositionLayoutsTest2d(u, domainDecompInfo);
		--OneToManyTests2d(u)
	end
]]

	----------------------------------------------------------
	print("    'setup_fetisolver.lua': returning FETI solver 'fetiSolver', ready for application!")
	return fetiSolver
end
