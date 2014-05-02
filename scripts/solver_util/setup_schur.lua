----------------------------------------------------------
--[[!
--   \addtogroup scripts_util_solver
--   \{
--   \file setup_schur.lua
--   \brief Lua - Script which creates and configures the Schur Preconditioner
-- 	 \author Martin Rupp
--
!]]--

util = util or {}
util.schur = util.schur or {}

-- local schurType = util.GetParam("-schurType", "Full")
-- local skeletonSolverType = util.GetParam("-schurSkeletonSolver", "ILU")	


--! creates a Schur Preconditioner object
--! parameters are handle via command line:
--! @param schurType 			how to handle the schur complement
--! @param schurSkeletonSolver  how to solve the schur complement
--! schurType can be [Full | AGG | BiCGStab ]
--!   Full: Full calculated Schur complement.
--!   AGG: A_GammaGamma as approximation of the schur complement
--!   BiCGStab: Solve the schur complement with BiCGStab. Schur Complement
--!   is only needed as a operator (non-matrix-based solver). 
--!     (no schurSkeletonSolver needed)
--! schurSkeletonSolver can be [ SparseLU | LapackLU | LU | Jacobi | SGS]
--!   LapackLU, SparseLU: direct Solvers. Only for schurType FULL 
--!   Jacobi, SGS, ILU: BiCGStab+Precond.
--! example usage:
--! \code
--! ug_load_script("solver_util/setup_schur.lua")
--! precond = util.schur.GetPreconditioner()
--! \endcode
--! note that FAMG is a preconditioner, you can use it e.g. in a Linear Solver:
--! \code
--! precond = util.schur.GetPreconditioner()
--! linSolver = BiCGStab()
--! linSolver:set_preconditioner(precond)
--! linSolver:set_convergence_check(ConvCheck(40, 1e-16, 1e-9))
--! -- use linSolver
--! \endcode
function util.schur.GetPreconditioner(schurType, skeletonSolverType, schur_dirichlet_solver)
	if schurType == nil then schurType = "Full" end
	if skeletonSolverType == nil then skeletonSolverType = "ILU" end
	if schur_dirichlet_solver == nil then
		schur_dirichlet_solver = LU()
	end	
	
	if schurType == "Full" or schurType == "AGG" then
		if 	skeletonSolverType == "SparseLU" then
			ug_assert(schurType == "Full", "not possible with AGG")
			skeletonSolver = AgglomeratingSolver(LU())
		elseif skeletonSolverType == "LapackLU" then
			ug_assert(schurType == "Full", "not possible with AGG")
			schurLU = LU()
			schurLU:set_minimum_for_sparse(-1) -- no sparse
			schurLU:set_info(true)
			skeletonSolver = AgglomeratingSolver(schurLU)
		elseif 	skeletonSolverType == "SuperLU" then
			skeletonSolver = AgglomeratingSolver(SuperLU())
		elseif skeletonSolverType == "Jacobi" then
			skeletonSolver = BiCGStab()
			skeletonSolver:set_preconditioner(Jacobi())
		elseif skeletonSolverType == "SGS" then
			skeletonSolver = BiCGStab()
			skeletonSolver:set_preconditioner(SymmetricGaussSeidel())
		elseif skeletonSolverType == "Schur" then
			skeletonSolver = BiCGStab()
			schur2 = SchurComplement()
			schur2:set_dirichlet_solver(LU())
			schur2:set_skeleton_solver(SchurInverseWithFullMatrix(AgglomeratingSolver(LU()) ) )
			skeletonSolver:set_preconditioner(schur2)	
		elseif skeletonSolverType == "ILU" then
			skeletonSolver = BiCGStab()
			skeletonSolver:set_preconditioner(ILU())
			skeletonSolver:set_convergence_check(ConvCheck(1000, 1e-12, 1e-10, false))	
		else
			error("unknown skeletonSolver "..skeletonSolverType)
		end
	
		if schurType == "AGG" then
			schur_skeleton_solver = SchurInverseWithAGammaGamma(skeletonSolver)
		else
			schur_skeleton_solver = SchurInverseWithFullMatrix(skeletonSolver)
		end	
	elseif schurType == "auto" then	
		skeletonSolver = AutoLinearSolver(1e-8, 1e-8)
		skeletonSolver:set_convergence_check(ConvCheck(200, 1e-12, 1e-12, true))
		if 	skeletonSolverType == "SparseLU" then
			skeletonSolver:set_preconditioner(AgglomeratingPreconditioner(ILUT(0.0)))		
		elseif skeletonSolverType == "ILU" then			
			skeletonSolver:set_preconditioner(ILU())	
		end			
		schur_skeleton_solver = SchurInverseWithAutoFullMatrix(skeletonSolver)
	elseif schurType == "BiCGStab" then
		schur_skeleton_solver = SchurInverseWithOperator(BiCGStab())
	else
		error("unknown schurType "..schurType..". known types are [Full | AGG | BiCGStab ]")
	end


	local schur = SchurComplement()
	schur:set_dirichlet_solver(schur_dirichlet_solver)
	schur:set_skeleton_solver(schur_skeleton_solver)
	print("The Schur Preconditioner:")
	print(schur:config_string())
	return schur
end



--! creates a Schur Solver object
--! see also util.schur.GetPreconditioner
--! example usage:
--! \code
--! ug_load_script("solver_util/setup_schur.lua")
--! solver = util.schur.GetSolver()
--! \endcode
--! \sa util.schur.GetPreconditioner
function util.schur.GetSolver(schurType, skeletonSolverType, schur_dirichlet_solver)
	local linSolver = LinearSolver()
	linSolver.set_preconditioner(util.schur.GetPreconditioner(schurType, skeletonSolverType, schur_dirichlet_solver))
	return linSolver
end

--[[!
\}
]]--

