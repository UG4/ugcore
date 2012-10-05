----------------------------------------------------------
--[[!
--   \file setup_famg.lua
--   \brief Lua - Script which creates and configures the FAMG solver
-- 	 \author Martin Rupp
--  \sa util.SetupFAMGPreconditioner
--
--  \note You have to enable the AMG plugin: cmake -Damg=ON ..
!]]--


--! helper function for SetupFAMGPreconditioner
function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
	local amgTestvector;
	amgTestvector = GridFunctionVectorWriter()
	amgTestvector:set_reference_grid_function(gridfunction)
	amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
	return amgTestvector	
end

--! helper function for SetupFAMGPreconditioner	
function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	local amgDirichlet0 = GridFunctionVectorWriterDirichlet0()
	amgDirichlet0:init(dirichletBND, approxSpace)
	return amgDirichlet0
end

util = util or {}

--! creates a FAMG Preconditioner object
--! \param base ILinearOperatorInverse to solve coarse level problem. if nil or omitted: LU 
--! \param presmoother ILinearIterator for presmoothing. if nil or omitted: SymmetricGaussSeidel 
--! \param postsmoother ILinearIterator for postsmoothing. if nil or omitted: SymmetricGaussSeidel
--!
--! example usage:
--! \code
--! ug_load_script("setup_famg.lua")
--! precond = util.SetupFAMGPrecondition(LU(), GaussSeidel(), GaussSeidel())
--! \endcode
--! note that FAMG is a preconditioner, you can use it e.g. in a Linear Solver:
--! \code
--! precond = util.SetupFAMGPrecondition()
--! linSolver = CG()
--! linSolver:set_preconditioner(precond)
--! linSolver:set_convergence_check(ConvCheck(40, 1e-16, 1e-9))
--! -- use linSolver
--! \endcode
--! \sa FAMG 
function util.SetupFAMGPreconditioner(base, presmoother, postsmoother)
	if base == nil then
		base = LU()
	end
	if presmoother == nil then
		presmoother = SymmetricGaussSeidel()
	end
	if postsmoother == nil then
		postsmoother = presmoother
	end
	amg = FAMGPreconditioner()
	
	amg:set_delta(0.1)
	amg:set_theta(0.9)

	local jac = Jacobi()
	jac:set_damp(0.8)
	amg:set_testvector_smooths(3)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
	amg:set_testvector_smoother(jac)
	amg:set_testvector_from_matrix_rows()
	
	amg:set_presmoother(presmoother)	
	amg:set_postsmoother(postsmoother)
	
	amg:set_num_presmooth(5)
	amg:set_num_postsmooth(5)	
	
	amg:set_base_solver(base)

	amg:set_min_nodes_on_one_processor(100)
	amg:set_preferred_nodes_on_one_processor(2000)

	amg:set_max_nodes_for_base(400)
	amg:set_max_fill_before_base(0.4)
	amg:set_fsmoothing(true)
	
	amg:set_galerkin_truncation(1e-9)
	amg:set_H_reduce_interpolation_nodes_parameter(0.0)
	amg:set_prereduce_A_parameter(0.0)
	
	amg:set_write_f_values(false)
	
	return amg	
end


