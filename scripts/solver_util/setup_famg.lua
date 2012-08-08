
-- create AMG ---
-----------------

function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
	local amgTestvector;
	amgTestvector = GridFunctionVectorWriter()
	amgTestvector:set_reference_grid_function(gridfunction)
	amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
	return amgTestvector	
end
	
function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	local amgDirichlet0 = GridFunctionVectorWriterDirichlet0()
	amgDirichlet0:init(dirichletBND, approxSpace)
	return amgDirichlet0
end

function SetupFAMGSolver(base, presmoother, postsmoother)
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