-- Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
-- Author: Martin Rupp
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

----------------------------------------------------------
--[[!
--   \addtogroup scripts_util_solver
--   \{
--   \file setup_famg.lua
--   \brief Lua - Script which creates and configures the FAMG solver
--   \author Martin Rupp
--   \sa util.SetupFAMGPreconditioner
--
--   \note You have to enable the AMG plugin: cmake -Damg=ON ..
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
--! ug_load_script("solver_util/setup_famg.lua")
--! precond = util.SetupFAMGPreconditioner(LU(), GaussSeidel(), GaussSeidel())
--! \endcode
--! note that FAMG is a preconditioner, you can use it e.g. in a Linear Solver:
--! \code
--! precond = util.SetupFAMGPreconditioner()
--! linSolver = CG()
--! linSolver:set_preconditioner(precond)
--! linSolver:set_convergence_check(ConvCheck(40, 1e-16, 1e-9))
--! -- use linSolver
--! \endcode
--! \sa FAMG 
function util.SetupFAMGPreconditioner(base, presmoother, postsmoother)
	PluginRequired("amg")
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
	jac:set_damp(0.6)
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
	
	amg:set_galerkin_truncation(0.0)
	amg:set_H_reduce_interpolation_nodes_parameter(0.0)
	amg:set_prereduce_A_parameter(0.0)
	
	amg:set_write_f_values(false)
	
	return amg	
end

--[[!
\}
]]--
