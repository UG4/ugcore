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
--   \file setup_rsamg.lua
--   \brief Lua - Script which creates and configures the RSAMG solver
-- 	 \author Martin Rupp
--  \sa util.SetupRSAMGPreconditioner
--
--  \note You have to enable the AMG plugin: cmake -Damg=ON ..
!]]--

util = util or {}

--! creates a RSAMG Preconditioner object
--! \param base ILinearOperatorInverse to solve coarse level problem. if nil or omitted: LU 
--! \param presmoother ILinearIterator for presmoothing. if nil or omitted: SymmetricGaussSeidel 
--! \param postsmoother ILinearIterator for postsmoothing. if nil or omitted: SymmetricGaussSeidel
--!
--! example usage:
--! \code
--! ug_load_script("solver_util/setup_rsamg.lua")
--! precond = util.SetupRSAMGPreconditioner(LU(), GaussSeidel(), GaussSeidel())
--! \endcode
--! note that FAMG is a preconditioner, you can use it e.g. in a Linear Solver:
--! \code
--! precond = util.SetupRSAMGPreconditioner()
--! linSolver = CG()
--! linSolver:set_preconditioner(precond)
--! linSolver:set_convergence_check(ConvCheck(40, 1e-16, 1e-9))
--! -- use linSolver
--! \endcode
--! \sa RSAMG
function util.SetupRSAMGPreconditioner(base, presmoother, postsmoother)
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
	amg = RSAMGPreconditioner()
	if IsDefinedUG_PARALLEL() then
		amg:set_parallel_coarsening(GetColorCoarsening())
	end
	
	amg:set_presmoother(presmoother)	
	amg:set_postsmoother(postsmoother)
	
	amg:set_num_presmooth(5)
	amg:set_num_postsmooth(5)	
	amg:set_base_solver(base)
	--amg:enable_aggressive_coarsening_A(2)
	
	amg:set_epsilon_strong(0.1)
	amg:set_max_nodes_for_base(100)
	amg:set_max_fill_before_base(0.4)
	amg:set_prolongation_truncation(1e-4)
	amg:set_min_nodes_on_one_processor(1000)
	amg:set_preferred_nodes_on_one_processor(1000)
	
	amg:set_fsmoothing(true)
	
	return amg	
end

--[[!
\}
]]--
