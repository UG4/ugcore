-- Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

util = util or {}

function util_std_zero_2(x, y, t)
	return 0
end

function util_std_zero_3(x, y, z, t)
	return 0
end

function util_std_dirichlet_zero_2(x, y, t)
	return true, 0
end

function util_std_dirichlet_zero_3(x, y, z, t)
	return true, 0
end

function util_std_dirichlet_one_2(x, y, t)
	return true, 1
end

function util_std_dirichlet_one_3(x, y, z, t)
	return true, 1
end

--! this function ensures that elements in T which are
--! not set (i.e. T[i]=nil) are set to T[i] = defaultT[i]
--! this way we can use T as a parameter Table
--! and all parameters not given are set to default values
--! e.g. SetDefaultParameters(T, {smoother="jac", maxLevels=2}) 
function SetDefaultParameters(T, defaultT)
	for i, v in pairs(defaultT) do
		if T[i] == nil then T[i] = v end
	end
end

--! checks in the table p that all elements of parameters
--! are non-nil. if some are nil, an error is thrown
function NeededParameters(p, parameters)
	for i, v in pairs(parameters) do
		if p[v] == nil then
			print("\ntable ")
			print(p)
			print("needs parameter \""..v.."\" !")
			ug_assert(false, "table needs parameter \""..v.."\" !")
		end	
	end
end

function util.EasyLoadGrid(p)
	NeededParameters(p, {"gridName", "requiredSubsets"})
	
	-- refinements:
	numPreRefs = util.GetParamNumber("-numPreRefs", p.numPreRefs or 0, "number of refinements before parallel distribution")
	numRefs    = util.GetParamNumber("-numRefs",    p.numRefs or 0, "number of refinements")
			
	if numPreRefs > numRefs then
		print("WARNING: numPreRefs = "..numPreRefs.." > numRefs = "..numRefs..". Setting now numPreRefs = numRefs = "..numRefs) 
	--	numPreRefs = numRefs 
	end
		
	dom = util.CreateAndDistributeDomain(p.gridName, numRefs, numPreRefs, p.requiredSubsets)
	
	print(dom:domain_info():to_string())
		
	local stretchX = util.GetParamNumber("-stretchX", 1.0)
	local stretchY = util.GetParamNumber("-stretchY", 1.0)
	local stretchZ = util.GetParamNumber("-stretchZ", 1.0)
	if stretchX ~= 1.0 or stretchY ~= 1.0 or stretchZ ~= 1.0 then
		print("scaling domain by ("..stretchX..", "..stretchY..")") 
		ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, stretchZ))
	end
	
	local distortDomain = util.GetParamNumber("-distortDomain", 0.0)
	if distortDomain > 0.0 then
		distortDomain = math.pow(distortDomain, numRefs+1)
		print("distortDomain ="..distortDomain)
		RandomizeDomain(dom, distortDomain, distortDomain, distortDomain)
	end			
		
	return dom
end

function util.EasyCreateApproxSpace(dom, fct)
	-- create Approximation Space
	local approxSpace = ApproximationSpace(dom)
	for i, v in pairs(fct) do
		approxSpace:add_fct(v, "Lagrange", 1)
	end    
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	--approxSpace:print_local_dof_statistic(2)
	--approxSpace:print_layout_statistic()
	--approxSpace:print_statistic()
	
	return approxSpace
end

function util.EasyDirichlet_0_Boundary(p)
	NeededParameters(p, {"dim", "boundaries"})
	
	if p.fct == nil then p.fct = {"c"} end
	local dirichletBND = DirichletBoundary()
	for i,v in pairs(p.boundaries) do
		for ifct,vfct in pairs(p.fct) do
			dirichletBND:add("util_std_dirichlet_zero_"..p.dim, vfct, v)
		end
	end
	return dirichletBND
end