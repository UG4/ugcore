-- Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
-- Author: Andreas Vogel
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

--[[!
\file domain_disc_util.lua
\defgroup scripts_util_domaindisc Domain Discretization Utility
\ingroup scripts_util
\{
\brief functions to create DomainDiscs using a string disc-type identifier
]]--

--!	Returns a NeumannBoundary Element-Disc of the requested type
--! @return Returns the element discreatization (of the Neumann BC)
--! @param fcts (String) names of symbolic function 
--! @param discType (String) discretizatin scheme 
function NeumannBoundary(fcts, discType)
	if discType == nil then discType = "fv1" end
	if 		discType == "fv1"  then return NeumannBoundaryFV1(fcts)
	elseif  discType == "fv"   then return NeumannBoundaryFV(fcts)
	elseif  discType == "fe"   then return NeumannBoundaryFE(fcts)
	else 
		print("NeumannBoundary: no disc type '"..discType.."' available. Aborting")
		exit();
	end
end

--! Returns a IConvectionShape class object of an upwind method (for FV1 discretization)
--! @return Returns the object of the upwind method
--! @param upwindType type of the upwind as string
--! @param upwindParam object-specific parameters
function UpwindFV1(upwindType, upwindParam)
	if upwindType == nil then upwindType = "partial" end
	if		upwindType == "no" then return NoUpwind()
	elseif	upwindType == "full" then return FullUpwind ()
	elseif	upwindType == "weighted" then
		local upwind = WeightedUpwind ()
		if upwindParam ~= nil then
			upwind:set_weight(upwindParam)
		end
		return upwind
	elseif	upwindType == "partial" then return PartialUpwind ()
	else
		print("UpwindFV1: no upwind type '"..upwindType.."' available. Aborting")
		exit();
	end
end

--[[!
\}
]]--
