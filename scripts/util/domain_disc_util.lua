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
--! @return Returns the domain discreatization
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

--!	Returns a NavierStokes Element-Disc of the requested type
--! @return Returns the domain discreatization
--! @param fcts (String) names of symbolic functions 
--! @param subsets (String) names of symbolic subsets 
--! @param discType (String) discretizatin scheme 
function NavierStokes(fcts, subsets, discType)
	if discType == nil then discType = "fv1" end
	if 		discType == "fv1"  then return NavierStokesFV1(fcts, subsets)
	elseif  discType == "fv"   then return NavierStokesFV(fcts, subsets)
	elseif  discType == "fvcr" then return NavierStokesFVCR(fcts, subsets)
	elseif  discType == "fe"   then return NavierStokesFE(fcts, subsets)
	elseif  discType == "fecr" then return NavierStokesFE(fcts, subsets)
	else 
		print("NavierStokes: no disc type '"..discType.."' available. Aborting")
		exit();
	end
end

--!	Returns a NavierStokesInflow Element-Disc of the requested type
--! @return Returns the domain discreatization
--! @param fcts (String) names of symbolic functions 
--! @param subsets (String) names of symbolic subsets 
--! @param discType (String) discretizatin scheme 
function NavierStokesInflow(spMaster)
	if spMaster == nil then
		print("NavierStokesInflow: master disc is not available. Aborting")
		exit();
	end

	local discType = spMaster:disc_type();	
	if 		discType == "fv1"  then return NavierStokesInflowFV1(spMaster)
	elseif  discType == "fv"   then return NavierStokesInflowFV(spMaster)
	elseif  discType == "fvcr" then return NavierStokesInflowFVCR(spMaster)
	elseif  discType == "fe"   then return NavierStokesInflowFE(spMaster)
	elseif  discType == "fecr" then return NavierStokesInflowFE(spMaster)
	else 
		print("NavierStokesInflow: no disc type '"..discType.."' available. Aborting")
		exit();
	end
end

--!	Returns a NavierStokesNoNormalStressOutflow Element-Disc of the requested type
--! @return Returns the domain discreatization
--! @param fcts (String) names of symbolic functions 
--! @param subsets (String) names of symbolic subsets 
--! @param discType (String) discretizatin scheme 
function NavierStokesNoNormalStressOutflow(spMaster)
	if spMaster == nil then
		print("NavierStokesNoNormalStressOutflow: master disc is not available. Aborting")
		exit();
	end

	local discType = spMaster:disc_type();	
	if 		discType == "fv1"  then return NavierStokesNoNormalStressOutflowFV1(spMaster)
	elseif  discType == "fvcr" then return NavierStokesNoNormalStressOutflowFVCR(spMaster)
	elseif  discType == "fv"   then return NavierStokesNoNormalStressOutflowFV(spMaster)
	else 
		print("NavierStokesNoNormalStressOutflow: no disc type '"..discType.."' available. Aborting")
		exit();
	end
end

--[[!
\}
]]--
