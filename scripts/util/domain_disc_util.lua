--[[!
\file scripts/util/domain_disc_util.lua
\brief	functions to create DomainDiscs using a string disc-type identifier
]]--


--!	Returns a ConvectionDiffusion Element-Disc of the requested type
--! @return Returns the domain discreatization
--! @param fcts (String) names of symbolic functions 
--! @param subsets (String) names of symbolic subsets 
--! @param discType (String) discretizatin scheme 
function ConvectionDiffusion(fcts, subsets, discType)
	if 		discType == "fv1"  then return ConvectionDiffusionFV1(fcts, subsets)
	elseif  discType == "fe"   then return ConvectionDiffusionFE(fcts, subsets)
	elseif  discType == "fvcr" then return ConvectionDiffusionFVCR(fcts, subsets)
	elseif  discType == "fv"   then return ConvectionDiffusionFV(fcts, subsets)
	else 
		print("ConvectionDiffusion: no disc type '"..discType.."' available. Aborting")
		exit();
	end
end

--!	Returns a ConstantEquation Element-Disc of the requested type
--! @return Returns the domain discreatization
--! @param fcts (String) names of symbolic functions 
--! @param subsets (String) names of symbolic subsets 
--! @param discType (String) discretizatin scheme 
function ConstantEquation(fcts, subsets, discType)
	if 		discType == "fv1"  then return ConstantEquationFV1(fcts, subsets)
	else 
		print("ConstantEquation: no disc type '"..discType.."' available. Aborting")
		exit();
	end
end