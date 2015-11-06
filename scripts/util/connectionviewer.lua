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
	
util.ConnectionViewerPath = ug_get_root_path().."/externals/ConnectionViewer/ConnectionViewer.app/Contents/Resources/java/"

-- example:
-- 	util.ConnectionViewerExport("Stiffness.mat", {height=800, width=1000, scaleZoom=0.99, 
--		drawDiffusion=true, drawConvection=true, arrowConnections=true, exportPDF="Stiffness.pdf"})

function util.ConnectionViewerExport(matFile, p)
	if ProcRank() ~= 0 then return end
	
	function addIf(str, p, name)
		if p[name] == nil then
			return str
		else
			return str .." -"..name.." "..p[name]
		end
	end
	
	function addIfB(str, p, name)
		if p[name] == nil then
			return str
		else
			if p[name] then
				return str .." -"..name.." 1"
			else
				return str .." -"..name.." 0"
			end
		end
	end

	local CVoptions = matFile.." "
	CVoptions = addIf(CVoptions, p, "scaleZoom")
	CVoptions = addIf(CVoptions, p, "height")
	CVoptions = addIf(CVoptions, p, "width")
	
	CVoptions = addIf(CVoptions, p, "arrowSize")
	CVoptions = addIf(CVoptions, p, "fontsize")
	CVoptions = addIf(CVoptions, p, "zcompression")	
	
	CVoptions = addIfB(CVoptions, p, "drawConvection")
	CVoptions = addIfB(CVoptions, p, "drawDiffusion")
	CVoptions = addIfB(CVoptions, p, "arrowConnections")
	CVoptions = addIfB(CVoptions, p, "automaticReload")
	CVoptions = addIfB(CVoptions, p, "showParallelNodes")
	
	CVoptions = addIfB(CVoptions, p, "drawConnections")
	if p.exportPDF then
		CVoptions = CVoptions.." -exportPDF "..p.exportPDF
	end
	if p.exportTex then
		CVoptions = CVoptions.." -exportTex "..p.exportTex
	end
	CVoptions = CVoptions.." -quit"
	local CVcall = "java -cp "..util.ConnectionViewerPath.."/ConnectionViewer.jar connectionviewer.ConnectionViewer "..CVoptions
	--print(CVcall)
	os.execute(CVcall)
	
	
end