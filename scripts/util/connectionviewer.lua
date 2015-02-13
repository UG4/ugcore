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