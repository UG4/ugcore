--[[!
\file gnuplot.lua
\defgroup scripts_util_gnuplot GNUPlot Utility
create a globally seen package (all non-local functions call then be called using gnuplot.xxx.)
\ingroup scripts_util
\{
]]--

gnuplot = gnuplot or {}

ug_load_script("util/gnuplot_rgb_color.lua")

gnuplot.terminals = nil
--! list of all available terminals
function gnuplot.available_terminals()

	if not gnuplot.terminals then
		gnuplot.terminals = {}
		local termFile = os.tmpname()
		os.execute("echo | gnuplot -e 'set terminal' &> "..termFile)
		for line in io.lines(termFile) do
			line = string.gsub (line, "Press return for more:", "")
		 	table.insert(gnuplot.terminals, string.match (line, "%s+(%w*)%s+"))
		end
		os.remove(termFile)
	end
	
	return gnuplot.terminals
end

--[[! writes the passed data array to file

filename			the output file
data				table containing the data
passRows		    flag indication if data array stores rows (default = false)
--]]
function gnuplot.write_data(filename, data, passRows)

	-- default is table data by row
	passRows = passRows or false
	forNil = "--"

	-- check data
	if type(data) ~= "table" then
		io.stderr:write("Gnuplot Error: Expect table as data.\n");
		return 1	
	end
	if data[1]==nil then return 0; end

	-- pure data is rowwise
	if type(data[1]) == "number" then passRows = true end

	-- open file
	local file = io.open(filename, "w+")
	if (not file) then
		io.stderr:write("Gnuplot Error: cannot open output file: '")
		io.stderr:write(filename .. " '\n");
		return 1
	end
	
	-- case 1: data is given column per column
	if not passRows then
	
		-- check columns
		local bHasMatrixData = false
		local minRow = math.huge
		local maxRow = -math.huge
		for _, column in ipairs(data) do

			-- check that tables passed as coulumns
			if not(type(column) == "table") then
				io.stderr:write("Gnuplot Error: Data array must contain only tables.\n");
				return 1						
			end
	
			-- check if column data is actually a matrix	
			-- or extract row range		
			if type(column[1]) == "table" then
				bHasMatrixData = true		
			else
				minRow = math.min(minRow, table.imin(column))
				maxRow = math.max(maxRow, table.imax(column))		
			end
		end

		-- write column
		if not(bHasMatrixData) then
			for row = minRow, maxRow do
				for col, _ in ipairs(data) do
					local item = data[col][row]
					if item ~= nil then file:write(item, " ")
					else 				file:write(forNil, " ") end
				end
				file:write("\n")
			end		
		else
			if type(data[1][1]) == "table" or type(data[2][1]) == "table" then
				io.stderr:write("Gnuplot: Format for matrix data: x, y, data");
				return 1						
			end			
			for i = 3,#data do
				if type(data[i][1]) ~= "table" then
					io.stderr:write("Gnuplot: Format for matrix data: x, y, data");
					return 1						
				end			
			end
			
			for i = 1, #data[1] do
				for j = 1, #data[2] do
					file:write(data[1][i], " ", data[2][j], " ")					
					for k = 3,#data do
						file:write(data[k][i][j], " ")		
					end
					file:write("\n")
				end
				file:write("\n")
			end
		end
		
	-- case 2: data is given column per column
	else
		local plainArray, rowSize
	
		for i, item in ipairs(data) do
			if type(item) == "table" then
				-- check for pure tables
				if plainArray == nil then plainArray = false 
				elseif plainArray == true then
					io.stderr:write("Gnuplot Error: Data array mixes tables and numbers.\n");
					return 1	
				end		
	
				-- check number data items in Row
				if rowSize == nil then rowSize = #item 
				elseif not (rowSize == #item) then
					io.stderr:write("Gnuplot Error: Data array of mixed sized rows.\n");
					return 1	
				end		
				
				-- write data row
				for i, value in ipairs(item) do
					file:write(value, " ");
				end		
				file:write("\n");
			elseif type(item) == "number" then
				-- check for pure numbers
				if plainArray == nil then plainArray = true 
				elseif plainArray == false then
					io.stderr:write("Gnuplot Error: Data array mixes tables and numbers.\n");
					return 1	
				end		
				
				file:write(i, " ", item, "\n");
			else
				io.stderr:write("Gnuplot Error: Data array must contain only tables or only numbers.\n");
				return 1				
			end
		end
	end
		
	-- close file
	file:close()
end

-- example execution
--[[

plaindata = {5, 6, 7, 8}
rowdata = {{3, 6, 4}, {4, 8, 5}, {5, 10, 6}, {6, 6, 7}}
coldata = {{1,2,3,4}, {2,3,4,2}, {3,4,5,3}}

gnuplot.write_data("somedata.txt", data, false)
gnuplot.plot("somedata.pdf", 
			{	
				{file="somedata.txt", 1, 2}
				--{file="somedata.txt", 1, 3} 
			})
			
datasource = {	
	{ label = "catsdogs",  file = "./TestData/dogdata", style = "points", 1, 5 },
	{ label = "catsdogs3", file = "./TestData/dogdata", style = "points", 1, 6 },
	{ label = "catsdogs4", file = "./TestData/dogdata", style = "lines", 1, 5 },
	{ label = "catsdogs5", file = "./TestData/dogdata", style = "lines", 1, 6 },
	{ label = "catsdogs5", file = "./TestData/dogdata", style = "lines", 1, 6 },
	{ label = "sinuskurve",func = "sin(x)", style = "lines"},
	{ label = "rowdata_13", data = rowdata, row = true, style = "lines", 1, 3},
	{ label = "rowdata_12", data = rowdata, row = true, style = "lines"},
	{ label = "col data 12", data = coldata, style = "lines"},
	{ label = "col data 13", data = coldata, style = "lines", 1, 3}
}

-- all options are optional
options = {	title =				"Title", 
			label = 			{x = "xAxis", y = "yAxis"}
			label = 			{"xAxis", "yAxis"}
			xlabel =			"xAxis",
			ylabel =			"yAxis",
			range = 			{ x = {-100, 600}, y = {-5, 5} },
			range = 			{ {-100, 600}, {-5, 5} },
			logscale = 			true,
			logscale = 			{x = true, y = false},
			logscale = 			{true, false},
			dim = 				2,
			grid = 				true,
			key = 				true,
			timestamp =    		true,
			timestamp =     	"My Timestamp string",
			font = 				"Arial",
			fontsize =			10,
			multiplot = 		true,
			multiplotrows = 	5,
			multiplotjoined =	true,
			"additional param #1",
			"additional param #2",
			"additional param #3",
			...
}

gnuplot.plot("vibration.eps", datasource, options)
gnuplot.plot("vibration.svg", datasource, options)
gnuplot.plot("vibration.pdf", datasource, options)
gnuplot.plot("vibration.tex", datasource, options)
gnuplot.plot(nil, datasource, options)

]]--
-- filename		output filename
-- datasource	table of data to be plotted
-- options		table of options
function gnuplot.plot(filename, datasource, options)

	if not datasource or type(datasource) ~= "table" then
		io.stderr:write("Gnuplot Error: a data source dictionary must be provided!\n");
		exit();
	end
	
	----------------------------------------------------------------------------
	-- Check passed options and set defaults
	----------------------------------------------------------------------------

	local options = options or {}
	local title = options.title or "" 
	local label = options.label or false
	local range = options.range or {}
	local padrange = options.padrange or {}
	local logscale = options.logscale or false
	local plotDim = options.dim
	
	local terminal = options.terminal
	local cairo = true
	if type(options.cairo) == "boolean" then cairo = options.cairo end
	local path = options.path or "./" -- output of data

	local size = options.size
	local sizeunit = options.sizeunit or "inch"
	local color = true
	if type(options.color) == "boolean" then color = options.color end
	local enhanced = true
	if type(options.enhanced) == "boolean" then enhanced = options.enhanced end
	local dashed = options.dashed or false
	local font = options.font or "Verdana"
	local fontsize = options.fontsize or 12
	local fontscale = options.fontscale or 1
	local linewidth = options.linewidth or 1
	local linestyle = options.linestyle
	local dashlength = options.dashlength or 1
	local add_term_opt = options.add_term_opt or ""

	local grid = options.grid or false
	local decimalsign = options.decimalsign or "."
	local tics = options.tics
	local mtics = options.mtics or false
	local key = options.key or "on"
	local border = options.border or ""
	
	local slope = options.slope or nil

	local timestamp = options.timestamp or false	
	local multiplot = options.multiplot or false
	local multiplotjoined = options.multiplotjoined or false	
	local multiplotrows = options.multiplotrows

	----------------------------------------------------------------------------
	-- Prepare params
	----------------------------------------------------------------------------
	
	-- check for 2d or 3d data
	for _, source in ipairs(datasource) do
		if #source == 2 then
			if plotDim and plotDim ~= 2 then
				io.stderr:write("Gnuplot Error: Mixed 2d/3d data.\n"); exit();
			end
			plotDim = 2 
		elseif #source == 3 then  
			if plotDim and plotDim ~= 3 then
				io.stderr:write("Gnuplot Error: Mixed 2d/3d data.\n"); exit();
			end
			plotDim = 3 
		elseif #source == 0 then
			-- function  
		else 
			io.stderr:write("Gnuplot Error: pass 0, 2 or 3 columns as data selection.\n"); exit();
		end		
	end
	if plotDim == nil then 
			io.stderr:write("Gnuplot Error: Cannot detect plot dimension.\n"); exit();		
	end
	local DimNames = {"x", "y"}
	if plotDim == 3 then DimNames = {"x", "y", "z"} end
	
	-- set labels
	if label then
		for d, dim in ipairs(DimNames) do
			if not(label[dim]) then 
				label[dim] = label[d] or ""
			end
		end		
	else
		label = {}
		label["x"] = options.xlabel or ""
		label["y"] = options.ylabel or ""
		label["z"] = options.zlabel or ""
	end
	
	
	-- set ranges (otherwise autoscale)
	if type(range) ~= "table" then
		io.stderr:write("Gnuplot Error: range must be a table.");
		exit();
	end
	for d, dim in ipairs(DimNames) do
		if not(range[dim]) then 
			if type(range[d]) == "table" then
				range[dim] = range[d] 
			end
		end
	end
	
	-- set logscales
	if logscale then
		if type(logscale) == "boolean" then
			logscale = {}
			for d, dim in ipairs(DimNames) do
				logscale[dim] = true
			end			
		else
			for d, dim in ipairs(DimNames) do
				if logscale[dim] == nil then 
					logscale[dim] = logscale[d] or false
				end
			end
		end
	end

	-- check for multiplot
	local MultiPlotRows, MultiPlotCols
	if multiplot then
		MultiPlotRows = multiplotrows or math.ceil(math.sqrt(#datasource))
		MultiPlotCols = math.ceil(#datasource / MultiPlotRows )
		fontsize = math.ceil(fontsize / math.max(MultiPlotRows, MultiPlotCols))		
	end

	----------------------------------------------------------------------------
	-- Detect terminal
	----------------------------------------------------------------------------
		
	-- available terminals
	local availTerms = gnuplot.available_terminals()
	local supportedTerms = {"x11", 
							"pdfcairo", "pdf", "epscairo", "postscript",
					   		"pngcairo", "png", "svg", 
					   		"tikz", "cairolatex", "epslatex"}
		
	-- if terminal not given explicitly, try to detect
	if terminal == nil then
		if filename == nil then
			if table.contains(availTerms, "x11") then
				terminal = "x11"	
				-- for interactive, add persitance; size is automatic
				add_term_opt = add_term_opt.." persist raise"
				size = nil
			else
				io.stderr:write("Gnuplot: no terminal for interactive found.\n")
				return 2			
			end		
		else
			local ending = string.sub(filename, -4)
			--- PDF 
			if     ending == ".pdf" then
				if 		table.contains(availTerms, "pdf") 			then terminal = "pdf"	
				elseif 	table.contains(availTerms, "pdfcairo")  	then terminal = "pdfcairo"	end
				if cairo and table.contains(availTerms, "pdfcairo") then terminal = "pdfcairo" end
								
			--- EPS 
			elseif ending == ".eps" then
				if 		table.contains(availTerms, "postscript") 	then terminal = "postscript"	
				elseif 	table.contains(availTerms, "epscairo")  	then terminal = "epscairo"	end
				if cairo and table.contains(availTerms, "epscairo") then terminal = "epscairo" end
				
			--- PNG 
			elseif ending == ".png" then
				if 		table.contains(availTerms, "png") 			then terminal = "png"	
				elseif 	table.contains(availTerms, "pngcairo")  	then terminal = "pngcairo"	end
				if cairo and table.contains(availTerms, "pngcairo") then terminal = "pngcairo" end
				
			--- SVG 
			elseif ending == ".svg" then
				if table.contains(availTerms, "svg")	then terminal = "svg" end
				
			--- TEX 
			elseif ending == ".tex" then
				if table.contains(availTerms, "tikz")			then terminal = "tikz"		
				else
					if table.contains(availTerms, "epslatex") 	then terminal = "epslatex"		
					elseif table.contains(availTerms, "cairolatex")	then terminal = "cairolatex" end		
					if cairo and table.contains(availTerms, "cairolatex") then terminal = "cairolatex" end			
				end
			end	
		end
		
		-- if still nil, error
		if terminal == nil then
			io.stderr:write("Gnuplot Error: cannot deduce terminal for: '"..filename.."'\n")
			io.stderr:write("Gnuplot: Supported endings: pdf, eps, png, svg, tex.\n")
			io.stderr:write("Gnuplot: Supported Terminals: "..table.concat(supportedTerms, ", ").."\n")
			io.stderr:write("Gnuplot: Available Terminals: "..table.concat(availTerms, ", ").."\n")
			return 2			
		end		
	end
		
	-- check valid term
	if not table.contains(supportedTerms, terminal) then
		io.stderr:write("Gnuplot Error: unsupprted terminal: '"..terminal.."'\n")
		io.stderr:write("Gnuplot Error: supported are: "..table.concat(supportedTerms, ", ").."\n")
		return 2				
	end	

	----------------------------------------------------------------------------
	-- terminal options
	----------------------------------------------------------------------------	

	local term = {}

	-- size of canvas
	-- Note: 1 inch = 72 pt = 2.54 cm
	local dpi = 72
	local inchINcm = 2.54
	term.size = "" -- on support: x11
	if size then
		if type(size) ~= "table" or #size ~= 2 then
			io.stderr:write("Gnuplot: specify size with table of {<x>,<y>}\n")
			return 2				
		end
	
		-- for vector graphics use cm 
		if table.contains({"pdfcairo", "pdf", "epscairo", "postscript",
			   			   "tikz","cairolatex","epslatex"}, terminal) 
		then
			local s = {size[1], size[2]}
			if sizeunit == "cm" then
				-- we express it in cm
			elseif sizeunit == "mm" then
				s[1] = s[1] * 0.1
				s[2] = s[2] * 0.1
			elseif sizeunit == "inch" or sizeunit == "in" then
				s[1] = s[1] * inchINcm
				s[2] = s[2] * inchINcm
			elseif sizeunit == "pt" or sizeunit == "pixel" then
				s[1] = s[1] * (inchINcm/dpi)			
				s[2] = s[2] * (inchINcm/dpi)			
			else
				io.stderr:write("Gnuplot: sizeunit invalid. use: cm, mm, in, inch, pt, pixel\n")
				return 2				
			end	
			
			term.size = "size "..s[1].."cm,"..s[2].."cm"
		end	
		
		-- for pictures (non-vector) use pixel
		if table.contains({"pngcairo", "png", "svg"}, terminal)
		then
			local s = {size[1], size[2]}
			if sizeunit == "pt" or sizeunit == "pixel" then
				-- we express it in pt
			elseif sizeunit == "inch" or sizeunit == "in" then
				s[1] = s[1] * dpi
				s[2] = s[2] * dpi
			elseif sizeunit == "cm" then
				s[1] = s[1] * (dpi/inchINcm)			
				s[2] = s[2] * (dpi/inchINcm)			
			elseif sizeunit == "mm" then
				s[1] = s[1] * (dpi/inchINcm) * 0.1
				s[2] = s[2] * (dpi/inchINcm) * 0.1
			else
				io.stderr:write("Gnuplot: sizeunit invalid. use: cm, mm, in, inch, pt, pixel\n")
				return 2				
			end		
			
			term.size = "size "..s[1]..","..s[2]
		end		
	end
	
	-- { enhanced |Ênoenhanced }
	term.enhanced = "" -- no support: tikz, cairolatex, epslatex
	if table.contains({"pdfcairo", "pdf", "epscairo", "postscript",
			   			"pngcairo", "png", "svg", "x11"}, terminal) then
		if enhanced then term.enhanced = "enhanced"
		else 			 term.enhanced = "noenhanced" end
	end
	
	-- {color |Êmono}
	term.color = "" -- no support: svg, png 
	if color then
		if table.contains({"pdfcairo","epscairo","pngcairo", 
							"pdf","postscript","tikz","cairolatex","epslatex"}, terminal) then
			term.color = "color"
		end
	else
		if table.contains({"pdfcairo","epscairo","pngcairo","cairolatex"}, terminal) then
			term.color = "mono"	
		elseif table.contains({"pdf","postscript","tikz","epslatex"}, terminal) then
			term.color = "monochrome"
		end
	end
	
	-- font + fontsize + fontscale
	-- NOTE: for tex, the font should be controled within latex
	
	-- NOTE: default fontscale is 0.5 for {pdf/eps}cairo [who invented that ...!? crazy!]
	if table.contains({"pdfcairo","epscairo"}, terminal) 	then
		fontscale = fontscale * 0.5
	end
	
	term.font = "" -- no support: tikz, cairolatex, epslatex
	if table.contains({"pdfcairo","epscairo","pngcairo", 
						"pdf","postscript","png"}, terminal) 
	then
		term.font = "font '"..font..","..fontsize.."' fontscale "..fontscale	
	elseif table.contains({"svg", "x11"}, terminal) then
		term.font = "font '"..font..","..fontsize.."'"	
	end	
	
	-- { solid |Êdashed}
	term.dashed = "" -- no support: svg, png 
	if table.contains({"pdfcairo","epscairo","pngcairo", 
						"pdf","postscript","tikz","cairolatex","epslatex", "x11"}, terminal) then
		if dashed or not color then
			term.dashed = "dashed"
		else
			term.dashed = "solid"
		end
	end
	
	-- { rounded |Êbutt}
	-- currently not supported, since not needed - butt seems always good
	
	-- linewidth
	term.linewidth = "" -- no support: 
	if table.contains({"pdfcairo","epscairo","pngcairo", 
						"pdf","postscript","png","svg","cairolatex","epslatex", "x11"}, terminal) 
	then
		term.linewidth = "linewidth "..linewidth
	end	
	
	-- dashlength
	term.dashlength = "" -- no support: svg
	if table.contains({"pdfcairo","epscairo","pngcairo", 
						"postscript","png","cairolatex","epslatex"}, terminal) 
	then
		term.dashlength = "dashlength "..dashlength
	elseif table.contains({"pdf"}, terminal) then
		term.dashlength = "dl "..dashlength
	end	
	
	-- fix postscript (that is actually only used for eps)
	if terminal == "postscript" then terminal = "postscript eps" end
		
	----------------------------------------------------------------------------
	-- Write script header
	----------------------------------------------------------------------------

	if filename then filename = string.gsub(tostring(filename), " ", "_" )
	else filename = "__name__" end
	
	local tmpPath = path.."tmp-gnuplot/" -- some tmp path
	if not(DirectoryExists(tmpPath)) then CreateDirectory(tmpPath) end
	local scriptName = tmpPath.."tmp_gnuplot_script_"..string.gsub(filename, "[./]", "_")..".gnu"
	local script = io.open(scriptName, "w+")
	if not script then
		io.stderr:write("Gnuplot Error: cannot open output file: '")
		io.stderr:write(scriptName .. " '\n");
		return 2
	end

	-- start the output file
	script:write("reset\n")
		
	-- set terminal currently only pdf
	script:write("set term "..terminal)
	script:write(" "..term.size)
	script:write(" "..term.enhanced)
	script:write(" "..term.color)
	script:write(" "..term.dashed)
	script:write(" "..term.font)
	script:write(" "..term.linewidth)
	script:write(" "..term.dashlength)
	script:write(" "..add_term_opt)
	
	script:write("\n")
		
	if terminal ~= "x11" then
		script:write("set output \"", path, filename,"\" \n")
	end

	-- encoding and decimalsign
	script:write("\nset encoding utf8\n")
	script:write("set decimalsign '"..decimalsign.."'\n")
	
		
	-- title and axis label
	script:write("set title '"..title.."'\n")
	
	-- labels
	for _, dim in ipairs(DimNames) do
		script:write("set "..dim.."label '"..label[dim].."'\n")
	end
	
	-- write timestamp
	if timestamp then
		if type(timestamp) == "string" then
			script:write("set timestamp '"..timestamp.."'\n")
		else
			script:write("set timestamp \"%H:%M:%S  %Y/%m/%d\" bottom\n")
		end
	else 
		script:write("unset timestamp\n"); 
	end

	-- tics
	script:write("set tics\n");
	if tics ~= nil then
		if type(tics) == "boolean" and tics == false then
			script:write("unset tics\n"); 		
		elseif type(tics) == "string" then 
			script:write("set tics "..tics.."\n"); 		
		elseif type(tics) == "table" then
			for dim, dimTic in pairs(tics) do
				if type(dimTic) == "boolean" and dimTic == false then 
					script:write("unset "..dim.."tics\n"); 		
				elseif type(dimTic) == "string" then 
					script:write("set "..dim.."tics "..dimTic.."\n"); 		
				end
			end
		end
	end

	-- mtics
	if mtics then
		if type(mtics) == "boolean" and mtics == true then
			for _, dim in ipairs(DimNames) do
					script:write("set m"..dim.."tics default\n"); 		
			end
		elseif type(mtics) == "number" then 
			for _, dim in ipairs(DimNames) do
					script:write("set m"..dim.."tics "..mtics.."\n"); 		
			end
		elseif type(mtics) == "table" then
			for dim, dimTic in pairs(mtics) do
				if type(dimTic) == "boolean" then
				 	if dimTic == false then script:write("unset m"..dim.."tics\n"); 		
				 	else 					script:write("set m"..dim.."tics default\n"); 
				 	end
				elseif type(dimTic) == "number" then 
					script:write("set m"..dim.."tics "..dimTic.."\n"); 		
				end
			end
		end
	else
		for _, dim in ipairs(DimNames) do
				script:write("unset m"..dim.."tics\n"); 		
		end
	end
	
	-- enable grid
	if grid then
		script:write("set grid ")
		for _, dim in ipairs(DimNames) do
			script:write(dim.."tics m"..dim.."tics "); 
		end
		if type(grid) == "string" then
			script:write(grid); 
		end
		script:write(" back\n"); 
	end

	-- logscale
	if logscale then
		for _, dim in ipairs(DimNames) do
			if logscale[dim] then
				script:write("set logscale "..dim.."\n");
			end
		end
	end

	-- set default linetypes
	if linestyle and linestyle.colors then
		local colors = linestyle.colors
		local linewidth = linestyle.linewidth or 1
		local pointsize = linestyle.pointsize or 1
		for i=1,#colors do
		 	script:write("set linetype "..i.." lc rgb \""..colors[i].."\" "
		 					.."lw "..linewidth.." ps "..pointsize.."\n")
		end		
		script:write("set linetype cycle "..#colors.."\n")		
	end

	if border then			
		script:write("set border "..border.."\n")
	end
	
	-- write additional options, passed as strings
	for _, userParam in ipairs(options) do
		if userParam and type(userParam) == "string" then
			script:write(userParam.."\n")
		end
	end	

	-- multiplot sizes
	if multiplot then
		script:write("set multiplot layout ", MultiPlotRows,", ", MultiPlotCols)
		script:write(" title '"..title.."' \n\n" )	
	end

	----------------------------------------------------------------------------
	-- Detect data sets sizes
	----------------------------------------------------------------------------

	-- find data ranges
	local stats = {}
	for s, source in ipairs(datasource) do
	
		stats[s] = {}
		local stat = stats[s]
		
		-- get the data column mapping
		stat.map = {};
		local map = stat.map
		if     #source == 2 then map[1] = source[1]; map[2] = source[2];
		elseif #source == 3 then map[1] = source[1]; map[2] = source[2]; map[3] = source[3];
		else
			if plotDim == 2 then map[1] = 1; map[2] = 2;
			else 				 map[1] = 1; map[2] = 2; map[3] = 3; end
		end

		-- extract values if in file
		if source.file then	
			stat.val = {}
			local val = stat.val
			for line in io.lines(source.file) do
				local numbers = {}
				if not string.match(line, "#") then
					for n in string.gmatch(line, "[^ ]+") do numbers[#numbers+1] = tonumber(n) end
	
					for d = 1,#map do 
						val[d] = val[d] or {}
						val[d][#(val[d])+1] = numbers[map[d]]
					end
				end
			end
		end

		-- extract values if in array
		if source.data then
			stat.val = {}
			if source.row then
				for row in 1,#source.data do 
					for d = 1,#map do 
						val[d] = val[d] or {}
						val[d][#(val[d])+1] = source.data[row][map[d]]
					end
				end
			else
				for row in 1,#source.data[1] do 
					for d = 1,#map do 
						val[d] = val[d] or {}
						val[d][#(val[d])+1] = source.data[map[d]][row]
					end
				end			
			end
		end

		-- find min/max
		stat.min = {}
		stat.max = {}
		for d = 1,#stat.val do
			stat.min[d] = math.huge
			stat.max[d] = -math.huge
			for n = 1,#stat.val[d] do
				stat.min[d] = math.min(stat.min[d], stat.val[d][n])	
				stat.max[d] = math.max(stat.max[d], stat.val[d][n])	
			end		
		end
		
		-- global min/max
		if not stats.min or not stats.max then
			stats.min = {}
			stats.max = {}
			stats.range = {}
			for d = 1,#stat.min do
				stats.min[d] = math.huge
				stats.max[d] = -math.huge
			end
		end
		for d = 1,#stat.val do
			stats.min[d] = math.min(stats.min[d], stat.min[d])	
			stats.max[d] = math.max(stats.max[d], stat.max[d])	
			stats.range[d] = stats.max[d] - stats.min[d]
		end
	end

	-- add slope triangle to every source
	for s, source in ipairs(datasource) do
		local function drawSlopTri(xo, yo, dy, p)
--			local dy = dx^p
			local dx = dy^(1/p)
			script:write("set arrow from "..xo..","..yo.." rto "..dx..","..dy.." nohead lw 2\n") -- slope
			script:write("set arrow from "..xo..","..yo*dy.." rto "..dx..", 1 nohead lw 2\n") -- waagerecht
			script:write("set arrow from "..xo..","..yo.." rto 1,"..dy.." nohead lw 2\n") -- vertical
			script:write("set label right '"..p.."' at "..(xo*1.05)..","..math.pow(dy, 2/3)*yo.. "font ',8'\n") -- label
			script:write("set label center '"..(1).."' at "..math.sqrt(dx)*xo..","..(yo*dy*1.8).. "font ',8'\n") -- label
		end
		
		if slope then
			drawSlopTri(stats[s].min[1], stats[s].min[2]*1.5, slope.dx, s)							
		end
	end
	
	-- add key (legend)
	if key then script:write ("set key "..key.."\n") end

	-- range
	script:write ("set autoscale\n")
	for _, dim in ipairs({"x", "y", "z"}) do	
		if range[dim] then
			script:write ("set "..dim.."range [",range[dim][1],":",range[dim][2],"]\n")
		end
	end
	for d, dim in ipairs(DimNames) do	
		if padrange[dim] then
			script:write("set "..dim.."range [")
			script:write(stats.min[d]*padrange[dim][1])
			script:write(":")
			script:write(stats.max[d]*padrange[dim][2])
			script:write("]\n")
		end
	end

	----------------------------------------------------------------------------
	-- Write data sets
	----------------------------------------------------------------------------

	for s, source in ipairs(datasource) do
	
		-- get the data column mapping
		local map = table.concat(stats[s].map, ":")

		-- special 3d treatment
		if plotDim == 3 then
			script:write("set surface\n") 
			script:write("unset contour\n") 
		end
		
		-- determine the plot style - data source table has priority
		local style = source.style
		if not style then
			if plotDim == 3 then style = "points"
			else				 style = "linespoints" end
		end
		
		-- check style		
		if not table.contains({"lines","points","linespoints","boxes","dots","vectors","yerrorbars"}, style) then
			io.stderr:write("Gnuplot Error: style=\""..style.."\" not supported.\n");
			exit()
		end
		
		-- assign default data label if needed
		if not source.label then source["label"] = "plot" .. s end

		-- Position plots in multiplot
		if multiplot then
			script:write("unset title \n" )

			if multiplotjoined then
				local spaceTop = 0.1
				local spaceBottom = 0.1
				local spaceLeft = 0.1
				local SpaceRight = 0.1
				
				local rowSize = (1-spaceBottom-spaceTop)/MultiPlotRows
				local colSize = (1-spaceLeft-SpaceRight)/MultiPlotCols
				
				local mpCol = ((s-1) % MultiPlotCols) + 1
				local mpRow = math.ceil(s / MultiPlotCols)
				
				script:write("set tmargin at screen "..spaceBottom + (MultiPlotRows-mpRow+1)*rowSize.."\n")
				script:write("set bmargin at screen "..spaceBottom + (MultiPlotRows-mpRow)*rowSize.."\n")
				script:write("set rmargin at screen "..spaceLeft + (mpCol)*colSize.."\n")
				script:write("set lmargin at screen "..spaceLeft + (mpCol-1)*colSize.."\n")
				
				-- left bnd
				if mpCol == 1 then 
					script:write("set ylabel '"..label["y"].."'\n\n")
					script:write("set format y\n")
				else 
					script:write("unset ylabel\n")
					script:write("set format y \"\"\n")
				end

				-- right bnd
				if mpCol == MultiPlotCols then end
				
				-- top bnd
				if mpRow == 1 then end

				-- bottom bnd
				if mpRow == MultiPlotRows then 
					script:write("set xlabel '"..label["x"].."'\n\n")
					script:write("set format x\n")
				else 
					script:write("unset xlabel\n")
					script:write("set format x \"\"\n")
				end
			end
		end

		-- build up the plot command, layer by layer in non-multiplotmode
		if s == 1 or multiplot then 	
			if plotDim == 2 then script:write ("plot  \\\n");
			else 				 script:write ("splot  \\\n"); end
		else 				
			script:write(", \\\n");	
		end
		
		if source.file then		
			script:write ("\"", source.file, "\"",
							  " using ", map, 
							  " title '", source.label, "'", 
							  " with ", style)
		end
		if source.func then
			script:write (" ", source.func," ",
							  " title '", source.label, "'", 
							  " with ", style)
		end 
		if source.data then
			tmpDataFile = tmpPath.."tmp_"..source.label.."_"..s..".dat"
			tmpDataFile = string.gsub(tmpDataFile, " ", "_" )
			gnuplot.write_data(tmpDataFile, source.data, source.row)
			script:write ("\"", tmpDataFile, "\"",
							  " using ", map, 
							  " title '", source.label, "'", 
							  " with ", style)
		end 
								  
		if multiplot then script:write("\n") end						  
	end
	
	-- finish script
	script:write ("\nunset multiplot\n");
	script:close()

	-- launch gnuplot and run in background
	os.execute("gnuplot "..scriptName.." --persist &")

	return 0
end


-- return a new array containing the concatenation of all of its 
-- parameters. Scaler parameters are included in place, and array 
-- parameters have their values shallow-copied to the final array.
-- Note that userdata and function values are treated as scalar.
function gnuplot.array_concat(...) 
    local t = {}
    for n = 1,select("#",...) do
        local arg = select(n,...)
        if type(arg)=="table" then
            for _,v in ipairs(arg) do
                t[#t+1] = v
            end
        else
            t[#t+1] = arg
        end
    end
    return t
end

--[[!
\}
]]--
