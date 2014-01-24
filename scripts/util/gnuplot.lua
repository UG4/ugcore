--[[!
\file gnuplot.lua
\defgroup scripts_util_gnuplot GNUPlot Utility
create a globally seen package (all non-local functions call then be called using gnuplot.xxx.)
\ingroup scripts_util
\{
]]--

gnuplot = gnuplot or {}

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
	local logscale = options.logscale or false
	local plotDim = options.dim
	
	local terminal = options.terminal
	local cairo = true
	if type(options.cairo) == "boolean" then cairo = options.cairo end
	local path = options.path or "./" -- output of data

	local size = options.size
	local color = true
	if type(options.color) == "boolean" then color = options.color end
	local enhanced = true
	if type(options.enhanced) == "boolean" then enhanced = options.enhanced end
	local dashed = options.dashed or false
	local font = options.font or "Verdana"
	local fontsize = options.fontsize or 12
	local fontscale = options.fontscale or 1
	local linewidth = options.linewidth or 1
	local dashlength = options.dashlength or 1
	local add_term_opt = options.add_term_opt or ""

	local grid = options.grid or false
	local key = true
	if type(options.key) == "boolean" then key = options.key end

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
				-- for interactive, add persitance
				add_term_opt = add_term_opt.." persist raise"
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
				elseif table.contains(availTerms, "epslatex") 	then terminal = "epslatex"		
				elseif table.contains(availTerms, "cairolatex")	then terminal = "cairolatex" end		
				if cairo and table.contains(availTerms, "cairolatex") then terminal = "cairolatex" end			
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
	

	local term = {}
	
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
	if size then script:write(" size "..size[1]..","..size[2]) end
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
		
	-- title and axis label
	script:write("\nset encoding utf8\n")
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

	-- enable grid
	if grid then
		script:write("set grid ")
		for _, dim in ipairs(DimNames) do
			script:write(dim.."tics "); 
		end
		script:write("back\n"); 
	end

	-- logscale
	if logscale then
		for _, dim in ipairs(DimNames) do
			if logscale[dim] then
				script:write("set logscale "..dim.."\n");
			end
		end
	end
	
	-- remove key (legend)
	if not(key) then script:write ("unset key\n") end

	-- range
	script:write ("set autoscale\n")
	for _, dim in ipairs({"x", "y", "z"}) do	
		if range[dim] then
			script:write ("set "..dim.."range [",range[dim][1],":",range[dim][2],"]\n")
		end
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
	-- Write data sets
	----------------------------------------------------------------------------

	for sourceNr, source in ipairs(datasource) do
	
		-- get the data column mapping
		local map = "";
		if     #source == 2 then map = source[1]..":"..source[2]
		elseif #source == 3 then map = source[1]..":"..source[2]..":"..source[3]
		else
			if plotDim == 2 then map = "1:2"
			else map = "1:2:3" end
		end

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
		if not source.label then source["label"] = "plot" .. sourceNr end

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
				
				local mpCol = ((sourceNr-1) % MultiPlotCols) + 1
				local mpRow = math.ceil(sourceNr / MultiPlotCols)
				
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
		if sourceNr == 1 or multiplot then 	
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
			tmpDataFile = tmpPath.."tmp_"..source.label.."_"..sourceNr..".dat"
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
