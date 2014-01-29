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
			
plots = {	
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

gnuplot.plot("vibration.eps", plots, options)
gnuplot.plot("vibration.svg", plots, options)
gnuplot.plot("vibration.pdf", plots, options)
gnuplot.plot("vibration.tex", plots, options)
gnuplot.plot(nil, plots, options)

]]--
-- filename		output filename
-- plots	table of data to be plotted
-- options		table of options
function gnuplot.plot(filename, data, options)

	local options = options or {}

	----------------------------------------------------------------------------
	-- Detect type of data
	----------------------------------------------------------------------------
	
	-- returns if a table item is a correct data item
	-- a dataset is 
	local function IsDataSet(dataset)
		-- must be a table
		if type(dataset) ~= "table" then return false end 
		
		-- must contain exactly one of: file, func, data
		if 	(dataset.file and dataset.func) or 
			(dataset.file and dataset.data) or
			(dataset.func and dataset.data) then return false end
			
		-- must not contain table in integral keys
		for _, v in ipairs(dataset) do
			if type(v) == "table" then return false end
		end
			
		-- check for source-specification
		if 	   dataset.file then
			return true
		elseif dataset.func then
			-- no column specs for function
			if #dataset == 0 then return true else return false end		
		elseif dataset.data then
			return true
		else 		
			-- not contained a source-spec
			return false
		end
	end
	
	-- prepare plot array 
	-- we allow three cases:
	-- a) table depth 1: single data: data = {file = "f1.dat", 1, 2}	
	-- b) table depth 2: multi-data:  data = {{file = "f1.dat", 1, 2}, {file = "f1.dat", 1, 2}, ..}	
	-- c) table depth 3: multiplot-multi-data:  
	--			data = { { {file = "f1.dat", 1, 2}, {file = "f1.dat", 1, 2}, ...}, 
	--					 { {file = "f3.dat", 1, 2}, {file = "f3.dat", 1, 2}, ...}, ...}	
	-- NOTE: plots will always be of structure c), i.e.
	--	plots = { { DataSet, DataSet, ...}, { DataSet, DataSet, ...}, ... }
	-- where DataSet is of type { [file | data |Êfunc ] = , ...}
	local plots = nil
	local multiplot = nil -- will be set based on data
	if data then
		-- a)
		if IsDataSet(data) then
			plots = {{table.deepcopy(data)}}
			multiplot = false
		elseif type(data) == "table" then 
		 	local bMulti = true
		 	for _,dataset in ipairs(data) do
		 		if not IsDataSet(dataset) then bMulti = false end
		 	end
		 	
	 		-- b)
		 	if bMulti then
				if options.multiplot ~= nil and options.multiplot then
					-- create several plots
					plots = {}
				 	for s,dataset in ipairs(data) do
				 		plots[s] = {table.deepcopy(data[s])}
					end
					multiplot = options.multiplot
				else
					-- default is no multiplot
					plots = {table.deepcopy(data)}
					multiplot = false
				end		
		 	else
			 	local bMultiMulti = true
			 	for _,plot in ipairs(data) do
			 		if type(plot) ~= "table" then bMultiMulti = false 
			 		else
					 	for _,dataset in ipairs(plot) do
					 		if not IsDataSet(dataset) then bMultiMulti = false end
					 	end
					end
				end
				
				-- c) 
				if bMultiMulti then 
					if type(options.multiplot) == "boolean" and not options.multiplot then
						-- create one plot
						plots = {{}}
					 	for _,plot in ipairs(data) do
						 	for _,dataset in ipairs(plot) do
					 			plots[1][#plots[1]+1] = table.deepcopy(dataset)
					 		end
						end
						multiplot = false
					else
						-- default is multiplot
						plots = table.deepcopy(data)
						multiplot = options.multiplot or true
					end		
				end
		 	end
		end
	end

	if not plots then
		io.stderr:write("Gnuplot Error: data source not correctly specifiec. Valid formats are:\n");
		io.stderr:write(" a) single data: data = {file = \"f1.dat\", 1, 2}\n");
		io.stderr:write(" b) multi-data:  data = {{file = \"f1.dat\", 1, 2}, {file = \"f1.dat\", 1, 2}, ...}	\n");
		io.stderr:write(" c) multiplot-multi-data: \n");		
		io.stderr:write("    data = { { {file = \"f1.dat\", 1, 2}, {file = \"f1.dat\", 1, 2}, ...}, \n");		
		io.stderr:write("             { {file = \"f3.dat\", 1, 2}, {file = \"f3.dat\", 1, 2}, ...}, ...} \n");		
		exit();
	end
	
	----------------------------------------------------------------------------
	-- Terminal options
	----------------------------------------------------------------------------
	
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

	local decimalsign = options.decimalsign or "."
	local timestamp = options.timestamp or false		
	
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
	term.sizes = {12.5, 8.75}
	term.sizeunit = "cm"
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
			
			term.sizes = {s[1], s[2]}
			term.sizeunit = "cm"
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
			term.sizes = {s[1], s[2]}
			term.sizeunit = ""
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
	-- Write terminal
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
	
	----------------------------------------------------------------------------
	-- Detect Plot Dimension
	----------------------------------------------------------------------------
	
	-- check for 2d or 3d data
	local plotDim = options.dim
	for _, plot in ipairs(plots) do
		for _, dataset in ipairs(plot) do
			if #dataset == 2 then
				if plotDim and plotDim ~= 2 then
					io.stderr:write("Gnuplot Error: Mixed 2d/3d data.\n"); exit();
				end
				plotDim = 2 
			elseif #dataset == 3 then  
				if plotDim and plotDim ~= 3 then
					io.stderr:write("Gnuplot Error: Mixed 2d/3d data.\n"); exit();
				end
				plotDim = 3 
			elseif #dataset == 0 then
				-- function  
			else 
				io.stderr:write("Gnuplot Error: pass 0, 2 or 3 columns as data selection.\n"); exit();
			end		
		end
	end
	if (plotDim ~= 2 and plotDim ~= 3) then 
		io.stderr:write("Gnuplot Error: Cannot detect plot dimension (only 2d or 3d).\n"); exit();		
	end
	
	-- dim fields
	local DimNames = {"x", "y"}
	if plotDim == 3 then DimNames = {"x", "y", "z"} end

	----------------------------------------------------------------------------
	-- Detect data sets sizes
	----------------------------------------------------------------------------

	-- find data ranges
	for _, plot in ipairs(plots) do
		for _, dataset in ipairs(plot) do
	
			-- get the data column mapping
			dataset.map = {}
			local map = dataset.map
			if #dataset == #DimNames then 
				for d, _ in ipairs(DimNames) do map[d] = dataset[d]; end
			else
				for d, _ in ipairs(DimNames) do map[d] = d; end
			end

			dataset.val = {}
			local val = dataset.val
	
			-- extract values if in file
			if dataset.file then	
				for line in io.lines(dataset.file) do
					local numbers = {}
					if not string.match(line, "#") then
						for n in string.gmatch(line, "[^ ]+") do numbers[#numbers+1] = tonumber(n) end
		
						for d, _ in ipairs(DimNames) do 
							val[d] = val[d] or {}
							val[d][#(val[d])+1] = numbers[map[d]]
						end
					end
				end
			end
	
			-- extract values if in array
			if dataset.data then
				if dataset.row then
					for row = 1, #(dataset.data) do 
						for d, _ in ipairs(DimNames) do 
							val[d] = val[d] or {}
							val[d][#(val[d])+1] = dataset.data[row][map[d]]
						end
					end
				else
					for row = 1, #dataset.data[1] do 
						for d, _ in ipairs(DimNames) do 
							val[d] = val[d] or {}
							val[d][#(val[d])+1] = dataset.data[map[d]][row]
						end
					end			
				end
			end
	
			-- find min/max
			dataset.min = {}
			dataset.max = {}
			for d, _ in ipairs(DimNames) do
				dataset.min[d] = math.huge
				dataset.max[d] = -math.huge
				for n = 1,#dataset.val[d] do
					dataset.min[d] = math.min(dataset.min[d], dataset.val[d][n])	
					dataset.max[d] = math.max(dataset.max[d], dataset.val[d][n])	
				end		
			end
		end
	end
	
	-- global min/max
	local stats = {}
	for _, plot in ipairs(plots) do
		for _, dataset in ipairs(plot) do
			if not stats.min or not stats.max then
				stats.min = {}
				stats.max = {}
				stats.range = {}
				for d, _ in ipairs(DimNames) do
					stats.min[d] = math.huge
					stats.max[d] = -math.huge
				end
			end
			for d, _ in ipairs(DimNames) do
				stats.min[d] = math.min(stats.min[d], dataset.min[d])	
				stats.max[d] = math.max(stats.max[d], dataset.max[d])	
				stats.range[d] = stats.max[d] - stats.min[d]
			end
		end
	end

	----------------------------------------------------------------------------
	-- special treatments
	----------------------------------------------------------------------------
	
	-- special 3d treatment
	if plotDim == 3 then
		script:write("set surface\n") 
		script:write("unset contour\n") 
	end
	
-- multiplot sizes
	local MultiPlotRows , MultiPlotCols
	if multiplot then
		if type(multiplot) == "table" then
			MultiPlotRows = multiplot.rows or math.ceil(math.sqrt(#plots))			
		end
		MultiPlotRows = MultiPlotRows or math.ceil(math.sqrt(#plots))		
		MultiPlotCols = math.ceil(#plots / MultiPlotRows )
		
		script:write("set term "..terminal.." size ")
		script:write(MultiPlotRows*term.sizes[1]..term.sizeunit)
		script:write(",")
		script:write(MultiPlotCols*term.sizes[2]..term.sizeunit)
		script:write("\n")
		
		local title = options.title or "" 
		script:write("set multiplot layout ", MultiPlotRows,", ", MultiPlotCols)
		script:write(" title '"..title.."' \n\n" )	
	end

	----------------------------------------------------------------------------
	-- loop plots
	----------------------------------------------------------------------------

	for p, plot in ipairs(plots) do

		------------------------------------------------------------------------
		-- Plot options
		------------------------------------------------------------------------
		
		local title = plot.title or options.title or "" 
		
		local label = plot.label or options.label or false
		for _, dim in ipairs(DimNames) do label[dim] = label[dim] or "" end
		
		local range = plot.range or options.range or {}
		local padrange = plot.padrange or options.padrange or {}
		
		local logscale = plot.logscale or options.logscale or false

		local tics = plot.tics or options.tics or true
		if tics and type(plot.tics) == "boolean" and not plot.tics then tics = false end
		if tics and type(options.tics) == "boolean" and not options.tics then tics = false end
		
		local mtics = plot.mtics or options.mtics or false
		local grid = plot.grid or options.grid or false
		
		local key = plot.key or options.key or "on"
		local border = plot.border or options.border or ""	
		
		local slope = plot.slope or options.slope or nil
		
		------------------------------------------------------------------------
		-- Write Plot options
		------------------------------------------------------------------------
					
		-- title
		script:write("set title '"..title.."'\n")
			
		-- labels
		script:write("unset label\n")
		for _, dim in ipairs(DimNames) do
			script:write("set "..dim.."label '"..label[dim].."'\n")
		end
		
		-- range
		script:write ("set autoscale\n")
		for _, dim in ipairs(DimNames) do	
			if range[dim] then
				script:write ("set "..dim.."range [",range[dim][1],":",range[dim][2],"]\n")
			end
		end
		
		-- padrange
		for d, dim in ipairs(DimNames) do	
			if padrange[dim] then
				script:write("set "..dim.."range [")
				script:write(stats.min[d]*padrange[dim][1])
				script:write(":")
				script:write(stats.max[d]*padrange[dim][2])
				script:write("]\n")
			end
		end

		-- logscale
		script:write("unset logscale\n");
		if logscale then
			if type(logscale) == "boolean" then
				logscale = {}
				for _, dim in ipairs(DimNames) do logscale[dim] = true end			
			end
			
			for _, dim in ipairs(DimNames) do
				if logscale[dim] then
					script:write("set logscale "..dim.."\n");
				end
			end
		end
	
		-- tics
		local sTic = {}
		for _, dim in ipairs(DimNames) do 
			sTic[dim] = "set "..dim.."tics border mirror norotate\n"
		end
		
		if type(tics) == "boolean" and not tics then
			for _, dim in ipairs(DimNames) do 
				sTic[dim] = "unset "..dim.."tics\n"
			end
		elseif type(tics) == "string" then 
			for _, dim in ipairs(DimNames) do 
				sTic[dim] = "set "..dim.."tics "..tics.."\n" 
			end 		
		elseif type(tics) == "table" then
			for dim, dimTic in pairs(tics) do
				if type(dimTic) == "boolean" and not dimTic then 
					sTic[dim] = "unset "..dim.."tics\n" 		
				elseif type(dimTic) == "string" then 
					sTic[dim] = "set "..dim.."tics "..dimTic.."\n" 		
				end
			end
		end
		for _, dim in ipairs(DimNames) do
			script:write(sTic[dim]); 		
		end
			
		-- mtics
		for _, dim in ipairs(DimNames) do
			script:write("unset m"..dim.."tics\n"); 		
		end
		if mtics then
			if type(mtics) == "boolean" then
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
		end
		
		-- grid
		if grid then
			script:write("set grid ")
			for _, dim in ipairs(DimNames) do
				script:write(dim.."tics m"..dim.."tics "); 
			end
			if type(grid) == "string" then
				script:write(grid); 
			end
			script:write(" back\n"); 
		else
			script:write("unset grid \n")
		end
	
		-- default linetypes
		if linestyle then
			local colors = linestyle.colors or {}
			local linewidth = linestyle.linewidth or 1
			local pointsize = linestyle.pointsize or 1
			for i=1,#colors do
			 	script:write("set linetype "..i.." lc rgb \""..colors[i].."\" "
			 					.."lw "..linewidth.." ps "..pointsize.."\n")
			end		
			script:write("set linetype cycle "..#colors.."\n")		
		end
	
		-- border
		script:write("set border 31\n")
		if border then			
			script:write("set border "..border.."\n")
		end
		
		-- write additional options, passed as strings
		for _, userParam in ipairs(options) do
			if userParam and type(userParam) == "string" then
				script:write(userParam.."\n")
			end
		end	
	
		-- add key (legend)
		script:write("unset key\n")
		if key then 
			script:write("set key "..key.."\n") 
		end
	
		------------------------------------------------------------------------
		-- Position plots in multiplot
		------------------------------------------------------------------------
		if multiplot then

			if type(multiplot) == "table" and multiplot.conjoined then
				
				-- only one global title --> unset plot title
				script:write("unset title \n" )
			
				local spaceTop = 0.1
				local spaceBottom = 0.1
				local spaceLeft = 0.1
				local SpaceRight = 0.1
				
				local rowSize = (1-spaceBottom-spaceTop)/MultiPlotRows
				local colSize = (1-spaceLeft-SpaceRight)/MultiPlotCols
				
				local mpCol = ((p-1) % MultiPlotCols) + 1
				local mpRow = math.ceil(p / MultiPlotCols)
				
				script:write("set tmargin at screen "..spaceBottom + (MultiPlotRows-mpRow+1)*rowSize.."\n")
				script:write("set bmargin at screen "..spaceBottom + (MultiPlotRows-mpRow)*rowSize.."\n")
				script:write("set rmargin at screen "..spaceLeft + (mpCol)*colSize.."\n")
				script:write("set lmargin at screen "..spaceLeft + (mpCol-1)*colSize.."\n")
				
				-- left bnd
				if mpCol == 1 then 
					script:write("set ylabel '"..label["y"].."'\n")
					script:write(sTic["y"])
				else 
					script:write("unset ylabel\n")
					script:write("set ytics format ''\n")
				end

				-- right bnd
				if mpCol == MultiPlotCols then end
				
				-- top bnd
				if mpRow == 1 then end

				-- bottom bnd
				if mpRow == MultiPlotRows then 
					script:write("set xlabel '"..label["x"].."'\n")
					script:write(sTic["x"])
				else 
					script:write("unset xlabel\n")
					script:write("set xtics format ''\n")
				end
			end
		end

		------------------------------------------------------------------------
		-- add slope triangle to every source
		------------------------------------------------------------------------
		if slope then
			local function drawSlopTri(xo, yo, dy, p)
				local dx = dy^(1/p)
				script:write("set arrow from "..xo..","..yo.." rto "..dx..","..dy.." nohead lw 2\n") -- slope
				script:write("set arrow from "..xo..","..yo*dy.." rto "..dx..", 1 nohead lw 2\n") -- waagerecht
				script:write("set arrow from "..xo..","..yo.." rto 1,"..dy.." nohead lw 2\n") -- vertical
				local align = "right"; if p < 0 then align = "left" end
				script:write("set label "..align.." '"..p.."' at "..(xo*1.0)..","..math.pow(dy, 1/2)*yo.. "font ',8'\n") -- label
				script:write("set label center '"..(1).."' at "..math.sqrt(dx)*xo..","..(yo*dy).. "font ',8'\n") -- label
			end
		
			local function round(num, quantum)
	  			return math.floor(num / quantum + 0.5) * quantum
			end
		
			for _, dataset in ipairs(plot) do
				local valx = dataset.val[1]
				local valy = dataset.val[2]
				local facx = valx[#valx-1]/valx[#valx]
				local facy = valy[#valy-1]/valy[#valy]
				local rate = math.log(facy) / math.log(facx)
				rate = round(rate, slope.quantum or 1)
				local xo,yo
				local at = slope.at or "last"
				if at == "min" then
					xo, yo = dataset.min[1], dataset.min[2]
				elseif at == "max" then
					xo, yo = dataset.max[1], dataset.max[2]
				elseif at == "last" then
					xo, yo = dataset.val[1][#dataset.val[1]], dataset.val[2][#dataset.val[2]]
				elseif at == "first" then
					xo, yo = dataset.val[1][1], dataset.val[2][1]
				else print("slope.at: only min,max,last,first"); exit(); end
									
				drawSlopTri(xo, yo*1.5, slope.dy, rate)							
			end -- end datasets
		end -- end slope
	
		------------------------------------------------------------------------
		-- Execute Plot
		------------------------------------------------------------------------
	
		-- build up the plot command, layer by layer in non-multiplotmode
		if plotDim == 2 then script:write ("plot  \\\n");
		else 				 script:write ("splot  \\\n"); end

		for s, dataset in ipairs(plot) do
		
			-- get the data column mapping
			local map = table.concat(dataset.map, ":")
				
			-- determine the plot style - data source table has priority
			local style = dataset.style
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
			if not dataset.label then dataset.label = "" end
					
			if s > 1 then script:write(", \\\n");	end
					
			if dataset.file then		
				script:write ("\"", dataset.file, "\"",
								  " using ", map, 
								  " title '", dataset.label, "'", 
								  " with ", style)
			end
			if dataset.func then
				script:write (" ", dataset.func," ",
								  " title '", dataset.label, "'", 
								  " with ", style)
			end 
			if dataset.data then
				tmpDataFile = table.concat({tmpPath.."tmp",p,s,"cols.dat"}, "_")
				gnuplot.write_data(tmpDataFile, dataset.data, dataset.row)
				script:write ("\"", tmpDataFile, "\"",
								  " using ", map, 
								  " title '", dataset.label, "'", 
								  " with ", style)
			end 
		end	-- end datasets		  
																  
		script:write("\n")	

		------------------------------------------------------------------------
		-- cleanup
		------------------------------------------------------------------------

		if slope then
			script:write("unset label\n")	
			script:write("unset arrow\n")	
		end			
		
	end -- end plots
	
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
