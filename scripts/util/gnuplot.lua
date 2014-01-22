--[[!
\file gnuplot.lua
\defgroup scripts_util_gnuplot GNUPlot Utility
create a globally seen package (all non-local functions call then be called using gnuplot.xxx.)
\ingroup scripts_util
\{
]]--

gnuplot = gnuplot or {}

function gnuplot.getArraySizes(t)

	local min = 1e1000
	local max = 0
	
	for k, v in pairs(t) do
		if type(k) == "number" and math.floor(k)==k and v ~= nil then
			min = math.min(min, k)
			max = math.max(max, k)
		end
	end
	
	if max == 0 and min == 1e1000 then return 0,-1 end
	
	local conseqMax = 0
	for i = min, max do
		if t[i] ~= nil then
			conseqMax = i
		end
	end

	return min, conseqMax
end

-- writes the passed data array to file
--
-- filename			the output file
-- data				table containing the data
-- passRows		   flag indication if data array stores rows (default = false)
function gnuplot.write_data(filename, data, passRows)
	-- open file
	local file = io.open(filename, "w+")
	if (not file) then
		io.stderr:write("Gnuplot Error: cannot open output file: '")
		io.stderr:write(filename .. " '\n");
		return 1
	end
	
	if type(data) ~= "table" then
		io.stderr:write("Gnuplot Error: Expect table as data.\n");
		return 1	
	end
	if data[1]==nil then
		return 0;
	end

	-- default is table data by row
	if passRows == nil then
		passRows = false
	end
	-- pure data is rowwise
	if type(data[1]) == "number" then
		passRows = true
	end

	local plainArray, rowSize
	local minRow = nil
	local maxRow = nil 
	local bHasMatrixData = false
	
	if not passRows then
		for i, item in ipairs(data) do
			if not(type(item) == "table") then
				io.stderr:write("Gnuplot Error: Data array must contain only tables.\n");
				return 1						
			end
		end

		for i = 1, #data do
			if type(data[i][1]) == "table" then
				bHasMatrixData = true		
			end
		end
		
		-- check column sizes
		if bHasMatrixData == false then
			for i, item in ipairs(data) do
				local min, max = gnuplot.getArraySizes(item)
				if minRow == nil then minRow = min
				elseif minRow ~= min then
					io.stderr:write("Gnuplot Error: Data array of mixed sized columns.");
					io.stderr:write("Expected min row: "..minRow..", got: "..min.."\n");
					return 1						
				end
				if maxRow == nil then maxRow = max
				elseif maxRow ~= max then
					io.stderr:write("Gnuplot Error: Data array of mixed sized columns.");
					io.stderr:write("Expected max row: "..maxRow..", got: "..max.."\n");
					return 1						
				end
			end
		end
		
		-- write column
		if bHasMatrixData == false then
			for i = minRow,maxRow do
				for j = 1,#data do
					file:write(data[j][i], " ")
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
	else
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

local tmpPath = "./tmp-gnuplot/" -- some tmp path
local plotImagePath = "./" -- output of data
local tmpScriptName = "tmp_gnuplot_script_" -- name for temporary script

-- filename		output filename
-- datasource	table of data to be plotted
-- options		table of options
function gnuplot.plot(filename, datasource, options)

	if (not datasource or type(datasource) ~= "table") then
		io.stderr:write("Gnuplot Error: a data source dictionary must be provided!\n");
		exit();
	end

	-- check for 2d or 3d data
	local is3d = nil
	for i, source in ipairs(datasource) do

		if ( #source ~= 0 ) then -- use default 
			if ( #source == 2 ) then
				if options.is3d ~= nil and options.is3d == true then
					io.stderr:write("Gnuplot Error: Explict 3d data request, but 2d data given.\n"); exit();
				end
				if is3d ~= nil and is3d == true then 
					io.stderr:write("Gnuplot Error: Mixed 2d/3d data.\n"); exit();
				else is3d = false 
				end
			elseif ( #source == 3 ) then  
				if options.is3d ~= nil and options.is3d == false then
					io.stderr:write("Gnuplot Error: Explict 2d data request, but 3d data given.\n"); exit();
				end
				if is3d ~= nil and is3d == false then 
					io.stderr:write("Gnuplot Error: Mixed 2d/3d data.\n"); exit();
				else is3d = true 
				end			
			else io.stderr:write("Gnuplot Error: pass 0, 2 or 3 columns as data selection.\n");
				 exit();
			end		
		end
	end
	if options.is3d == nil then 
		if is3d ~= nil then 
			options.is3d = is3d
		else 
			options.is3d = false 
		end
	end
	
	-- select sensible defaults
	if ( type(options) ~= "table" ) then
		options = { title = "", xlabel = "", ylabel = "", zlabel = ""}
	end
	if options.title == nil then options.title = "" end
	if options.xlabel == nil then options.xlabel = "" end
	if options.ylabel == nil then options.ylabel = "" end
	if options.zlabel == nil then options.zlabel = "" end

	-- check output file name
	if filename ~= nil then
		if ( type(filename) ~= "string") then 
			io.stderr:write("Gnuplot Error: Filename must be a string.");
			exit();
		else
			filename = string.gsub(filename, " ", "_" )
		end
	end
	
	-- if sensible ranges are specified use them otherwise autoscale
	if options.range ~= nil then
		if type(options.range) == "table" then
			if options.range.xrange == nil then 
				if type(options.range[1]) == "number" and
				   type(options.range[2]) == "number" then
					options.range.xrange = {options.range[1], options.range[2]}
				elseif type(options.range[1]) == "table" then
					options.range.xrange = options.range[1] 
				end
			end
			if options.range.yrange == nil then 
				if type(options.range[3]) == "number" and
				   type(options.range[4]) == "number" then
					options.range.yrange = {options.range[3], options.range[4]}
				elseif type(options.range[2]) == "table" then
					options.range.yrange = options.range[2] 
				end
			end
			if options.range.zrange == nil then 
				if type(options.range[5]) == "number" and
				   type(options.range[6]) == "number" then
					options.range.zrange = {options.range[5], options.range[6]}
				elseif type(options.range[3]) == "table" then
					options.range.zrange = options.range[3] 
				end
			end
		else
			io.stderr:write("Gnuplot Error: range must be a table.");
			exit();
		end
	else
		options.range = {};
	end

	-- open a temporary script file
	plotScriptName = tmpPath..tmpScriptName..string.gsub(filename, "[./]", "_")..".gnu"
	local plotScript = io.open(plotScriptName, "w+")
	if (not plotScript) then
		os.execute("mkdir " .. tmpPath)
		plotScript = io.open(plotScriptName, "w+")
		if (not plotScript) then		
			io.stderr:write("Gnuplot Error: cannot open output file: '")
			io.stderr:write(plotScriptName .. " '\n");
			return 2
		end
	end

	-- check for multiplot
	local multiplot = false; 
	local multiplotjoined = false;
	local MultiPlotRows, MultiPlotCols
	if type(options.multiplot) == "boolean" then multiplot = options.multiplot end
	if ( multiplot == true ) then

		if type(options.multiplotrows) == "number" then 
			MultiPlotRows = options.multiplotrows 
		else
			MultiPlotRows = math.ceil(math.sqrt(#datasource))
		end
		MultiPlotCols = math.ceil(#datasource / MultiPlotRows )
		
		if type(options.multiplotjoined) == "boolean" then multiplotjoined = options.multiplotjoined end	
	end

	-- specify the output file
	plotScript:write("reset\n")

	-- find font and fontsize
	local font = "Verdana"
	local fontsize = 10
	if options.font then font = options.font end
	if options.fontsize then fontsize = options.fontsize end
	if ( multiplot == true ) then
		fontsize = math.ceil(fontsize / math.max(MultiPlotRows, MultiPlotCols))
	end
	
	-- set terminal currently only pdf
	if filename == nil then
		plotScript:write("set terminal x11 persist raise\n\n")
	elseif string.find(filename, ".pdf", -4) ~= nil then
		if options.font == nil then font = "Times-Roman" end
		plotScript:write("set term pdfcairo enhanced font '"..font..","..fontsize.."'\n\n")
		plotScript:write("set output \"", plotImagePath, filename,"\"\n\n")
	elseif string.find(filename, ".eps", -4) ~= nil then
		plotScript:write("set term postscript eps enhanced color font '"..font..","..fontsize.."'\n\n")
		plotScript:write("set output \"", plotImagePath, filename,"\"\n\n")
	elseif string.find(filename, ".tex", -4) ~= nil then
		plotScript:write("set term epslatex color colortext font '"..font..","..fontsize.."'\n\n")
		plotScript:write("set output \"", plotImagePath, filename,"\"\n\n")
	elseif string.find(filename, ".svg", -4) ~= nil then
		plotScript:write("set term svg enhanced fname '"..font.."' fsize "..fontsize.."\n\n")
		plotScript:write("set output \"", plotImagePath, filename,"\"\n\n")
	else
		io.stderr:write("Gnuplot Error: wrong file type: '"..filename.."'\n")
		io.stderr:write("Supported endings: pdf, eps, svg.\n")
		return 2		
	end
	
	-- title and axis label
	plotScript:write("set title '"..options.title.."' textcolor lt 2\n")
	plotScript:write("set xlabel '"..options.xlabel.."'\n")
	plotScript:write("set ylabel '"..options.ylabel.."'\n\n")
	if options.is3d == true then
	plotScript:write("set zlabel '"..options.zlabel.."'\n\n")
	end
	
	-- write timestamp
	if options.timestamp then
		if type(options.timestamp) == "string" then
			plotScript:write("set timestamp "..options.timestamp.."\n")
		elseif type(options.timestamp) == "boolean" then
			if options.timestamp == true then
				plotScript:write("set timestamp \"%H:%M:%S  %Y/%m/%d\" bottom\n")
			else 
				plotScript:write("unset timestamp\n"); 
			end
		else
			io.stderr:write("Gnuplot Error: timestamp option not recognized.\n");
			exit();
		end
	else 
		plotScript:write("unset timestamp\n"); 
	end

	-- enable grid
	if options.grid == true then
		if options.is3d == false then plotScript:write("set grid xtics ytics back\n"); 
		else plotScript:write("set grid xtics ytics ztics back\n"); 
		end
	end

	-- locscale
	if options.logscale == true or options.xlogscale == true then
		plotScript:write("set logscale x\n");
	else
		plotScript:write("unset logscale x\n");
	end
	if options.logscale == true or options.ylogscale == true then
		plotScript:write("set logscale y\n"); 
	else
		plotScript:write("unset logscale y\n");
	end
	if options.is3d == true then
	if options.logscale == true or options.zlogscale == true then
		plotScript:write("set logscale z\n"); 
	else
		plotScript:write("unset logscale z\n");
	end
	end
	
	-- write additional options, passed as strings
	for i, source in ipairs(options) do
		if source and type(source) == "string" then
			plotScript:write(source.."\n")
		end
	end	

	-- multiplot sizes
	if ( multiplot == true ) then
		plotScript:write("set multiplot layout ", MultiPlotRows,", ", MultiPlotCols)
		plotScript:write(" title '"..options.title.."' \n\n" )	
	end

	-- remove key (legend)
	if options.key ~= nil and options.key == false then plotScript:write ("unset key\n") end

	-- range
	plotScript:write ("set autoscale\n")
	if options.range.xrange ~= nil then
		plotScript:write ("set xrange [",options.range.xrange[1],":",options.range.xrange[2],"]\n")
	end
	if options.range.yrange ~= nil then
		plotScript:write ("set yrange ["..options.range.yrange[1]..":"..options.range.yrange[2].."]\n")
	end
	if options.range.zrange ~= nil then
		plotScript:write ("set zrange ["..options.range.zrange[1]..":"..options.range.zrange[2].."]\n")
	end
	

	-- loop data sets
	for i, source in ipairs(datasource) do
	
		-- get the data column mapping
		local map = "";
		if options.is3d == false and #source == 2 then
			map = source[1]..":"..source[2]
		elseif options.is3d == true and #source == 3 then
			map = source[1]..":"..source[2]..":"..source[3]
		else
			if options.is3d == false then map = "1:2"
			else map = "1:2:3"
			end
		end

		-- determine the plot type - data source table has priority
		if options.is3d == true then
			plotScript:write("set surface\n") 
			plotScript:write("unset contour\n") 
			local style = "points"
		else
			local style = "linespoints"
		end
		if source.style then
			style = source.style
		end
		
		-- check style		
		if 	style ~= "lines" and
			style ~= "points" and
			style ~= "linespoints" and
			style ~= "boxes" and
			style ~= "dots" and
			style ~= "vectors" and
			style ~= "yerrorbars"
		then
			io.stderr:write("Gnuplot Error: style=\""..style.."\" not supported.\n");
			exit()
		end
		
		-- assign default data label if needed
		if ( not source.label ) then source["label"] = "plot" .. i end

		-- Position plots in multiplot
		if ( multiplot == true ) then
			plotScript:write("unset title \n" )
			
			if multiplotjoined then
				-- left bnd
				if i % MultiPlotCols == 1 then 
					plotScript:write("unset lmargin\n") 
					plotScript:write("set ylabel '"..options.ylabel.."'\n\n")
					plotScript:write("set format y\n")
				else 
					plotScript:write("set lmargin 0\n")
					plotScript:write("unset ylabel\n")
					plotScript:write("set format y \"\"\n")
				end

				-- right bnd
				if i % MultiPlotCols == 0 then 
					plotScript:write("unset rmargin\n") 
				else 
					plotScript:write("set rmargin 0\n")
				end
				
				-- top bnd
				if math.ceil(i / MultiPlotRows) == 1 then 
					plotScript:write("unset tmargin\n") 
				else 
					plotScript:write("set tmargin 0\n")
				end

				-- bottom bnd
				if math.ceil(i / MultiPlotRows) == MultiPlotRows then 
					plotScript:write("unset bmargin\n") 
					plotScript:write("set xlabel '"..options.xlabel.."'\n\n")
					plotScript:write("set format x\n")
				else 
					plotScript:write("set bmargin 0\n")
					plotScript:write("unset xlabel\n")
					plotScript:write("set format x \"\"\n")
				end
			end
		end

		-- build up the plot command, layer by layer in non-multiplotmode
		if (i == 1 or multiplot == true) then 	
			if options.is3d == false then 
				 plotScript:write ("plot  \\\n");
			else plotScript:write ("splot  \\\n"); end
		else 				
			plotScript:write(", \\\n");	
		end
		
		if source.file ~= nil then		
			plotScript:write ("\"", source.file, "\"",
							  " using ", map, 
							  " title '", source.label, "'", 
							  " with ", style)
		end
		if source.func ~= nil then
			plotScript:write (" ", source.func," ",
							  " title '", source.label, "'", 
							  " with ", style)
		end 
		if source.data ~= nil then
			tmpDataFile = tmpPath.."tmp_"..source.label.."_"..i..".dat"
			tmpDataFile = string.gsub(tmpDataFile, " ", "_" )
			write_data(tmpDataFile, source.data, source.row)
			plotScript:write ("\"", tmpDataFile, "\"",
							  " using ", map, 
							  " title '", source.label, "'", 
							  " with ", style)
		end 
								  
		if (multiplot == true) then
			plotScript:write("\n")
		end						  
	end
	plotScript:write ("\nunset multiplot\n");
	plotScript:close()

	-- launch gnuplot and run in background
	os.execute("gnuplot "..plotScriptName.." --persist &")
	return 0
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
--	{ label = "catsdogs",  file = "./TestData/dogdata", style = "points", 1, 5 },
--	{ label = "catsdogs3", file = "./TestData/dogdata", style = "points", 1, 6 },
--	{ label = "catsdogs4", file = "./TestData/dogdata", style = "lines", 1, 5 },
--	{ label = "catsdogs5", file = "./TestData/dogdata", style = "lines", 1, 6 },
--	{ label = "catsdogs5", file = "./TestData/dogdata", style = "lines", 1, 6 },
--	{ label = "sinuskurve",func = "sin(x)", style = "lines"},
	{ label = "rowdata_13", data = rowdata, row = true, style = "lines", 1, 3},
	{ label = "rowdata_12", data = rowdata, row = true, style = "lines"},
	{ label = "col data 12", data = coldata, style = "lines"},
	{ label = "col data 13", data = coldata, style = "lines", 1, 3}
}

options = {	title =			"Title", 
			xlabel =		"xAxis",
			ylabel =		"yAxis",
--			range = 		{ xrange = {-100, 600}, yrange = {-5, 5} },
--			range = 		{ {-100, 600}, {-5, 5} },
--			font = 			"Arial",
--			fontsize =		10 / 1,
--			grid = 			true,
--			multiplot = 	true,
--			multiplotrows = 5,
--			multiplotjoined=true,
--			key = 			true,
--			logscale = true,
--			"set some user defined lines to be added to script"
}

--gnuplot.plot("vibration.eps", datasource, options)
--gnuplot.plot("vibration.svg", datasource, options)
--gnuplot.plot("vibration.pdf", datasource, options)
--gnuplot.plot("vibration.tex", datasource, options)
--gnuplot.plot(nil, datasource, options)

]]--

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
