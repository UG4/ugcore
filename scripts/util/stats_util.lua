--[[!
\addtogroup scripts_util_stats
\{
]]--
-- created by Martin Rupp, 06.02.2012, martin.rupp@gcsc.uni-frankfurt.de
-- for documentation: stats_util.doxygen
util = util or {}

--! util.getStats
--! internal method
function util.getStats(stats, bHeader, seperator, bStats, seperator2)
	local output=""
	for i,v in ipairs(stats) do
		if bHeader then output = output..v[1]..seperator end
		if bStats then
			if v[2] == nil then
				print("value for item "..v[1].." is nil!")
				output = output.." "..seperator2
			else				
				output = output..tostring(v[2])..seperator2
			end
		end
	end
	output = output.."\n"
	return output
end


--! util.writeFileStats
--! writes stats to a file
--! @param stats
--! @param filename
--! @param seperator (default " \t")
function util.writeFileStats(stats, filename, seperator)
	if seperator == nil then
		seperator = " \t"
	end
	local output = io.open(filename, "a")
	
	if output == null then 
	   print("Could not open"..filename)
	else
		if fsize(output) == 0 then
			-- print("file is empty, writing header...")
			output:write("#"..util.getStats(stats, true, seperator, false, seperator))		
		end
		output:write(util.getStats(stats, false, seperator, true, seperator))
	end
end


--! util.fill
--! returns a string consisting of N times the character c
--! @param N number of times c is to be repeated
--! @param c character to repeat (if omitted, " ")
function util.fill(N, c)
	local s=""
	if c == nil then c = " " end
	for i=1,N do
		s=s..c
	end
	return s
end

--! util.adjuststring
--! returns a string with whitespace left and right so that total string length
--! is len
--! @param str   string used
--! @param len   total length of resulting string
--! @param type  padding type: "l" = string is on left, "c" = centered, "r" = right. default "l"
function util.adjuststring(str, len, type)
	if type == nil then
		type = "l"
	end
	local l=string.len(str)
	local m=0
	if type == "c" then
		m=math.ceil((len-l)/2)
	elseif type == "r" then
		m=len-l
	end
	return util.fill(m)..str..util.fill(len-m-l)
end

--! util.printStats
--! prints stats to the console
--! @param stats
function util.printStats(stats)
	local maxlen=2
	for i,v in ipairs(stats) do
		maxlen = math.max(maxlen, string.len(v[1]))		
	end

	for i,v in ipairs(stats) do
		print(util.adjuststring(v[1], maxlen)..": "..tostring(v[2]))		
	end	
end


function util.statsUtilGetHeader(header, tab)
	if header == nil then
		local h = {}
		header = {}		
		for i, v in pairs(tab) do
			for j, col in pairs(v) do
				if h[j] == nil then
					h[j] = j		
					table.insert(header, j)
				end		 
			end
		end
	end
	return header
end

--[[
A={}
A.col1 = 1
A.col2 = "test"
A.col3 = 3

B={}
B.col1 = 5
B.col3 = 6
C = {col1="hey", col4="ho"}
util.printFormattedTable({A, B, C}, true, "c", {"col1", "col3"})
output:
# | col1 | col3 |
-----------------
1 |   1  |   3  |
2 |   5  |   6  |
3 |  hey |      |
]]--

--! util.printFormattedTable
--! prints stats to the console
--! @param tab		table in the form t={{["colA"
--! @param bNumbers	if true, print a colum with the key of the table (default false)
--! @param type		type of column padding ("r", "l" or "c", default "l")
--! @param header	colums to be printed, like {"colA", "colB"}. if nil, all.

function util.printFormattedTable(tab, bNumbers, type, header)
	local length = {}
	local j, col, i, v

	header = util.statsUtilGetHeader(header, tab)
	
	for j, col in ipairs(header) do
		length[col] = string.len(tostring(col))
		for i, v in pairs(tab) do
			if v[col] ~= nil then
				length[col] = math.max(string.len(tostring(v[col])), length[col])				
			end
		end	
		length[col] = length[col]+2	
	end	
	
	local t={}
	out=""
	
	local numberslen = 0
	if bNumbers ~=nil and bNumbers == true then
		for i, v in pairs(tab) do
			numberslen = math.max(numberslen, string.len(tostring(i)))
		end
		numberslen=numberslen+1
		out=out..util.adjuststring("#", numberslen, "l")
		
	end
	out=out.."|"
	
	for j, col in ipairs(header) do
		if type ~=nil and string.len(type) >= j then
			t[j] = string.sub(type, j, 1)
		else
			t[j] = type
		end
		out=out..util.adjuststring(tostring(col), length[col], t[j]).."|"		
	end
	print(out)
	out=util.fill(numberslen+1, "-")
	for j, col in ipairs(header) do
		out=out..util.fill(length[col]+1, "-")				
	end
	print(out)
	for i, v in pairs(tab) do
		out=""
		if bNumbers ~=nil and bNumbers == true then
			out=out..util.adjuststring(tostring(i), numberslen, "l")
		end
		out=out.."|"
		for j, col in ipairs(header) do
			if v[col] ~= nil then
				out=out..util.adjuststring(tostring(v[col]), length[col], t[j])				
			else
				out=out..util.fill(length[col])
			end
			out=out.."|"	
		end	
		print(out)
		
	end
end
--[[
another example:
jacResult = {}
jacResult.avgConv = 0.2
jacResult.last = 0.4

sgsResult = {}
sgsResult.avgConv = 0.1
sgsResult.last = 0.3
tab = {}
tab["jac"] = jacResult
tab.sgs = sgsResult

util.printFormattedTable(tab, true, "c")
]]-- 

--! util.printFormattedTableSideways
--! prints stats to the console, sideways
--! @param tab		table in the form t={{["colA"
--! @param type		type of row padding ("r", "l" or "c", default "l")
--! @param bNumbers	if true, print a colum with the number of the table (default false)
--! @param header	colums to be printed, like {"colA", "colB"}, if nil: all
function util.printFormattedTableSideways(tab, bNumbers,  header)
	local length = {}
	local j, col, i, v
	local t={}
	local tcol = nil
	
	header = util.statsUtilGetHeader(header, tab)

	for i, v in pairs(tab) do
		length[i] = 0
				
		for j, col in ipairs(header) do
			if v[col] ~= nil then
				length[i] = math.max(string.len(tostring(v[col])), length[i])				
			end
		end	
		length[i] = length[i]+2	
	end	
	
	local headerlength = 0
	for j, col in ipairs(header) do
		headerlength = math.max(string.len(tostring(col)), headerlength)
	end
	
	if bNumbers ~= nil and bNumbers == true then
		out = util.adjuststring("#", headerlength, tcol)..": | "		
		for i, v in pairs(tab) do
			out=out..util.adjuststring(tostring(i), length[i], "c").." | "
		end
		print(out)
		local totallen = headerlength+4
		for i, v in pairs(tab) do
			totallen = totallen+length[i]+3
		end
		print(util.fill(totallen, "-"))
	end
	
	for j, col in ipairs(header) do
		out = util.adjuststring(tostring(col), headerlength, tcol)..": | "		
		for i, v in pairs(tab) do
			if v[col] ~= nil then
				out=out..util.adjuststring(tostring(v[col]), length[i], "c")				
			else
				out=out..util.fill(length[i])
			end
			out=out.." | "	
		end	
		print(out)
		
	end
end


function util.StringTableFromTable(tab, header, s)
	ug_assert(tab ~= nil)
	local length = {}
	local j, col, i, v
	if s == nil then
		s = StringTable()
	end

	header = util.statsUtilGetHeader(header, tab)	
	local x=1
	local y=0
	for j, col in ipairs(header) do
		s:set(0, x, col)
		x = x+1		
	end
	
	y =1
	for j, col in pairs(tab) do
		s:set(y, 0, j)
		y = y +1
	end
	
	y=1
	for i, v in pairs(tab) do
		x=1
		for j, col in ipairs(header) do
			if v[col] ~= nil then
				s:set(y, x, tostring(v[col]))				
			end				
			x=x+1
		end
		y=y+1
	end
	return s
end



--! inserts the table tab into resTable in a "flat" way
--! that means, subtables are getting indices so that
--! the resulting table has no subtables
--! example:
--! tab = {a=2, b="hi", c={word=3, nose=4} }
--! res = {}
--! flattenTable(tab, "tab", res)
--! res = {
--!  {"tab.a",  2},
--!  {"tab.b", "hi"},
--!  {"tab.c.word", 3},
--!  {"tab.c.nose", 4}
--! }
--! especially usefull for util.writeFileStats
--! note: syntax {name, value} used instead of name=value
--! because this way the order is stable
function util.flattenTable(tab, name, resTable)
	local tName = ""
	if name ~= nil then
		tName = name.."."
	end
	if type(tab) == "table" then
		for i, v in pairs(tab) do
			util.flattenTable(v, tName..i, resTable)
		end
	elseif type(tab) ~= "userdata" then		
		table.insert(resTable, {name, tab})
	end
	return resTable
end

-- end group scripts_util_stats
--[[!
\}
]]--
