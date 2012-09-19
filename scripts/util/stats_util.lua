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
			print("file is empty, writing header...")
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

--[[
A={}
A.col1 = 1
A.col2 = "test"
A.col3 = 3

B={}
B.col1 = 5
B.col3 = 6
C = {col1="hey", col4="ho"}
util.printFormattedTable({A, B, C}, {"col1", "col3"}, "c", true)
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
--! @param header	colums to be printed, like {"colA", "colB"}
--! @param type		type of column padding ("r", "l" or "c", default "l")
--! @param bNumbers	if true, print a colum with the number of the table (default false)
function util.printFormattedTable(tab, header, type, bNumbers)
	local length = {}
	local j, col, i, v
	for j, col in ipairs(header) do
		length[col] = string.len(tostring(col))
		for i, v in ipairs(tab) do
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
		for i, v in ipairs(tab) do
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
	for i, v in ipairs(tab) do
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
util.printFormattedTable({A, B, C}, {"col1", "col3"}, "c", true)

  # : |  1  |  2  |   3   | 
----------------------------
col1: |  1  |  5  |  hey  | 
col3: |  3  |  6  |       | 

]]--

--! util.printFormattedTable
--! prints stats to the console, sideways
--! @param tab		table in the form t={{["colA"
--! @param header	colums to be printed, like {"colA", "colB"}
--! @param type		type of row padding ("r", "l" or "c", default "l")
--! @param bNumbers	if true, print a colum with the number of the table (default false)
function util.printFormattedTableSideways(tab, header, type, bNumbers)
	local length = {}
	local j, col, i, v
	local t={}
	local tcol = nil
	if type ~= nil then
		tcol = string.sub(type, 0, 1)
	end
	for i, v in ipairs(tab) do
		length[i] = 0
		
		if type ~=nil and string.len(type) >= i+1 then
			t[i] = string.sub(type, i+1, 1)
		else
			t[i] = type
		end
		
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
		for i, v in ipairs(tab) do
			out=out..util.adjuststring(tostring(i), length[i], t[i]).." | "
		end
		print(out)
		local totallen = headerlength+4
		for i, v in ipairs(tab) do
			totallen = totallen+length[i]+3
		end
		print(util.fill(totallen, "-"))
	end
	
	for j, col in ipairs(header) do
		out = util.adjuststring(tostring(col), headerlength, tcol)..": | "		
		for i, v in ipairs(tab) do
			if v[col] ~= nil then
				out=out..util.adjuststring(tostring(v[col]), length[i], t[i])				
			else
				out=out..util.fill(length[i])
			end
			out=out.." | "	
		end	
		print(out)
		
	end
end

