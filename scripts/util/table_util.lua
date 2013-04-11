--[[!
-- \defgroup scripts_util_table Table Utility
-- \ingroup scripts_util
-- Utility functions for easy but advanced table access.
-- \{
]]--

table = table or {}

function table.getArraySizes(t)

	local min = 1e1000
	local max = 0
	
	for k, v in pairs(t) do
		if type(k) == "number" and math.floor(k)==k and v ~= nil then
			min = math.min(min, k)
			max = math.max(max, k)
		end
	end
	
	if max == 0 then return 0,0 end
	
	local conseqMax = 0
	for i = min, max do
		if t[i] ~= nil then
			conseqMax = i
		end
	end

	return min, conseqMax
end

function table.print(data, style)

	local title = style.title
	local format = style.format
	local vline = false
	if type(style.vline) == "boolean" then vline = style.vline end
	local hline = false
	if type(style.hline) == "boolean" then hline = style.hline end


	local numCols = #data
	local maxRows = 0
	local minRows = 1e1000
	local maxSize = {}
	for col = 1,numCols do
		local minRow, maxRow = table.getArraySizes(data[col])
		minRows = math.min(minRows, minRow)
		maxRows = math.max(maxRows, maxRow)		
		
		maxSize[col] = 0
		for row = minRow, maxRow do
			local s = data[col][row];
			if format ~= nil and format[col] ~= nil 
			   and s ~= nil and type(s) == "number" then 
				s = string.format(format[col], s)
			end
			maxSize[col] = math.max(maxSize[col], #tostring(s))
		end
	end
	
	if title ~= nil then
		local vsize = 0
		for col = 1, #title do
			maxSize[col] = math.max(maxSize[col], #tostring(title[col]))
			local s = string.format("%"..maxSize[col].."s", tostring(title[col]))
			write(" "..s.." ");
			vsize = vsize + maxSize[col]+2;
			if hline == true then
				vsize = vsize + 1;
				write("|");
			end
		end
		write("\n")
		if vline == true then
			write(string.rep("-", vsize), "\n");
		end
	end
	
	for row = minRows, maxRows do
		for col = 1, numCols do
			local s = data[col][row];
			if format ~= nil and format[col] ~= nil 
				and s ~= nil and type(s) == "number"  then 
				s = string.format(format[col], s)
			end
			
			if s ~= nil then 
				s = string.format("%"..maxSize[col].."s", tostring(s))
			else 
				s = string.format("%"..maxSize[col].."s", " ")
			end
			write (" "..s.." ")
			if hline == true then
				write("|");
			end
		end
		write("\n");
	end		
end

-- end group scripts_util_table
--[[!
\}
]]--
