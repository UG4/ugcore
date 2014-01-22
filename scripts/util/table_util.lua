--[[!
-- \defgroup scripts_util_table Table Utility
-- \ingroup scripts_util
-- Utility functions for easy but advanced table access.
-- \{
]]--

table = table or {}

--! appends a table to a table
function table.append(t1, t2)
	if type(t1) ~= "table" then
		print("table.append called on non-table"); exit();
	end
	
	if type(t2) == "table" then
       for _,v in ipairs(t2) do
            t1[#t1+1] = v
       end
    else
       t1[#t1+1] = t2
    end	
end

--! returns the smallest integer key of the table (even negative or null if present)
function table.imin(t)
	local min = math.huge
	for k, _ in pairs(t) do
		if type(k) == "number" and math.floor(k)==k then
			min = math.min(min, k)
		end
	end
	return min
end

--! returns the largest integer key of the table (even in non-consecutive arrays)
function table.imax(t)
	local max = -math.huge
	for k, _ in pairs(t) do
		if type(k) == "number" and math.floor(k)==k then
			max = math.max(max, k)
		end
	end
	return max
end

--! prints the table
function table.print(data, style)

	-- get format
	local heading = style.heading
	local format = style.format or {}
	local vline = style.vline or false
	local hline = style.hline or false
	local forNil = style.forNil or " "

	-- compute row range	
	local minRow = math.huge
	local maxRow = -math.huge
	for _, column in ipairs(data) do
		if type(column) == "table" then
			minRow = math.min(minRow, table.imin(column))
			maxRow = math.max(maxRow, table.imax(column))		
		else 
			print("table.print: expect consecutive lua-table-entries as column")
		end
	end		

	-- compute column width (in data)
	local width = {}
	for col, _ in ipairs(data) do
		width[col] = 0
		for row = minRow, maxRow do
			local s = data[col][row];
			if s ~= nil then
				if type(s) == "number" and format[col] ~= nil then 
					s = string.format(format[col], s)
				end
				width[col] = math.max(width[col], #tostring(s))
			end
		end
	end
	
	-- compute column width (in header)
	-- write header and vline
	if heading ~= nil then
		local linesize = 0
		for col = 1, #heading do
			width[col] = math.max(width[col], #tostring(heading[col]))
			local s = string.format("%"..width[col].."s", tostring(heading[col]))
			write(" "..s.." ");
			linesize = linesize + width[col] + 2;
			if hline == true then
				linesize = linesize + 1;
				write("|");
			end
		end
		write("\n")
		if vline == true then
			write(string.rep("-", linesize), "\n");
		end
	end
	
	-- write data
	for row = minRow, maxRow do
		for col, _ in ipairs(data) do
			local s = data[col][row];
			
			if s ~= nil then 
				if type(s) == "number" and format[col] ~= nil then
					s = string.format(format[col], s)
				else
					s = tostring(s)
				end
			else		
				s = forNil
			end

			write (" ".. string.format("%"..width[col].."s", s) .." ")
			
			if hline == true then write("|"); end
		end
		write("\n");
	end		
end

-- end group scripts_util_table
--[[!
\}
]]--
