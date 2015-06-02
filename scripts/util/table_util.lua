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

--! returns shallow-copy of table (i.e., only top level values)
function table.shallowcopy(t)
    local copy
    if type(t) == "table" then
        copy = {}
        for k, v in pairs(t) do
            copy[k] = v
        end
    else -- number, string, boolean, etc
        copy = t
    end
    return copy
end

--! returns shallow-copy of integer-key part of table (i.e., only top level values)
function table.ishallowcopy(t)
    local copy
    if type(t) == "table" then
        copy = {}
        for k, v in ipairs(t) do
            copy[k] = v
        end
    else -- number, string, boolean, etc
        copy = t
    end
    return copy
end


--! returns deep-copy of table (i.e., copies recursive all contained tables as well)
function table.deepcopy(t)
    local copy
    if type(t) == "table" then
        copy = {}
        for k, v in next, t, nil do
            copy[table.deepcopy(k)] = table.deepcopy(v)
        end
        setmetatable(copy, table.deepcopy(getmetatable(t)))
    else -- number, string, boolean, etc
        copy = t
    end
    return copy
end

--! returns deep-copy of integer-key part of table (i.e., copies recursive all contained tables as well)
function table.ideepcopy(t)
    local copy
    if type(t) == "table" then
        copy = {}
        for k, v in ipairs(t) do
            copy[k] = table.deepcopy(v)
        end
    else -- number, string, boolean, etc
        copy = t
    end
    return copy
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

--! returns iterator function iterating through first consecutive integer-key range (even starting at negative or null)
function iipairs(tab)
	local n = table.imin(tab)-1
	local t = tab
	return function (i)
      	 n = n + 1
      	 local v = t[n]
     	 if v then
        	return n, v
      	else
      		return nil
    	end
    end
end


--! checks if a value is contained in a table
function table.contains(t, value)
	for _,v in pairs(t) do
		if v == value then return true end
	end
	return false
end

--! prints the table
function table.print(data, style)

	-- get format
	if style == nil then style = {} end
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
		for col, _ in ipairs(data) do
			local s = 	heading[col]
			if s ~= nil then s = tostring(s)
			else			 s = tostring(forNil) end
		
			width[col] = math.max(width[col], #s)
			write(" ".. string.format("%"..width[col].."s", s) .." ");
			
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


--! traverses a table recursively and returns the length of the longest identifyer
function table.len_of_longest_identifyer(t)
	local maxLen = 0
	for n, v in pairs(t) do
		if type(v) == "table" then
			maxLen = math.max(maxLen, table.len_of_longest_identifyer(v))
		else
			maxLen = math.max(maxLen, string.len(n))
		end
	end
	return maxLen
end


function table.print_flat(t, prefix, maxLen)
	if prefix == nil then prefix = "" end
	if maxLen == nil then
		maxLen = table.len_of_longest_identifyer(t)
	end

	local namesSorted = {}
	local c = 1
	for n, _ in pairs(t) do
		namesSorted[c] = n
		c = c + 1
	end
	table.sort(namesSorted)

	for i = 1, c do
		local n = namesSorted[i]
		if n ~= nil then
			local v = t[n]
			
			local l = string.len(n)
			local pn = n
			if l > maxLen then
				pn = string.sub(pn, 1, maxLen)
				l = maxLen
			end
			pn = pn .. " " .. string.rep(".", 2 + maxLen - l)

			if type(v) == "table" then
				print(prefix, pn, " = {")
				table.print_flat(t, prefix .. "  ", maxLen)
				print(prefix, "}")
			else
				print(prefix, pn, " = ", v)
			end
		end
	end
end


--! Recursively adds entries from table t2 to table t1, if those entries were
--! not already contained in t1.
--! @returns t1
function table.merge_inplace(t1, t2)
	for n, v in pairs(t2) do
		if type(v) == "table" then
			if t1[n] == nil then
				t1[n] = table.deepcopy(v)
			elseif type(t1[n]) == "table" then
				table.merge_inplace(t1[n], v)
			end
		elseif t1[n] == nil then
			t1[n] = v
		end
	end
	return t1
end

--! Creates a new table which contains the entries from tables t1 and t2.
--! key-value pairs from table t1 have higher priority than those from t2,
--! meaning that if a key is found in both tables, the key-value pair from
--! t1 is considered only.
--! @returns A new table resulting from merging t1 and t2
function table.merge(t1, t2)
	local t = table.deepcopy(t1)
	table.merge_inplace(t, t2)
	return t
end

-- end group scripts_util_table
--[[!
\}
]]--
