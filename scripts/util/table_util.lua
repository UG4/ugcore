table = table or {}

function table.print(data, style)

	local title = style.title
	local format = style.format
	local vline = false
	if type(style.vline) == "boolean" then vline = style.vline end
	local hline = false
	if type(style.hline) == "boolean" then hline = style.hline end


	local numCols = #data
	local numRows = 0
	local maxSize = {}
	for col = 1,numCols do
		numRows = math.max(numRows, #data[col])
		
		maxSize[col] = 0
		for row = 1, #data[col] do
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
	
	for row = 1, numRows do
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

