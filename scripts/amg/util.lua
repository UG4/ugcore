
printf = function(s,...)
	print(s:format(...))
end -- function

formatf = function(s, ...)
	return s:format(...)
end

function fsize (file)
	local current = file:seek()      -- get current position
    local size = file:seek("end")    -- get file size
    file:seek("set", current)        -- restore position
    return size
end

function getStats(stats, bHeader, seperator, bStats, seperator2)
	local output=""
	for i,v in ipairs(stats) do
		if bHeader then output = output..v[1]..seperator end
		if bStats  then output = output..v[2]..seperator2 end
	end
	output = output.."\n"
	return output
end

function writeFileStats(stats, filename)
	local output = io.open(filename, "a")
	
	if fsize(output) == 0 then
		output:write(getStats(stats, true, "\t", false, "\t"))		
	end
	output:write(getStats(stats, false, "\t", true, "\t"))
end

function printStats(stats)
	print(getStats(stats, true, ": ", true, "\n"))
end

function PrintParallelProfileNode(name)
	pn = GetProfileNode(name)
	t = pn:get_avg_total_time_ms()/to100 * 100
	tmin = ParallelMin(t)
	tmax = ParallelMax(t)
	printf("%s:\n%.2f %%, min: %.2f %%, max: %.2f %%", name, t, tmin, tmax)
end