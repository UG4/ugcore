--[[!
-- \defgroup scripts_util_profiler Lua Utility Scripts
-- \ingroup scripts_util
-- \brief Helpers for profiling output
-- \{
]]--

util = util or {}

--[[
Prints the total time spent in the profile-node with the given name on the local proc.
nodeName may contain spaces to ease formatting.
]]--
function util.PrintProfile_TotalTime(nodeName)
	if(GetProfilerAvailable() == true) then
		rawName = string.gsub(nodeName, "^%s*(.-)%s*$", "%1")
		pn = GetProfileNode(rawName)
		if(pn:is_valid() == true) then
			print(nodeName .. " " .. pn:get_avg_total_time_ms() / 1000 .. " s")
			return
		end
	end
	print(nodeName .. " ---")
end

--[[
Prints the maximum total time spent in the profile-node with the given name on any proc.
nodeName may contain spaces to ease formatting.
]]--
function util.PrintProfile_MaxTotalTime(nodeName)
	if(GetProfilerAvailable() == true) then
		local rawName = string.gsub(nodeName, "^%s*(.-)%s*$", "%1")
		local pn = GetProfileNode(rawName)
		local t = 0
		if(pn:is_valid() == true) then
			t = pn:get_avg_total_time_ms() / 1000
		end
		t = ParallelMax(t)
		print(nodeName .. " " .. t .. " s")
	else
		print(nodeName .. " ---")
	end
end

--[[
Prints the maximum total time spent in the given profile nodes. profName has to
be a table with profilenames starting at index 1.
One global allreduce is performed to obtain the max-times on all procs.
]]--
function util.PrintProfiles_MaxTotalTimes(profNames)
	local maxLen = 0
	local profTimes = {}

	for i, name in ipairs(profNames) do
		local len = string.len(name)
		if len > maxLen then maxLen = len; end

		rawName = string.gsub(name, "^%s*(.-)%s*$", "%1")
		pn = GetProfileNode(rawName)
		if pn:is_valid() == true then
			profTimes[i] = pn:get_avg_total_time_ms() / 1000
		else
			profTimes[i] = 0
		end
	end

	maxLen = maxLen + 3
	profTimes = ParallelVecMax(profTimes)

	for i, time in ipairs(profTimes) do
		print(profNames[i] .. ":" .. string.rep("-", maxLen - string.len(profNames[i])) .. profTimes[i])
	end
end

--[[!
\}
]]--
