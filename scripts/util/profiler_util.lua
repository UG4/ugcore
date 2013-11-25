--[[!
-- \defgroup scripts_util_profiler Lua Utility Scripts
-- \ingroup scripts_util
-- \brief Helpers for profiling output
-- \{
]]--

util = util or {}

--[[
Prints the total time spent in the profile-node with the given name.
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

--[[!
\}
]]--
