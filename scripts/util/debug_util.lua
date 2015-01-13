
function ug_error(msg, otherMsg, backtraceSkipLevel)
	print("\n")
	print("[->                             ERROR:")
	print("==============================================================================\n!")		
	if msg ~= nil then print(msg) end
	print("!\n==============================================================================")
	if otherMsg ~= nil then print(otherMsg) end
	local f, l = test.getSourceAndLine()
	print("     File:      "..f)
	print("     Line:      "..l)
	if msg ~= nil then print("     Message:   "..msg) end
	print("LUA BACKTRACE:")
	if backtraceSkipLevel == nil then
		backtraceSkipLevel = 2
	end
	DebugBacktrace(backtraceSkipLevel)		
	print("==============================================================================")
	print("                                 ERROR                                     <-]\n")
	
	--assert(false)
	exit()
end

--! use it like ug_assert(numPreRefs <= numRefs, "It must be choosen: numPreRefs <= numRefs")
--! @param condition the condition to assert
--! @param msg (optional) message to be printed if condition is not fulfilled
function ug_assert(condition, msg)
	if condition then
		return
	else
		ug_error(msg, "ASSERTION FAILED:", 3)		
	end
end

function ug_warning(t)
	print(t)
	err_log(t)
end

function ug_cond_warning(condition, text)
	if condition then
		ug_warning(text)
	end
end



--! @param backtraceLevel the number of levels to go up
--! for backtraceLevel = 0, it returns file and line of the
--! calling function. for  backtraceLevel = 1 the
--! file and line of the caller of the calling function and so on.
function util.GetLUAFileAndLine (backtraceLevel)
	local level = 2+backtraceLevel
	local info = debug.getinfo(level, "Sl")
	if not info then return "" end
	if info.what == "C" then   -- is a C function?
		return "C function"
	else
		return string.format("[%s]:%d", info.short_src, info.currentline)
	end
end
