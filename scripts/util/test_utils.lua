--[[!
-- \defgroup scripts_util_test Tests Utility
-- \ingroup scripts_util
-- Utility functions for unit testing in Lua scripts.
-- \{
]]--

test = test or {}

function test.getSourceAndLine()
	-- gets the line of caller function
	local line = debug.getinfo(3).currentline
	-- gets caller script name, skip '@' char at beginning
	local file = string.sub(debug.getinfo(3).source, 2)
	
	return file, line
end



if 	boost_require_impl == nil or
	boost_check_impl == nil or
	boost_check_small_impl == nil or
	boost_message_impl == nil
then

	function test.require(expression, msg)
		file, line = test.getSourceAndLine()
		if not expression then 
			print("In File '" .. file .. "', Line " .. line ..": " .. msg); exit();
		end
	end

	function test.check(expression, msg)
		file, line = test.getSourceAndLine()
		if not expression then 
			print("In File '" .. file .. "', Line " .. line ..": " .. msg)
		end
	end
	
	function test.check_small(number, epsilon, msg)
		file, line = test.getSourceAndLine()
		if math.abs(number) < epsilon then
			print("In File '" .. file .. "', Line " .. line ..": " .. msg)	
		end
	end
	
	function test.message(msg) 
		file, line = test.getSourceAndLine()
		print("In File '" .. file .. "', Line " .. line ..": " .. msg)
	end

else

	function test.require(expression, msg)
		msg = msg or ""
		file, line = test.getSourceAndLine()
		-- break script execution if require_impl is false
		if boost_require_impl(expression, msg, file, line) == false then
			exit()
		end
	end

	function test.check(expression, msg)
		msg = msg or ""
		file, line = test.getSourceAndLine()
		
		-- pass result of expression to boost framework
		boost_check_impl(expression, msg, file, line)
	end
	
	function test.check_small(number, epsilon, msg)
		msg = msg or ""
		file, line = test.getSourceAndLine()
	
		-- pass numbers to boost framework
		boost_check_small_impl(number, epsilon, msg, file, line)
	end
	
	function test.message(msg) 
		file, line = test.getSourceAndLine()
	
		-- pass message to boost framework
		boost_message_impl(msg, file, line)
	end

end

-- end group scripts_util_test
--[[!
\}
]]--
