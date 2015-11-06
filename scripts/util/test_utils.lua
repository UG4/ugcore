-- Copyright (c) 2011-2013:  G-CSC, Goethe University Frankfurt
-- Author: Martin Scherer
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.

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
		boost_require_impl(expression, msg, file, line)
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
