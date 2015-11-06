-- Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
-- Author: Sebastian Reiter
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

-- This script contains some simple helper methods that assist in parsing
-- descriptor tables. Use-examples are e.g. solver_util_2.lua and
-- load_balancing_util_2.lua

util = util or {}
util.tableDesc = util.tableDesc or {}

function util.tableDesc.CondAbort(condition, message)
	if condition == true then
		print("ERROR in util.tableDesc: " .. message)
		exit()
	end
end

function util.tableDesc.IsPreset(desc)
	if type(desc) == "userdata" then
		return true
	else
		return false
	end
end

function util.tableDesc.ToNameAndDesc(descOrName)
	if type(descOrName) == "string" then
		return descOrName, nil
	elseif type(descOrName) == "table" then
		if descOrName.type then
			util.tableDesc.CondAbort(type(descOrName.type) ~= "string",
									 "'type' entry in table has to be of type 'string'")
			return descOrName.type, descOrName
		elseif descOrName.name then
			util.tableDesc.CondAbort(type(descOrName.name) ~= "string",
									 "'name' entry in table has to be of type 'string'")
			return descOrName.name, descOrName
		else
			print("Either 'type' or 'name' have to be specified in a descriptor-table")
			exit()
		end
	end
	util.tableDesc.CondAbort(true, "Invalid name or descriptor specified!")
end
