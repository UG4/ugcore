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
