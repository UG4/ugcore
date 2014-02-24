--[[!
-- \defgroup scripts_util Lua Utility Scripts
-- \ingroup scripts
-- General Lua utility scripts for ug4.
-- \{
]]--

-- Create util namespace
util = util or {}

ug_load_script("util/meta_util.lua")
ug_load_script("util/test_utils.lua")
ug_load_script("util/domain_distribution_util.lua")
ug_load_script("util/stats_util.lua")
ug_load_script("util/user_data_util.lua")
ug_load_script("util/gnuplot.lua")
ug_load_script("util/table_util.lua")
ug_load_script("util/time_step_util.lua")
ug_load_script("util/solver_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("util/domain_util.lua")
ug_load_script("util/math_util.lua")
ug_load_script("util/command_line_util.lua")



--------------------------------------------------------------------------------

--! use it like ug_assert(numPreRefs <= numRefs, "It must be choosen: numPreRefs <= numRefs")
--! @param condition the condition to assert
--! @param msg message to be printed if condition is not fulfilled
function ug_assert(condition, msg)
	if condition then
		return
	else
		print("BACKTRACE:")
		DebugBacktrace()
		print("ASSERTION FAILED:")
		local f, l = test.getSourceAndLine()
		print("     File:      "..f)
		print("     Line:      "..l)
		print("     Message:   "..msg)
		assert(false)
	end
end

--------------------------------------------------------------------------------

--! returns the standard path at which grids are stored
function util.GetGridPath()
	return ug_get_data_path().."/grids/"
end

--------------------------------------------------------------------------------
-- lua script functions
--------------------------------------------------------------------------------

function util.PrintTableHelper(indexPar, valuePar)
	if type(valuePar) == "table" then
		print(util.PrintTableHelperIntend .. tostring(indexPar)  .. " = {")
		util.PrintTableHelperIntend = util.PrintTableHelperIntend .. "  "
		
		for i,v in pairs(valuePar) do util.PrintTableHelper(i, v) end
		
		util.PrintTableHelperIntend = string.sub(util.PrintTableHelperIntend, 3)
		print(util.PrintTableHelperIntend .. "}")
	else
		if type(valuePar) == "string" or type(valuePar) == "number" then
			print(util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. valuePar )
		elseif type(valuePar) == "boolean" then
			print(util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. tostring(valuePar) )
		else
			print(util.PrintTableHelperIntend .. "type(" .. tostring(indexPar) .. ") = " .. type(valuePar) )
		end
	end
end

--! to print tables
function util.PrintTable(tablePar)
	util.PrintTableHelperIntend = ""
	util.PrintTableHelper("", tablePar)
end



--! pairsSortedByKeys
--! the normal pairs(table) function returns elements unsorted
--! this function goes through elements sorted.
--! see http://www.lua.org/pil/19.3.html
--! use it like e.g. 
--! for name,arg in pairsSortedByKeys(util.args) do
--! f can be nil (= std sort)
function pairsSortedByKeys (t, f)
    local a = {}
    for n in pairs(t) do 
        table.insert(a, n) 
    end
    table.sort(a, f)
    local i = 0      -- iterator variable
    local function iter()   -- iterator function
    	i = i + 1
        if a[i] == nil then return nil
        else return a[i], t[a[i]]
        end
    end
	return iter
end

--------------------------------------------------------------------------------
-- basic functions missing lua
--------------------------------------------------------------------------------

--! adds writeln 
function writeln(...)
	write(...)
	write("\n")
end


function formatf(s, ...)
	return s:format(...)
end

function printf(s,...)
	print(formatf(...))
end

--! fsize
--! returns the filesize of a file (http://www.lua.org/pil/21.3.html)
--! @param file
--! @return filesize
function fsize (file)
	local current = file:seek()      -- get current position
    local size = file:seek("end")    -- get file size
    file:seek("set", current)        -- restore position
    return size
end

function bool2string(boolB)
	if boolB then
		return "true"
	else
		return "false"
	end
end

--------------------------------------------------------------------------------
-- list and free user data
--------------------------------------------------------------------------------

function ListUserDataInTable(t, name)
   	for n,v in pairs(t) do
      if type(v) == "userdata" then
		 print(name.."["..n.."]")
      end
  
	  if type(v) == "table" then
		if(n ~= "_G" and n ~= "io" and n ~= "package" and n ~= "gnuplot") then 
			ListUserDataInTable(v, name.."["..n.."]")
		end
  	  end
    end
end

--! Lists all user data (even in tables)
function ListUserData()
   	for n,v in pairs(_G) do
	   -- all userdata
	   if type(v) == "userdata" then
		 print(n)
	   end
    
	    -- userdata in table
		if type(v) == "table" then
			if(n ~= "_G" and n ~= "io" and n ~= "package" and n ~= "gnuplot") then 
				ListUserDataInTable(_G[n], n)
			end
		end
    end
end

function FreeUserDataInTable(t)
   	for n,v in pairs(t) do
      if type(v) == "userdata" then
      	 t[n] = nil
      end
      
      if type(v) == "table" then
		if(n ~= "_G" and n ~= "io" and n ~= "package" and n ~= "gnuplot") then
			FreeUserDataInTable(v)
		end
  	  end
      
    end
end

--! sets all userdata to nil (even in tables) and calls garbage collector
function FreeUserData()
   -- set user data to nil
   for n,v in pairs(_G) do
      if type(v) == "userdata" then
		 _G[n] = nil
      end
      
      if type(v) == "table" then
		if(n ~= "_G" and n ~= "io" and n ~= "package" and n ~= "gnuplot") then
 	     	FreeUserDataInTable(_G[n])
 	     end
      end
   end
   
   -- call garbage collector
   collectgarbage("collect")
end

--! 
--! @param pluginNamesList a list like {"amg", "d3f"} of plugins to check
function RequiredPlugins(pluginNamesList)
	local notLoadedNames = ""
	local cmakePluginString = ""
	for i,v in pairs(pluginNamesList) do
		if PluginLoaded(v) == false then
			notLoadedNames=notLoadedNames..v.." "
			cmakePluginString = cmakePluginString.." -D"..v.."=ON"
		end
	end	
	if notLoadedNames:len() > 0 then
		print("Plugin(s) needed but not loaded: "..notLoadedNames)		
		print("Please use \n   cmake "..cmakePluginString.." ..; make\nin your build directory to add the plugin(s).")
		exit()
	end
	
end

function AssertPluginsLoaded(pluginNamesList)
	RequiredPlugins(pluginNamesList)
end




-- end group scripts_util
--[[!  
\} 
]]--
