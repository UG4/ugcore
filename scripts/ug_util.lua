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

if print_all == nil then
function print_all(...)
	local la = GetLogAssistant()
	local opp = la:get_output_process()
	la:set_output_process(-1)
	print(unpack(arg))
	la:set_output_process(opp) 
end
end    
    

util._original_io_open = io.open

--! WARNING: Parallel File open is REALLY slow on clusters
--! this function overwrite io.open and prints a warning
--! if you use it on a core which is not 0
--! to remove this warning
--! 1. check if you want to open the file on ALL cores
--!  if not, use   if ProcRank()==0    open, write,close    end
--! 2. if you're really sure you want to do that, use io.open_ALL. 
function util.safe_io_open(filename, model)	
	if ProcRank() == 1 then
		print_all("--- WARNING: opening a file not from proc 0 may harm performance (see util.IOOpen) ! "..util.GetLUAFileAndLine(1).." ---")
	end
	return util._original_io_open(filename, model)
end

function io.open_ALL(filename, model)
	return util._original_io_open(filename, model)
end

util.FileDummy = {}
function util.FileDummy.write(...) end
function util.FileDummy.read(...) error("io.open_0 does not support read.") end
function util.FileDummy.close() end
--! io.open_ONE opens the file on exactly one core
--! all other cores get dummy file objects (FileDummy)
function io.open_ONE(filename, model, rootNode)
	if rootNode == nil then rootNode = 0 end
	if ProcRank() == rootNode then
		return util._original_io_open(filename, model)
	else
		return util.FileDummy
	end
end

io.open = util.safe_io_open


-- end group scripts_util
--[[!  
\} 
]]--
