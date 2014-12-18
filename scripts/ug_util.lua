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
ug_load_script("util/vec_util.lua")
ug_load_script("util/gnuplot.lua")
ug_load_script("util/table_util.lua")
ug_load_script("util/time_step_util.lua")
ug_load_script("util/solver_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("util/domain_util.lua")
ug_load_script("util/math_util.lua")
ug_load_script("util/output_util.lua")
ug_load_script("util/checkpoint_util.lua")
-- ug_load_script("util/command_line_util.lua") is called below after the definition
-- of ug_assert. This should be changed! ug_assert should be in a separate file
-- that is included by command_line_util.lua


--------------------------------------------------------------------------------

--! use it like ug_assert(numPreRefs <= numRefs, "It must be choosen: numPreRefs <= numRefs")
--! @param condition the condition to assert
--! @param msg (optional) message to be printed if condition is not fulfilled
function ug_assert(condition, msg)
	if condition then
		return
	else
		print("\n")
		print("[->                             ERROR:")
		print("==============================================================================\n!")		
		if msg ~= nil then print(msg) end
		print("!\n==============================================================================")
		print("ASSERTION FAILED:")
		local f, l = test.getSourceAndLine()
		print("     File:      "..f)
		print("     Line:      "..l)
		if msg ~= nil then print("     Message:   "..msg) end
		print("LUA BACKTRACE:")
		DebugBacktrace()		
		print("==============================================================================")
		print("                                 ERROR                                     <-]\n")
		
		assert(false)
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

ug_load_script("util/command_line_util.lua")

--------------------------------------------------------------------------------

--! returns the standard path at which grids are stored
function util.GetGridPath()
	return ug_get_data_path().."/grids/"
end

--------------------------------------------------------------------------------
-- lua script functions
--------------------------------------------------------------------------------

function util.TableToTextLongHelper(indexPar, valuePar)
	local str=""
	if type(valuePar) == "table" then
		str = str..util.PrintTableHelperIntend .. tostring(indexPar)  .. " = {\n"
		util.PrintTableHelperIntend = util.PrintTableHelperIntend .. "  "
		
		for i,v in pairs(valuePar) do str = str..util.TableToTextLongHelper(i, v) end
		
		util.PrintTableHelperIntend = string.sub(util.PrintTableHelperIntend, 3)
		str = str..util.PrintTableHelperIntend .. "}\n"
	else
		if type(valuePar) == "string" or type(valuePar) == "number" then
			str = str..util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. valuePar .."\n" 
		elseif type(valuePar) == "boolean" then
			str = str..util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. tostring(valuePar) .."\n"
		else
			str = str..util.PrintTableHelperIntend .. " " .. tostring(indexPar) .. " = " .. tostring(valuePar) .."\n"
		end
	end
	return str
end

function util.TableToTextLong(tablePar)
	util.PrintTableHelperIntend = ""
	return util.TableToTextLongHelper("", tablePar)
end

--! to print tables
function util.PrintTable(tablePar)
	print(util.TableToTextLong(tablePar))
end

function util.TableToText(var)
	local out= ""
	local i
	local v
	if type(var) == "table" then
		out = out.." {"
		
		local count = 0
		for _ in pairs(var) do count = count + 1 end
		if count == #var then count = 0 else count = 1 end
		
		local bfirst = true		
		for i,v in pairs(var) do
			if bfirst then bfirst = false else out=out .. ", " end
			if count == 1 then out=out .. "["..tostring(i).."] = " end 
			out=out..util.TableToText(v)
			 
		end
		
		out = out.. "}"
	else
		out = out ..tostring(var)		
	end
	return out
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
    

util._original_io_open = util._original_io_open or io.open

function io.open_ALL(filename, model)
	ug_assert(filename ~= nil) 		
	return util._original_io_open(filename, model)
end
--! WARNING: Parallel File open is REALLY slow on clusters
--! this function overwrite io.open and prints a warning
--! if you use it on a core which is not 0
--! to remove this warning
--! 1. check if you want to open the file on ALL cores
--!  if not, use   if ProcRank()==0    open, write,close    end
--! 2. if you're really sure you want to do that, use io.open_ALL. 
function util.safe_io_open(filename, model)	
	if ProcRank() == 1 then
		ug_warning("--- WARNING: opening a file not from proc 0 may harm performance (see util.IOOpen) ! "..util.GetLUAFileAndLine(1).." ---")
	end
	return io.open_ALL(filename, model)
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

--! ParallelMaxMinAvg prints "min: (minimum), max: (maximum), avg: (average)"
--! for parallel different integers s
function util.ParallelMaxMinAvg(s)
	return "min: "..ParallelMin(s)..", max: "..ParallelMax(s)..", avg: ".. ParallelSum(s)/NumProcs()
end

_tostring = _tostring or tostring
function tostring(Val)
	if type(Val) == "table" then
   		return util.TableToTextLong(Val)
   	elseif type(Val) == "boolean" then
   		if Val then
   			return "true"
   		else
   			return "false"
   		end
   	else
   		return _tostring(Val)
   	end
end


-- end group scripts_util
--[[!  
\} 
]]--
