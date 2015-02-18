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
ug_load_script("util/debug_util.lua")
ug_load_script("util/command_line_util.lua")
ug_load_script("util/common_util.lua")

--------------------------------------------------------------------------------

--! returns the standard path at which grids are stored
function util.GetGridPath()
	return ug_get_data_path().."/grids/"
end

--! perhaps move this to bridge, and use function
function util.DirSeperator()
	if GetOperatingSystem() == "cygwin" then return "\\" else return "/" end
end

--! perhaps move this to file_util*
function util.GetBaseName(filename)
	
	for i=1,string.len(filename) do
		if string.char(string.byte(filename, i)) == util.DirSeperator() then
			return string.sub(filename, i+1)
		end 
	end
	return filename
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

    
------ PARALLEL FILE OPEN FUNCTIONS ------ 
-- the following code inserts checks to prevent unnecessary i/o --
  
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



-- end group scripts_util
--[[!  
\} 
]]--
