-- Create util namespace
util = util or {}

ug_load_script("util/test_utils.lua")
ug_load_script("util/partition_maps.lua")
ug_load_script("util/stats_util.lua")
ug_load_script("util/user_data_util.lua")

--! returns the standard path at which grids are stored
function util.GetGridPath()
	return ug_get_data_path().."/grids/"
end

function util.GlobalRefineParallelDomain(domain)
	local dim = domain:get_dim()
	if dim == 1 then
		return GlobalRefineParallelDomain1d(domain)
	elseif dim == 2 then
		return GlobalRefineParallelDomain2d(domain)
	elseif dim == 3 then
		return GlobalRefineParallelDomain3d(domain)
	end
	return false
end

--------------------------------------------------------------------------------
-- User Data utils
--------------------------------------------------------------------------------

--! creates a Const User Matrix 2d 
function util.CreateConstUserMatrix2d(m00, m01, m10, m11)
	local mat = ConstUserMatrix2d()
	mat:set_entry(0, 0, m00)
	mat:set_entry(0, 1, m01)
	mat:set_entry(1, 0, m10)
	mat:set_entry(1, 1, m11)	
	return mat
end

--! creates a Const User Matrix 3d
function util.CreateConstUserMatrix3d(m00, m01, m02, m10, m11, m12, m20, m21, m22)
	local mat = ConstUserMatrix3d()
	mat:set_entry(0, 0, m00)
	mat:set_entry(0, 1, m01)
	mat:set_entry(0, 2, m02)
	mat:set_entry(1, 0, m10)
	mat:set_entry(1, 1, m11)
	mat:set_entry(1, 2, m12)
	mat:set_entry(2, 0, m20)
	mat:set_entry(2, 1, m21)
	mat:set_entry(2, 2, m22)	
	return mat
end

--! creates a Const User Vector 2d
function util.CreateConstUserVector2d(v0, v1)
	local vec = ConstUserVector2d()
	vec:set_entry(0, v0)
	vec:set_entry(1, v1)
	return vec
end

--! creates a Const User Vector 3d
function util.CreateConstUserVector3d(v0, v1, v2)
	local vec = ConstUserVector3d()
	vec:set_entry(0, v0)
	vec:set_entry(1, v1)
	vec:set_entry(2, v2)
	return vec
end

--------------------------------------------------------------------------------
-- Subset utils
--------------------------------------------------------------------------------

--! util.CheckSubsets
--! checks if all required subsets are contained in the SubsetHandler
--! @param dom Domain
--! @param neededSubsets List of subsets the SubsetHandler must contain
--! @return true if all subsets are contained, false else
function util.CheckSubsets(dom, neededSubsets)
	sh = dom:subset_handler()
	for i, tval in ipairs(neededSubsets) do
		if sh:get_subset_index(tval) == -1 then
			print("Domain does not contain subset '"..tval.."'.")
			return false
		end
	end
	
	return true
end


function util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

	-- create Instance of a Domain
	local dom = Domain()
	
	-- load domain
	LoadDomain(dom, gridName)
	
	-- create Refiner
	if numPreRefs > numRefs then
		print("numPreRefs must be smaller than numRefs. Aborting.");
		exit();
	end
	
	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner = GlobalDomainRefiner(dom)
	
	-- Performing pre-refines
	for i=1,numPreRefs do
		refiner:refine()
	end
	
	-- Distribute the domain to all involved processes
	if DistributeDomain(dom) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	
	-- Perform post-refine
	for i=numPreRefs+1,numRefs do
		refiner:refine()
	end
	
	-- Now we loop all subsets an search for it in the SubsetHandler of the domain
	if neededSubsets ~= nil then
		if util.CheckSubsets(dom, neededSubsets) == false then 
			print("Something wrong with required subsets. Aborting.");
			exit();
		end
	end
	
	--clean up
	delete(refiner)
	
	-- return the created domain
	return dom
end

--------------------------------------------------------------------------------
-- some auxiliary functions
--------------------------------------------------------------------------------
--! function returns true if the number is a power of two
function util.IsPowerOfTwo(n)
	local number compare = 1

	while (compare < n) do
		compare = compare*2
	end

	return compare==n
end

--! function returns true if the number is a natural number
function util.IsNaturalNumber(n)
	if n-math.floor(n) == 0 then
		return true
	else
		return false
	end
end

--! function to factorise number which has to be a power of 2 in two factors
--! which differ at most by a factor of 2 and returns both
--! (first the smaller one, then the larger one).
function util.FactorizeInPowersOfTwo(n)
	if not util.IsPowerOfTwo(n) then
		print("Number to factorise must be a power of 2. Aborting.")
		exit()
	end

	local number firstFactor = n
	local number secFactor = 1

	while (firstFactor > 2*secFactor) do
		firstFactor = firstFactor/2
		secFactor = secFactor*2
	end

	return secFactor, firstFactor
end

--------------------------------------------------------------------------------
-- Command line functions
--------------------------------------------------------------------------------

util.args = util.args or {}
util.argDescription = util.argDescription or {}
util.argsDefault = util.argsDefault or {}
util.argUsed = util.argUsed or {}

--! util.GetParam
--! returns parameter in ugargv after ugargv[i] == name
--! @param name parameter in ugargv to search for
--! @param return_if_unavailable returned value if 'name' is not present (default nil)
--! @param description description for 'name' (default nil)
--! @return parameter in ugargv after ugargv[i] == name or return_if_unavailable if 'name' was not present
function util.GetParam(name, return_if_unavailable, description)
	local i
	if description or util.argDescription[name] == nil then
		util.argDescription[name]=" (string) "..(description or "")
	end
	util.argsDefault[name]=return_if_unavailable	
	for i = 1, ugargc-1 do
		if ugargv[i] == name then			
			util.argUsed[i]=true
			util.argUsed[i+1]=true
			util.args[name] = ugargv[i+1]			
			return ugargv[i+1]
		end
	end
	util.args[name] = " (string)"
	return return_if_unavailable; 
end


--! util.GetParamNumber
--! use with CommandLine to get option, like -useAMG
--! if parameter is not a number, returns return_if_unavailable
--! @param name parameter in ugargv to search for
--! @param return_if_unavailable returned value if 'name' is not present (default nil)
--! @param description description for 'name' (default nil)
--! @return the number after the parameter 'name' or return_if_unavailable if 'name' was not present/no number
function util.GetParamNumber(name, return_if_unavailable, description)
	if description or util.argDescription[name] == nil then
		util.argDescription[name]=" (number) "..(description or "")
	end
	local param = util.GetParam(name, nil, nil)
	
	util.argsDefault[name]=return_if_unavailable
	if param == nil then
		util.args[name]=" (number)" 
		return return_if_unavailable
	else
		local number = tonumber(param)
		if number == nil then
			print("WARNING: Parameter "..name.." is not a number, using "..return_if_unavailable.." instead\n")
			util.args[name]=" (number) " 
			return return_if_unavailable
		else
			return number
		end
	end	
end

--! util.HasParamOption
--! use with CommandLine to get option, like -useAMG
--! @param name option in argv to search for
--! @param description description for 'name' (default nil)
--! @return true if option found, else false
function util.HasParamOption(name, description)
	util.argsDefault[name]="false"
	if description or util.argDescription[name] == nil then
		util.argDescription[name]=" [option] "..(description or "")
	end
	for i = 1, ugargc do
		if ugargv[i] == name then
			util.argUsed[i]=true
			util.args[name] = "true"
			return true
		end		
	end	
	util.args[name]=" [option]"	
	return false 
end

--! returns all arguments from the command line
function util.GetCommandLine()
	local pline = ""
	for i=1, ugargc do
		pline = pline..ugargv[i].." "
	end
	return pline
end

--! lists all the command line arguments which where used or could
--! have been used with util.GetParam, util.GetParamNumber and util.HasParamOption
function util.PrintArguments()
	local pUsedLine=""
	local pOtherLine=""
	for name,value in pairsSortedByKeys(util.args) do
		if string.sub(value,1,1) ~= " " then
			pUsedLine=pUsedLine..name.." = "..value.."\n"
		else
			pOtherLine=pOtherLine..name..value..", "
		end	
	end
	if pUsedLine ~= "" then
		print("Used arguments:\n"..pUsedLine)		
	end	
	if pOtherLine ~= "" then
		print("Available arguments:\n"..string.sub(pOtherLine, 1, string.len(pOtherLine)-2).."\n")
	end	
end

--! prints out the description to each GetParam-parameter so far called
function util.PrintHelp()
	local length=0
	for name,value in pairs(util.argDescription) do
		length = math.max(string.len(name), length)
	end	
	for name,value in pairsSortedByKeys(util.argDescription) do
		print(util.adjuststring(name, length, "l").." : "..value..". default="..(util.argsDefault[name] or ""))
	end
end

--! calls util.PrintHelp() and exits if HasParamOption("-help")
--! @param optional description of the script
function util.CheckAndPrintHelp(desc)
	if util.HasParamOption("-help", "print this help") then
		if desc ~= nil then print(desc) end
		print()			
		util.PrintHelp()
		exit()
	end
end

--! lists all command line arguments which were provided but could not be used.
function util.PrintIgnoredArguments()
	local pline = ""
	for i=1, ugargc do
		if (util.argUsed == nil or util.argUsed[i] == nil) and 
			string.sub(ugargv[i], 1,1) == "-" then
			local imin=10
			local namemin=""
			for name in pairs(util.args) do
				if string.sub(util.args[name],1,1) == " " then
					local d = LevenshteinDistance(name, ugargv[i])
					if d < imin then
						imin = d
						namemin = name
					end
				elseif name == ugargv[i] then
					imin = 0
				end
			end
			if imin == 0 then
				pline=pline..ugargv[i].." is a doubled argument! Already set to "..util.args[ugargv[i]]..".\n"
			elseif imin < string.len(ugargv[i])/2 then
				pline=pline..ugargv[i].." (did you mean "..namemin..util.args[namemin].."?)\n"
			else
				pline=pline..ugargv[i].." "				
			end
		end
	end
	if pline ~= "" then
		print("Ignored arguments:\n"..pline.."\n")
	end
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
		print(util.PrintTableHelperIntend .. tostring(indexPar) .. " = " .. valuePar)
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
function pairsSortedByKeys (t, f)
	local a = {}
    for n in pairs(t) do table.insert(a, n) end
    table.sort(a, f)
    local i = 0      -- iterator variable
    local iter = function ()   -- iterator function
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


formatf = function(s, ...)
	return s:format(...)
end

printf = function(s,...)
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
		if(n ~= "_G" and n ~= "io" and n ~= "package") then 
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
			if(n ~= "_G" and n ~= "io" and n ~= "package") then 
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
		if(n ~= "_G" and n ~= "io" and n ~= "package") then
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
		if(n ~= "_G" and n ~= "io" and n ~= "package") then
 	     	FreeUserDataInTable(_G[n])
 	     end
      end
   end
   
   -- call garbage collector
   collectgarbage("collect")
end

--! 
--! @param pluginNamesList a list like {"amg", "d3f"} of plugins to check
function AssertPluginsLoaded(pluginNamesList)
	local notLoadedNames = ""
	local cmakePluginString = ""
	for i,v in pairs(pluginNamesList) do
		if PluginLoaded(v) == false then
			notLoadedNames=notLoadedNames..v.." "
			cmakePluginString = cmakePluginString.." -D"..v
		end
	end	
	if notLoadedNames:len() > 0 then
		print("Plugin(s) needed but not loaded: "..notLoadedNames)		
		print("Please use \n   cmake "..cmakePluginString.." ..; make\nin your build directory to add the plugin(s).")
		exit()
	end
	
end

			
