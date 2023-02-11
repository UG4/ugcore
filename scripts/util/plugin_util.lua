-- Copyright (c) 2011-2020:  G-CSC, Goethe University Frankfurt
-- Author: Tim Schön
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

if util.LuaPluginsLoaded == nil then

    util.LuaPluginsLoaded = true
    --- Load Plugin Lua files
    loadedPlugins = GetLoadedPlugins()
    pluginPath = ug_get_root_path()..common:path_sep().."plugins"..common:path_sep()

    for i, name in ipairs(loadedPlugins) do
        local filename = pluginPath..name..common:path_sep().."lua"..common:path_sep().."lua-include.lua"
        if FileExists(filename) then
            print("Loading Plugin Lua Content from "..filename)
            ug_load_script(filename)
        end
    end

end