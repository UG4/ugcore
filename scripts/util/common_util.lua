-- Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
-- Author: Stephan Grein
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
-- \defgroup scripts_util_common common Utility
-- \ingroup scripts_util
--
-- Utility functions for common operations
-- Questions? mailto:stephan.grein@gcsc.uni-frankfurt.de
-- \{
]]--

common = {}

--! unpacks a list of lists (trivial recursion to unpack all levels if needed)
local function prepare(s, ...)
   args = {}
   for i, v in ipairs(arg) do
      if (type(v) == "table") then
         for j, w in ipairs(v) do
            if (type(w) == "table") then
               for k, x in ipairs(w) do
                  args[#args+1] = x
               end
            else
               args[#args+1] = w
            end
         end
      else
         args[#args+1] = v
      end
   end

   return string.format(s, unpack(args))
end

--! emulates printf
function common:printf(s, ...) 
   print(prepare(s, ...))
end

--! emulates printf with newline appended
function common:printfn(s, ...)
   print(prepare(s .. "\n", ...))
end

--! emulates sprintf
function common:sprintf(s, ...)
   return prepare(s, ...)
end

--! emulates sprintfn
function common:sprintfn(s, ...)
   return prepare(s .. "\n" , ...)
end

--! tail
function common:tail(list, first)
    return { select(first or 10, unpack(list)) }
end

--! head
function common:head(list, last)
   return { table.remove({select(1, unpack(list))}, last) }
end

--! split filename by delimiter (defaults to slash) and takes the n-th group 
function string:split_path(sep, n)
   local sep, fields = sep or "/", {}
   local pat = string.format("([^%s]+)", sep)
   self:gsub(pat, function(z) fields[#fields+1] = z end)
   return fields[n or #fields]
end

-- end group scripts_util_common
--[[!
\}
]]--
