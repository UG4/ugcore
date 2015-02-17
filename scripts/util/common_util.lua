--[[!
-- \defgroup scripts_util_common common Utility
-- \ingroup scripts_util
--
-- Utility functions for common operations
-- Questions? mailto: stephan.grein@gcsc.uni-frankfurt.de
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
