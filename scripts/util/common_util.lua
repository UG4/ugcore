--[[!
-- \defgroup scripts_util_common common Utility
-- \ingroup scripts_util
-- Utility functions for common operations
-- \{
]]--

common = {}

--! emulates printf
function common:printf(s, ...) 
   print(string.format(s, ...))
end

--! emulates printf with newline appended
function common:printfn(s, ...)
   print(string.format(s .. "\n", ...))
end

--! emulates sprintf
function common:sprintf(s, ...)
   return string.format(s, ...)
end

--! emulates sprintfn
function common:sprintfn(s, ...)
   return string.format(s .. "\n", ...)
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
