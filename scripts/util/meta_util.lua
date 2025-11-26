-- Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
-- Author: Torbjörn Klatt
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
\file meta_util.lua
\defgroup scripts_util_meta Metaprogramming Utility
\ingroup scripts_util
\brief Meta Programming Functions for Lua Scripts
\{
]]--

util = util or {}

--[[!
\brief Create a function with named parameters
With this function it is easy to create Lua functions with default values for
their arguments and defining names for them to circumvent the position 
dependence of the arguments.
This is especially useful for functions with a lot of arguments.

\attention When calling with named parameters, you have to give a table with
  the named arguments as the only parameter to the created functions.

Use as:
\code
myfunction = util.CreateFancyFunction(
  {{"a"},{"b",7},{"c",5}},
  function(a, b, c)
    print(a, b, c)
  end
)

myfunction( 3 )            -->> prints: 3 7 5
myfunction( 4, 2 )         -->> prints: 4 2 5
myfunction{ a=3, c=1 }     -->> prints: 3 7 1
myfunction( { a=3, c=1 } ) -->> prints: 3 7 1
\endcode
\param[in] arg_def table defining the arguments of the function and their default values
\param[in] f the function itself
\return function with named parameters
\note Source: https://gist.github.com/stuartpb/975399
]]--
function util.CreateFancyFunction( arg_def, f )
  return function( args )
    local params = {}
    for i=1, #arg_def do
      local paramname = arg_def[i][1] -- the name of the first parameter to the function
      local default_value = arg_def[i][2]
			if #arg_def[i] > 1 then
			else
			end
      params[i] = args[i] or args[paramname] or default_value
    end
    return f( unpack( params, 1, #arg_def ) )
  end
end

--[[!
\}
]]--
