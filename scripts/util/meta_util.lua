--[[!
\file scripts/util/meta_util.lua
\brief Meta Programming Functions for Lua Scripts
]]--

util = util or {}

--[[!
\brief Create a function with namend parameters
With this function it is easy to create Lua functions with default values for
their arguments and defining names for them to circumvent the position 
dependentness of the arguments.
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
\return function with namend parameters
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
