--------------------------------------------------------------------------------
--	tut01_lua_basics.lua
--
--	This tutorial will show you some of the basics of the lua syntax.
--	Detailed information can be found at
--	http://www.lua.org/manual/5.1/manual.html
--
--	Note that lua has some powerful features that go far beyond this tutorial.
--	A comprehensive and structured introduction is given at
--	http://www.lua.org/pil/index.html	
--
--	This short tutorial contains several small sections:
--	* Section 0: comments
--	* Section 1: variables and strings
--	* Section 2: control structures
--	* Section 3: tables (arrays and much more)
--	* Section 4: functions
--	* Section 5: scopes
--  * Section 6: standard libraries 
--------------------------------------------------------------------------------




--------------------------------------------------------------------------------
-- Section 0: comments
--------------------------------------------------------------------------------
-- This is for sure the shortest section of this tutorial.
-- You may already have realized that lines starting with a -- are regarded
-- as comments in lua. Of course a comment does not have to start at the
-- beginning of a line but may be appended to a regular code statement.




--------------------------------------------------------------------------------
-- Section 1: variables and strings
--------------------------------------------------------------------------------

-- lua variables do not have a type. A variable is declared by simply using it:
-- we're using the print command to output values to the console.
var = 3
print(var)
var = "now I'm a string"
print(var)

-- One can use the .. operator to concatenate strings
-- Note that known types are implicitly converted to strings on the fly
var1 = "A concatenated"
var2 = " message and 5 is "
var3 = 5
str = var1 .. var2 .. var3
print(str)

-- uninitialized variables are set to nil by default.
-- Note that nil is a special value defined by lua and can not be transformed
-- to a number. You can invalidate a variable by setting it to nil.
print(varUninitialized)


-- to convert a number in a string to a real number, simply use tonumber
a = "5"
b = tonumber(a)
print(b + b)


-- There are a lot of great tools to manipulate, split and tokenize strings
-- in lua, including comprehensive pattern matching. It is definitley
-- a good idea to take a closer look into this if you intend to do any of
-- the above.




--------------------------------------------------------------------------------
-- Section 2: control structures
--------------------------------------------------------------------------------

-- lua supports the common control structures like for, if etc.
-- Here are some examples:
if varUninitialized == nil then
	print("varUninitialized is still uninitialized.")
else 
	print("someone initialized varUninitialized!")
end

-- the negation operator in lua is ~
if varUninitialized ~= nil then
	print("varUninitialized is still initialized...")
else 
	print("varUninitialized is still uninitialized...")
end

-- We now use a for loop to assemble a string.
-- We iterate from 1 to 5. The default step size is 1.
str = ""
for i = 1, 5 do
	str = str .. i .. " "
end
print("assembled string in for loop: " .. str)

-- Now lets use a step size of 2
str = ""
for i = 1, 5, 2 do
	str = str .. i .. " "
end
print("assembled sting in for loop with step size 2: " .. str)

-- There are many more control structures in lua. Please check the lua reference
-- for more details.




--------------------------------------------------------------------------------
-- Section 3: tables (arrays and much more)
--------------------------------------------------------------------------------

-- Lua uses tables for everything that goes beyond basic types. Especially arrays
-- and associative containers are modeled using tables. Note that tables in lua
-- can even be used to model classes for object oriented programming. However
-- we won't go into that in this tutorial.

-- A simple array. Note that indices start at 1.
print("An array:")
a = {-1, -2, -3, -4, -5}
i = 1
while a[i] ~= nil do
	print(a[i])
	i = i + 1
end


-- Arrays can also be filled dynamically. Note that the variable first has to
-- contain a table.
-- Another interesting thing about lua tables (and thus arrays) is that different
-- entries can contain values of different types.
dynArray = {}
dynArray[1] = "a string "
dynArray[2] = "and a number: "
dynArray[3] = 13

-- we're using the ipairs method here, which allows us to iterate over all entries
-- of an array until the value nil is reached.
str = ""
for i, tval in ipairs(dynArray) do
	str = str .. tval
end

print("And another array: " .. str)


-- tables can of course also contain other tables:
dynArray[4] = {33, 35, 1100, "end of subtable"}


-- Note that arrays are not copied during assignment. Instead both variables
-- then point to the same table
tmpArray = dynArray
-- tmpArray and dynArray now both point to the same table.


-- tables can also be used as associative containers (indeed, this is what they are):
con = {}
con["name"] = "Sebastian"
con["occupation"] = "..."
con["age"] = "private!"

-- There's an interesting alternative syntax to access those fields
-- which are identified by a string:
print("Name: " .. con.name)
print("Occupation: " .. con.occupation)
-- of course the original array-like access is also valid:
print("Age: " .. con["age"])


-- Be sure to check the reference documentation on tables.

--------------------------------------------------------------------------------
-- Section 4: functions
--------------------------------------------------------------------------------

-- In this section we'll introduce functions. They can be used similarily to
-- functions in nearly every other programming language. However there is one
-- important difference: A function can return multiple return values.

function AddFunc(a, b)
	return a + b
end

print("48 and 89 sums up to " .. AddFunc(48, 89))


-- Now a method with two return values

function Scale(x, y, a)
	return x*a, y*a
end

xScaled, yScaled = Scale(3, 4, 2)
print("3 and 4 scaled by 2 results in " .. xScaled .. " and " .. yScaled)


-- NOTE: omitted parameters are set to nil:
function myprint(a, b)
	if(a ~= nil) then print(a) end
	if(b ~= nil) then print(a) end	
end

print("myprint(\"hello:\") ")
myprint("hello") -- omitted parameter 'b'
myprint() -- omitted both parameters

-- (this is also true for wrapped C/C++-functions)

-- functions can also be part of tables (this can be used to create 'namespaces')
mytable = mytable or {} 
function mytable.print_stuff()
	print("hello!")
end
mytable.print_stuff();

-- to 'extend namespaces' use
mytable = mytable or {}
-- because if you use
mytable = {}
-- mytable.print_stuff is gone.

--------------------------------------------------------------------------------
-- Section 5: scopes
--------------------------------------------------------------------------------

-- Note that lua has special scope rules. A variable used in a method is
-- by default a global variable and exists even after the method was terminated.
-- If you want to use temporary variables only visible inside the function
-- itself, you have to use the 'local' keyword:

-- Here an example with a non-local variable first used in a function:
function SetGlobalVar(a)
	globalVar = a
end

print("Value of globalVar before call to SetGlobalVar(77): ")
print(globalVar)
SetGlobalVar(77)
print("Value of globalVar after call to SetGlobalVar(77): ")
print(globalVar)


-- And now an example with a local variable
function SetLocalVar(a)
	local localVar = a
end

print("Value of localVar before call to SetLocalVar(33): ")
print(localVar)
SetLocalVar(33)
print("Value of localVar after call to SetLocalVar(33): ")
print(localVar)


--------------------------------------------------------------------------------
-- Section 6: standard libraries (http://www.lua.org/manual/5.1/manual.html#5)
--------------------------------------------------------------------------------

print("LUA Version is " .. _VERSION) -- might be nice to know

-- string manipulation (http://www.lua.org/manual/5.1/manual.html#5.4)
-- sometimes you want to format your strings like this
print("normal: " .. math.pi)
print("formatted: " .. string.format ("%.40f", math.pi) )

-- math functions (http://www.lua.org/manual/5.1/manual.html#5.6)
-- math.abs     math.acos    math.asin    math.atan    math.atan2
-- math.ceil    math.cos     math.deg     math.exp     math.floor
-- math.log     math.log10   math.max     math.min     math.mod
-- math.pow     math.rad     math.sin     math.sqrt    math.tan
-- math.frexp   math.ldexp   math.random  math.randomseed

-- i/o . (http://www.lua.org/manual/5.1/manual.html#5.7)
-- with io.open, you get access to files. io.open returns a handle which is actually a FILE* pointer
file = io.open("tut01_lua_basics_output.txt", "a")
-- use file:write to write data to the file
file:write("hello world")  -- this is the same as file:write(tostring(s))
io.close(file)

-- operating system facilities (http://www.lua.org/manual/5.1/manual.html#5.8)
-- we only cover time and date functions here:
tBefore = os.clock()
-- do something
for i = 1, 10000, 1 do
str = str .. i .. " "
end
tAfter = os.clock()
print("took " .. tAfter-tBefore .. " seconds!")

print("it is " .. os.date("%X on %x"))
filenamedate = os.date("y%ym%dd%d")
print("you can use " .. filenamedate .. " as prefix for your files")
