-- Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
-- Author: Andreas Vogel
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
-- \file user_data_util.lua
-- \defgroup scripts_util_userdata UserData Utility
-- \ingroup scripts_util
-- Some usage info:
-- FreeUserDataInTable, FreeUserData.
--
-- Name of instance of ug4-object: ug_class_name(obj) (returns "" if not a ug4 class)
-- Check if class is base class: ug_is_base_class("BaseClass", "DerivClass")
-- Check if dimension compiled in:  ug_dim_compiled(dim)
-- Returning metatable for a ug4-class: ug_get_metatable("ClassName")
-- \{
!]]--

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



function __ug__CheckUserDataArgType(r, l)
	local rType = ug_class_name(r)
	local lType = ug_class_name(l)
	local rDim = -1
	local lDim = -1
	local Dim = -1
	local rData = ""
	local lData = ""
	local Data = ""

	---------------------------------------------------------
	-- Check type of operands (and read dim+type if UserData)
	---------------------------------------------------------
	-- one of the arguments must be exported by registry
	if rType == "" and lType == "" then
		error("Internal Error: No ug4 class as operand")
	end

	-- check that userdata are UserData
	if rType ~= "" then
		if not ug_is_base_class("UserDataInfo", rType) then
			error("Error: Can only operate on UserData, but summand is not derived from it.")
		end
		rDim = r:get_dim()
		Dim = rDim
		rData = r:type()
		Data = rData
	end
	if lType ~= "" then
		if not ug_is_base_class("UserDataInfo", lType) then
			error("Error: Can only operate on UserData, but summand is not derived from it.")
		end
		lDim = l:get_dim()
		Dim = lDim
		lData = l:type()
		Data = lData
	end

	if rData ~= "" and lData ~= "" then
		if rData == "Number" then
			Data = lData
		else
			Data = rData
		end
	end

	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------

	-- both operands are UserData
	-- check for same dimension
	if rType ~= "" and lType ~= "" then
		if lDim ~= rDim then
			error("Error: Dimensions of UserData does not match")
		end
	end

	return rType, lType, rDim, lDim, Dim, rData, lData, Data
end

--------------------------------------------------------------------------------
-- Function to add/subtract UserData
--------------------------------------------------------------------------------

--! functions user when '+/-' is called on an UserData (or a derived implementation)
function __ug__UserNumber_sum(lScale, l, rScale, r)
	local rType, lType, rDim, lDim, Dim, rData, lData, Data = __ug__CheckUserDataArgType(r, l)

	if l == nil or r == nil then
		error("Error in '+' or '-': Added a nil value (possibly uninitialized value?)")
	end

	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------

	-- both operands are UserData
	-- check for same data type
	if rType ~= "" and lType ~= "" then
		if rData ~= lData then
			error("Error in '+' or '-': Data Type of UserData does not match")
		end
	-- one operand is scalar value from lua
	else
		if rType == "" then
			if tonumber(r) then
				if lData ~= "Number" then
					error("Error in '+' or '-': Cannot add Number and "..lData)
				end
			else
				error("Error in '+' or '-': Summand must be scalar number or UserData.")
			end
		end
		if lType == "" then
			if tonumber(l) then
				if rData ~= "Number" then
					error("Error in '+' or '-': Cannot add Number and "..rData)
				end
			else
				error("Error in '+' or '-': Summand must be scalar number or UserData.")
			end
		end
	end

	-- All checks passed: Can create name for Class to return
	local ScaleAddLinkerName = "ScaleAddLinker"..Data..Dim.."d"

	---------------------------------------------------------
	-- a) Case: An operand is ScaleAddLinker
	---------------------------------------------------------
	-- left operand is linker
	if lType == ScaleAddLinkerName then
		-- copy linker and add also other operand
		local linker =  _G[ScaleAddLinkerName](l)
		linker:add(rScale, r)
		return linker
	end
	-- right operand is linker and add-operation
	-- NOTE: if right operand is subtracted (i.e. rScale == -1) we must use
	--       the general case
	if rType == ScaleAddLinkerName then
	 	if rScale == 1.0 then
			-- copy linker and add also other operand
			local linker =  _G[ScaleAddLinkerName](r)
			linker:add(lScale, l)
			return linker
		end
	end

	---------------------------------------------------------
	-- b) Case: No operand is ScaleAddLinker
	---------------------------------------------------------
	local linker =  _G[ScaleAddLinkerName]()
	linker:add(lScale, l)
	linker:add(rScale, r)
	return linker
end

--! functions used when '+' is called on an UserData (or a derived implementation)
function __ug__UserNumber_add(l,r)
	return __ug__UserNumber_sum(1.0, l, 1.0, r)
end

--! functions used when '-' is called on an UserData (or a derived implementation)
function __ug__UserNumber_sub(l,r)
	return __ug__UserNumber_sum(1.0, l, -1.0, r)
end

--------------------------------------------------------------------------------
-- Function to Multiply/Devide UserData
--------------------------------------------------------------------------------

--! functions user when '*' is called on an UserData (or a derived implementation)
function __ug__UserNumber_mul(l, r)
	local rType, lType, rDim, lDim, Dim, rData, lData, Data = __ug__CheckUserDataArgType(r, l)

	if l == nil or r == nil then
		error("Error in '*' : nil value (possibly uninitialized value?)")
	end
	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------

	-- if operand is not UserData it must be a number
	if rType == "" then
		if not tonumber(r) then
			error("Error in '*': Operand must be scalar number or UserData.")
		end
	end
	if lType == "" then
		if not tonumber(l) then
			error("Error in '*': Operand must be scalar number or UserData.")
		end
	end

	-- if plain lua-numbers, just multiply
	if rType == "" and lType == "" then
		return l*r
	end

	-- if one operand is a scalar or number
	if rData == "Number" or lData == "Number" or rData == "" or lData == "" then
		local ScaleAddLinkerName = "ScaleAddLinker"..Data..Dim.."d"
		local linker =  _G[ScaleAddLinkerName]()
		if rData == "" or rData == "Number" then linker:add(r, l)
		else linker:add(l, r) end
		return linker
	end

	if lData == "Matrix" and rData == "Vector" then
		local ScaleAddLinkerName = "ScaleAddLinkerVectorMatrix"..Dim.."d"
		local linker =  _G[ScaleAddLinkerName]()
		linker:add(l, r)
		return linker
	end

	if lData == "Vector" and rData == "Vector" then
		local ScaleAddLinkerName = "ScaleAddLinkerVectorVector"..Dim.."d"
		local linker =  _G[ScaleAddLinkerName]()
		linker:add(l, r)
		return linker
	end

	error("Error in '*': Cannot multiply "..lData.."*"..rData)
end

--! functions user when '/' is called on an UserData (or a derived implementation)
function __ug__UserNumber_div(l, r)
	local rType, lType, rDim, lDim, Dim, rData, lData, Data = __ug__CheckUserDataArgType(r, l)

	if l == nil or r == nil then
		error("Error in '/': nil value (possibly uninitialized value?)")
	end

	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------

	-- check r is number
	if rType ~= "" and rData ~= "Number" then
		error("Error in '/': Cannot divide by "..rData)
	end

	-- case: dividing by lua number
	if rType == "" then
		if not(tonumber(r)) then
			error("Error in '/': Cannot only divide by Number")
		end
		return __ug__UserNumber_mul(l, (1/r))
	end

	local InverseLinkerName = "InverseLinker"..Dim.."d"
	-- case: number by user-data-number
	if lType == "" or lData == "Number" then
		---------------------------------------------------------
		-- a) Case: An operand is InverseLinker
		---------------------------------------------------------
		-- left operand is linker
		if lType == InverseLinkerName then
			-- copy linker and add also other operand
			local linker =  _G[InverseLinkerName](l)
			linker:divide(1.0, r)
			return linker
		end

		-- right operand is linker
		if rType == InverseLinkerName then
		 		local linker =  _G[InverseLinkerName](r)
				linker:divide(l, 1.0)
				return linker
		end

		---------------------------------------------------------
		-- b) Case: No operand is InverseLinker
		---------------------------------------------------------
		local linker =  _G[InverseLinkerName]()
		linker:divide(l, r)
		return linker
	end

	-- case: user-data by user-data-number
	local linker = _G[InverseLinkerName]()
	linker:divide(1, r)
	return __ug__UserNumber_mul(l, linker)

end

--------------------------------------------------------------------------------
-- Function to Multiply/Devide UserData
--------------------------------------------------------------------------------

--! functions user when '^' is called on an UserData (or a derived implementation)
function __ug__UserNumber_pow(l, r)
	if l == nil or r == nil then
		error("Error in '^': nil value (possibly uninitialized value?)")
	end

	-- check that exponent is integer
	if not tonumber(r) or not (r%1==0) then
		error("Error in '^': Currently exponent must be plain lua integer")
	else
		-- check that type is number
		if ug_class_name(l) == "" then
			error("Error in '^': base must be UserData.")
		else
			if l:type() ~= "Number" then
				error("Error in '^': base must be a number-UserData.")
			end
		end


		-- case 1
		if r == 1 then return l end

		-- case > 1
		if r > 1 then
			local value = l
			for i = 2,r do
			value = __ug__UserNumber_mul(value, l)
			end
			return value
		end
		if r < 0 then
			local value = __ug__UserNumber_div(1.0, l)
			t = -r
			if r < -1 then
				for i=2,t do
				value = __ug__UserNumber_div(value, l)
				end
			end
			return value
		end
	end
end

--------------------------------------------------------------------------------
-- Loop to set the __add functions for UserData
--------------------------------------------------------------------------------

function set_user_data_overloads(name)
	-- request metatable for the classname
	mt = ug_get_metatable(name)
	if mt == nil then return end

	-- set __add function in metatable
	mt.__add = _G["__ug__UserNumber_add"]

	-- set __sub function in metatable
	mt.__sub = _G["__ug__UserNumber_sub"]

	-- set __mul function in metatable
	mt.__mul = _G["__ug__UserNumber_mul"]

	-- set __div function in metatable
	mt.__div = _G["__ug__UserNumber_div"]

	-- set __pow function in metatable
	mt.__pow = _G["__ug__UserNumber_pow"]
end

-- loop some kinds of UserData implementation
for k, class in ipairs({"CplUser", "ConstUser", "LuaUser", "ScaleAddLinker", "LuaUserFunction"}) do
	-- loop some kind of data types
	for k, type in ipairs({"Number", "Vector", "Matrix"}) do
		-- loop dimensions
		for dim =  1,3 do
			-- only set if dimension is compiled (otherwise metatable does not exist)
			if ug_dim_compiled(dim) then
				set_user_data_overloads(class..type..dim.."d")
			end
		end
	end
end

-- add GridFunctionData
for k, class in ipairs({"GridFunctionGradientData", "GridFunctionNumberData", "EffectiveUserTensor"}) do
	-- loop algebra
	for j, algebra in ipairs({"CPU1", "CPU2", "CPU3", "CPU4", "CPUVAR"}) do
		-- loop dimensions
		for dim =  1,3 do
			-- only set if dimension is compiled (otherwise metatable does not exist)
			if ug_dim_compiled(dim) and ug_algebra_compiled(algebra) then
				set_user_data_overloads(class..dim.."d"..algebra)
			end
		end
	end
end

for dim =  1,3 do
	-- only set if dimension is compiled (otherwise metatable does not exist)
	if ug_dim_compiled(dim) then
		set_user_data_overloads("DarcyVelocityLinker"..dim.."d")
		set_user_data_overloads("BinghamViscosityLinker"..dim.."d")
		set_user_data_overloads("UserTensor"..dim.."d")
		set_user_data_overloads("UserTensorEntryAdapter"..dim.."d")
		set_user_data_overloads("UserTensorColumnAdapter"..dim.."d")
    set_user_data_overloads("UserTensorMatrix"..dim.."d")
  set_user_data_overloads("TensorEntryAdapter"..dim.."d")
  set_user_data_overloads("TensorColumnAdapter"..dim.."d")
	end
end

--[[!
\}
]]--
