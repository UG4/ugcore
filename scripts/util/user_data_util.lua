--[[!
-- \file scripts/util/use_data_util.lua
-- 
-- Some usage info:
-- Name of instance of ug4-object: ug_class_name(obj) (returns "" if not a ug4 class)
-- Check if class is base class: ug_is_base_class("BaseClass", "DerivClass")
-- Check if dimension compiled in:  ug_dim_compiled(dim)
-- Returning metatable for a ug4-class: ug_get_metatable("ClassName")
!]]--

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
		rDim = r:get_dim()
		if not ug_is_base_class("IUserData"..rDim.."d", rType) then
			error("Error: Can only operate on UserData, but summand is "..rType)
		end
		Dim = rDim
		rData = r:type()
		Data = rData
	end
	if lType ~= "" then
		lDim = l:get_dim()
		if not ug_is_base_class("IUserData"..lDim.."d", lType) then
			error("Error: Can only operate on UserData, but summand is "..lType)
		end
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
	-- check for same dimesion
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
		
	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------
	
	-- both operands are UserData
	-- check for same data type
	if rType ~= "" and lType ~= "" then
		if rData ~= "Number" and lData ~= "Number" then
			error("Error in '*' or '/': One Data Type of UserData must be number")
		end
	-- one operand is scalar value from lua
	else
		if rType == "" then
			if not tonumber(r) then
				error("Error in '*' or '/': Operand must be scalar number or UserData.")
			end
		end
		if lType == "" then
			if not tonumber(l) then
				error("Error in '*' or '/': Operand must be scalar number or UserData.")
			end
		end		
	end	

	-- All checks passed: Can create name for Class to return
	local ScaleAddLinkerName = "ScaleAddLinker"..Data..Dim.."d"
	
	---------------------------------------------------------
	-- Multiply
	---------------------------------------------------------
	local linker =  _G[ScaleAddLinkerName]()
	if rData == "" or rData == "Number" then linker:add(r, l)
	else linker:add(l, r) end
	return linker				
end

--! functions user when '/' is called on an UserData (or a derived implementation)
function __ug__UserNumber_div(l, r)
	local rType, lType, rDim, lDim, Dim, rData, lData, Data = __ug__CheckUserDataArgType(r, l)
			
		
	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------
	
	-- both operands are UserData
	-- check for same data type
	if rType ~= "" and lType ~= "" then
		if rData ~= lData then
			error("Error in '/': Data Type of UserData does not match")
		end
	-- one operand is scalar value from lua
	else
		if rType == "" then
			if tonumber(r) then
				if lData ~= "Number" then
					error("Error in '/': Cannot divide Number and "..lData)
				end
			else
				error("Error in '/': Divisor must be scalar number or UserData.")
			end
		end
		if lType == "" then
			if tonumber(l) then
				if rData ~= "Number" then
					error("Error in '/': Cannot divide Number and "..rData)
				end
			else
				error("Error in '/': Dividend must be scalar number or UserData.")
			end
		end		
	end	

	-- All checks passed: Can create name for Class to return
	local InverseLinkerName = "InverseLinker"..Dim.."d"
	
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
	
--	if not tonumber(r) then
--		error("Error in '/': Currently divisor must be plain lua number")
--	else
--		return __ug__UserNumber_mul(l, (1/r))
--	end
end

--------------------------------------------------------------------------------
-- Function to Multiply/Devide UserData
--------------------------------------------------------------------------------

--! functions user when '^' is called on an UserData (or a derived implementation)
function __ug__UserNumber_pow(l, r)
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
for k, class in ipairs({"User", "ConstUser", "LuaUser", "ScaleAddLinker"}) do
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
for k, class in ipairs({"GridFunctionGradientData", "GridFunctionNumberData"}) do
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

