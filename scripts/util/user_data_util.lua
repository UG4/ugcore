

--! functions user when '+' is called on an IPData (or a derived implementation)
function __ug__UserNumber_sum(lScale, l, rScale, r)
	local rType = ug_class_name(r)
	local lType = ug_class_name(l)
	local rDim = -1
	local lDim = -1
	local Dim = -1
	local rData = ""
	local lData = ""
	local Data = ""
	
	---------------------------------------------------------
	-- Check type of operands (and read dim+type if IPData)
	---------------------------------------------------------
	-- one of the arguments must be exported by registry
	if rType == "" and lType == "" then
		error("Internal Error in __ug__UserNumber_add: No ug4 class as operand")
	end

	-- check that userdata are IPData	
	if rType ~= "" then
		if not ug_is_base_class("IIPData", rType) then
			error("Error in '+': Can only add UserData, but summand is "..rType)
		end
		rDim = r:get_dim()
		Dim = rDim
		rData = r:type()
		Data = rData
	end
	if lType ~= "" then
		if not ug_is_base_class("IIPData", lType) then
			error("Error in '+': Can only add UserData, but summand is "..lType)
		end
		lDim = l:get_dim()
		Dim = lDim
		lData = l:type()		
		Data = lData
	end
		
	---------------------------------------------------------
	-- Check match of types for operands
	---------------------------------------------------------
	
	-- both operands are IPData
	if rType ~= "" and lType ~= "" then
		-- check for same dimesion
		if lDim ~= rDim then
			error("Error in '+': Dimensions of UserData does not match")
		end
	
		-- check for same data type
		if rData ~= lData then
			error("Error in '+': Data Type of UserData does not match")
		end
	-- one operand is scalar value from lua
	else
		if rType == "" then
			if tonumber(r) then
				if lData ~= "Number" then
					error("Error in '+': Cannot add Number and "..lData)
				end
			else
				error("Error in '+': Summand must be scalar number or UserData.")
			end
		end
		if lType == "" then
			if tonumber(l) then
				if rData ~= "Number" then
					error("Error in '+': Cannot add Number and "..rData)
				end
			else
				error("Error in '+': Summand must be scalar number or UserData.")
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

--! functions used when '+' is called on an IPData (or a derived implementation)
function __ug__UserNumber_add(l,r)
	return __ug__UserNumber_sum(1.0, l, 1.0, r)
end

--! functions used when '-' is called on an IPData (or a derived implementation)
function __ug__UserNumber_sub(l,r)
	return __ug__UserNumber_sum(1.0, l, -1.0, r)
end

--------------------------------------------------------------------------------
-- Loop to set the __add functions for IPData
--------------------------------------------------------------------------------

-- loop some kinds of IPData implementation
for k, class in ipairs({"User", "ConstUser", "LuaUser", "ScaleAddLinker"}) do
-- loop some kind of data types
for k, type in ipairs({"Number", "Vector", "Matrix"}) do
-- loop dimensions
for dim =  1,3 do
	-- only set if dimension is compiled (otherwise metatable does not exist)
	if ug_dim_compiled(dim) then
		
		-- set __add function in metatable
		mt = ug_get_metatable(class..type..dim.."d")
		mt.__add = _G["__ug__UserNumber_add"]
	
		-- set __sub function in metatable
		mt = ug_get_metatable(class..type..dim.."d")
		mt.__sub = _G["__ug__UserNumber_sub"]
		
	end
end
end
end

-- Some usage info:
-- Name of instance of ug4-object: ug_class_name(obj) (returns "" if not a ug4 class)
-- Check if class is base class: ug_is_base_class("BaseClass", "DerivClass")
-- Check if dimension compiled in:  ug_dim_compiled(dim)
-- Returning metatable for a ug4-class: ug_get_metatable("ClassName")