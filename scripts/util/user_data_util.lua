

--! functions user when '+' is called on an IPData (or a derived implementation)
function __ug__UserNumber_add(l,r)
	local rType = ug_class_name(r)
	local lType = ug_class_name(l)
	local rDim = -1
	local lDim = -1
	local rData = ""
	local lData = ""
	
	-- one of the arguments must be exported by registry
	if rType == "" and lType == "" then
		error("Internal Error in __ug__UserNumber_add: No UserData as operand")
	end

	-- check that userdata are IPData	
	if rType ~= "" then
		if not ug_is_base_class("IIPData", rType) then
			error("Error in '+': Can only add UserData, but summand is "..rType)
		end
		rDim = r:get_dim()
		rData = r:type()
	end

	if lType ~= "" then
		if not ug_is_base_class("IIPData", lType) then
			error("Error in '+': Can only add UserData, but summand is "..lType)
		end
		lDim = l:get_dim()
		lData = l:type()		
	end
		
	-- r is numeric
	if rType == "" then
		if tonumber(r) then
			if lData ~= "Number" then
				error("Error in '+': Cannot add Number and "..lData)
			end
			
			local linker =  _G["ScaleAddLinkerNumber"..lDim.."d"]()
			linker:add(l, r)
			return linker			
		else
			error("Error in '+': Summand must be numeric of UserData.")
		end
	end

	-- l is numeric
	if lType == "" then
		if tonumber(l) then
			if rData ~= "Number" then
				error("Error in '+': Cannot add Number and "..rData)
			end
			
			local linker =  _G["ScaleAddLinkerNumber"..rDim.."d"]()
			linker:add(l, r)
			return linker			
		else
			error("Error in '+': Summand must be numeric of UserData.")
		end
	end
		
	-- check for same dimesion
	if lDim ~= rDim then
		error("Error in '+': Dimensions of UserData does not match")
	end

	-- check for same data type
	if rData ~= lData then
		error("Error in '+': Data Type of UserData does not match")
	end
	
	-- add
	local linker =  _G["ScaleAddLinker"..rData..rDim.."d"]()
	linker:add(1.0, l)
	linker:add(1.0, r)
	return linker				
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
	end
end
end
end

-- Some usage info:
-- Name of instance of ug4-object: ug_class_name(obj)
-- Check if class is base class: ug_is_base_class("BaseClass", "DerivClass")
-- Check if dimension compiled in:  ug_dim_compiled(dim)
-- Returning metatable for a ug4-class: ug_get_metatable("ClassName")