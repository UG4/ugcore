

--[[
b = GridFunction(approxSpace)
ev[1][1] = Eval(ev[1][1] / VecNorm(ev[1][1]))
print(ev[1][1])
print(VecProd(ev[1][1], ev[1][1]))
print(ev[1][1]*(ev[1][1]+ev[2][2]))
Assign(b, ev[1][1] + 3.0*ev[2][2])
SaveVectorForConnectionViewer(u, "u.vec")
SaveVectorForConnectionViewer(b, "b.vec")
]]--

function __ug__TemplateExpressions_add(a, b)
	return VecScaleAddClass(1.0, a, 1.0, b)
end
function __ug__TemplateExpressions_sub(a, b)	
	return VecScaleAddClass(1.0, a, -1.0, b)
end

function __ug__TemplateExpressions_mul(a, b)
	if type(a) == "number" or type(b) == "number" then
		return VecScaleAddClass(a, b)
	else
		local aa = a
		local bb = b
		if ug_class_group(a) == "VecScaleAddClass" then
			aa = a:eval()
		end
		if ug_class_group(b) == "VecScaleAddClass" then
			bb = b:eval()
		end
		return VecProd(aa, bb)
	end		
end

function __ug__TemplateExpressions_div(a, b)
	return VecScaleAddClass(1/b, a)		
end
-------------------------------------------------------

function set_user_data_overloads2(name)
	-- request metatable for the classname
	mt = ug_get_metatable(name)
	if mt == nil then return end
	
	-- set __add function in metatable
	mt.__add = _G["__ug__TemplateExpressions_add"]

	-- set __sub function in metatable
	mt.__sub = _G["__ug__TemplateExpressions_sub"]
	
	-- set __mul function in metatable
	mt.__mul = _G["__ug__TemplateExpressions_mul"]
	
	-- set __mul function in metatable
	mt.__div = _G["__ug__TemplateExpressions_div"]
end

-- add GridFunctionData
for k, class in ipairs({"GridFunction"}) do
	-- loop algebra
	for j, algebra in ipairs({"CPU1", "CPU2", "CPU3", "CPU4", "CPUVAR"}) do
		-- loop dimensions
		for dim =  1,3 do
			-- only set if dimension is compiled (otherwise metatable does not exist)
			if ug_dim_compiled(dim) and ug_algebra_compiled(algebra) then
				set_user_data_overloads2(class..dim.."d"..algebra)
			end
		end
	end
end

for j, algebra in ipairs({"CPU1", "CPU2", "CPU3", "CPU4", "CPUVAR"}) do
	if ug_algebra_compiled(algebra) then
		set_user_data_overloads2("Vector"..algebra)
	end
end

