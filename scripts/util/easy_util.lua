util = util or {}

--! this function ensures that elements in T which are
--! not set (i.e. T[i]=nil) are set to T[i] = defaultT[i]
--! this way we can use T as a parameter Table
--! and all parameters not given are set to default values
--! e.g. SetDefaultParameters(T, {smoother="jac", maxLevels=2}) 
function SetDefaultParameters(T, defaultT)
	for i, v in pairs(defaultT) do
		if T[i] == nil then T[i] = v end
	end
end

--! checks in the table p that all elements of parameters
--! are non-nil. if some are nil, an error is thrown
function NeededParameters(p, parameters)
	for i, v in pairs(parameters) do
		if p[v] == nil then
			print("\ntable ")
			print(p)
			print("needs parameter \""..v.."\" !")
			ug_assert(false, "table needs parameter \""..v.."\" !")
		end	
	end
end

function util.EasyLoadGrid(p)
	NeededParameters(p, {"gridName", "requiredSubsets"})
	
	-- refinements:
	numPreRefs = util.GetParamNumber("-numPreRefs", p.numPreRefs or 0, "number of refinements before parallel distribution")
	numRefs    = util.GetParamNumber("-numRefs",    p.numRefs or 0, "number of refinements")
			
	if numPreRefs > numRefs then
		print("WARNING: numPreRefs = "..numPreRefs.." > numRefs = "..numRefs..". Setting now numPreRefs = numRefs = "..numRefs) 
	--	numPreRefs = numRefs 
	end
		
	dom = util.CreateAndDistributeDomain(p.gridName, numRefs, numPreRefs, p.requiredSubsets)
	
	print(dom:domain_info():to_string())
		
	local stretchX = util.GetParamNumber("-stretchX", 1.0)
	local stretchY = util.GetParamNumber("-stretchY", 1.0)
	if stretchX ~= 1.0 or stretchY ~= 1.0 then
		print("scaling domain by ("..stretchX..", "..stretchY..")") 
		ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))
	end
	
	local distortDomain = util.GetParamNumber("-distortDomain", 0.0)
	if distortDomain > 0.0 then
		distortDomain = math.pow(distortDomain, numRefs+1)
		print("distortDomain ="..distortDomain)
		RandomizeDomain(dom, distortDomain, distortDomain, distortDomain)
	end			
		
	return dom
end