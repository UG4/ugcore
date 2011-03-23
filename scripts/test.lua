function GetParam(name)
	for i = 1, ugargc-1 do
		if ugargv[i] == name then
			return ugargv[i+1]
		end
	end
	return nil; 
end


dim = GetParam("-dim")
if dim ~= nil then
	print("dim is "..dim.."\n")
else
	print("missing parameter -dim\n")
end

blocksize = GetParam("-blocksize")
if blocksize ~= nil then
	print("blocksize is "..blocksize.."\n")
else
	print("missing parameter -blocksize\n")
end