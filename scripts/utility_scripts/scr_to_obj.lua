-- scaling analyzer
-- created by Sebastian Reiter
-- s.b.reiter@googlemail.com

ug_load_script("ug_util.lua")

local inFile = util.GetParam("-in", "")
local outFile = util.GetParam("-out", "")
local subsetName = util.GetParam("-subset ", "subset")


-- open files
local fi = io.open(inFile, "r")
if fi == nil then
	print("Invalid file specified for parameter -in: "..inFile)
	print("Aborting...")
	exit()
end

local fo = io.open(outFile, "w")
if fo == nil then
	print("Invalid file specified for parameter -out: "..outFile)
	print("Aborting...")
	exit()
end


-- write the object
fo:write("o "..subsetName.."\n")

-- this is used to match a number: (-?%d*%.*%d*).
-- No scientific format currently supported (feel free to imrove it).
local NUMMATCH = "(-?%d*%.*%d*)"

for line in fi:lines() do
-- try to match either bn, in or ie
	local v1, v2 = string.match(line, "bn%s+"..NUMMATCH.."%s+"..NUMMATCH)
	if v1 ~= nil then
		fo:write("v "..v1.." "..v2.." 0\n")
	else
		local v1, v2 = string.match(line, "in%s+"..NUMMATCH.."%s+"..NUMMATCH)
		if v1 ~= nil then
			fo:write("v "..v1.." "..v2.." 0\n")
		else
			local v1, v2, v3, v4 = string.match(line, "ie%s+(%d*)%s+(%d*)%s+(%d*)%s+(%d*)")
			if v1 ~= nil then
				fo:write("f "..(v1+1))
				if v2 ~= nil then fo:write(" "..(v2+1)) end
				if v3 ~= nil then fo:write(" "..(v3+1)) end
				if v4 ~= nil then fo:write(" "..(v4+1)) end
				fo:write("\n")
			end
		end
	end
end
