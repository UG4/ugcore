--! this function creates a RAM disk on Mac.
--! it can be used for debug output, in case you want to have it faster or don't want to stress your SSD
--! if the volume is already there, nothing is created
--! keep in mind that RAM disks disappear after sleep or power off, so
--! only use it as a debug output.
--! @param name a std filename (don't use special characters or whitespace)
--! @param sizeMB the size of the RAM Disk in MB 
function util.createMacRAMDisk(name, sizeMB)
	ug_assert(GetOperatingSystem() == "apple", "only on apple systems")
	if ProcRank() == 0 then
		os.execute("if [ ! -e \"/Volumes/"..name.."\" ]; then\n"..
			"echo \"Creating ug4 RAM Disk with "..sizeMB.." MB..\"\n".. 
			"diskutil erasevolume HFS+ '"..name.."' `hdiutil attach -nomount ram://"..(1024*1024*sizeMB/512).."`\n"..
			"echo \"Access the RamDisk in /Volumes/"..name.."\"\n"..
			"fi\n")
	end
	SynchronizeProcesses()
	return "/Volumes/"..name	
end