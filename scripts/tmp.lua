--	A small script to test SmartPtrs in lua.
sp = SmartTestImpl()

SmartTestFunc(sp)
ConstSmartTestFunc(sp)

sp = ConstSmartTestImpl()
collectgarbage("collect")

ConstSmartTestFunc(sp)

sp = nil
collectgarbage("collect")


