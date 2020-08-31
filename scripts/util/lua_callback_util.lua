
util = util or {}

-- As the LuaCallbackObserver class currently needs a globally defined function as callback, this is a workaround.


util.LuaCallbackHelper = { count = 0, instances = {} }

function LuaCallbackHelperCallback(step, time, dt, luaid)
    local cb = util.LuaCallbackHelper.instances[luaid]
    if cb.callback ~= nil then
        u = cb.CPPCallback:get_current_solution()
        local res = cb.callback(u, step, time, dt)
        if type(res) == "boolean" and res == false then
            return 0
        end
    end
    return 1
end

function util.LuaCallbackHelper:create(func)
    local cb = {}    
    cb.callback = func
    cb.CPPCallback = LuaCallbackObserver(util.LuaCallbackHelper.count)
    cb.CPPCallback:set_callback("LuaCallbackHelperCallback")

    util.LuaCallbackHelper.instances[util.LuaCallbackHelper.count] = cb
    util.LuaCallbackHelper.count = util.LuaCallbackHelper.count+1

    return cb
end
