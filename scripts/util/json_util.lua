util = util or {}

-- Load third party module
if IsDefinedUG_JSON() then

  local myPath = ug_get_current_path()
  print("JSON.UTIL:" ..myPath)
  package.path = package.path..";".. myPath.."../richards_app/?.lua"
  print(package.path)

  util.json = require("json")
end