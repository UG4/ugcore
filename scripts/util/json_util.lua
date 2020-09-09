util = util or {}

-- Load third party module. 
-- If UG_JSON is active, we load into 'util.json'
if IsDefinedUG_JSON() then

  local jsonPath = ug_get_current_path().."../../../externals/JSONForUG4/json-lua/"
  
  package.path = package.path..";".. jsonPath.."?.lua"
  -- print("JSON.UTIL:" ..jsonPath)
  -- print(package.path)
  util.json = require("json")
end