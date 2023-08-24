util = util or {}

-- If UG_JSON is active, load third party module as util.json
if IsDefinedUG_JSON() then

  local jsonPath = ug_get_root_path().."/externals/JSONForUG4/json-lua/"
  package.path = package.path..";".. jsonPath.."?.lua"
  util.json = require("json")
  
end
