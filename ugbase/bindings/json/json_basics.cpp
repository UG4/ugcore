

#include "json_basics.hpp"
#include "common/util/path_provider.h"
// Config for solvers.
namespace ug {

#define _JSON_DEFAULTS_SOLVER_ R"(\
#include "../scripts/util/solver.defaults.json"\
)"

#define _JSON_SCHEMAS_SOLVER_ R"(\
#include "../scripts/util/solver.defaults.json"\
)"
// Hack for gcc compiler bug where no raw string literals over multiple lines in macros are allowed
const char *defaults =
#include "../scripts/util/solver.defaults.json"
;

const nlohmann::json json_predefined_defaults::solvers
		= nlohmann::json::parse(defaults);//= nlohmann::json::parse(solver_defaults_path);

const nlohmann::json json_predefined_schemas::solvers
        = nlohmann::json::parse(defaults);//= nlohmann::json::parse(solver_defaults_path);


#undef json_defaults_solver
#undef json_schemas_solver

};



