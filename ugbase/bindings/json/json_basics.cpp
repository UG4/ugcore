

#include "json_basics.hpp"



// Config for solvers.
namespace ug {

#define _JSON_DEFAULTS_SOLVER_ R"(
#include "../scripts/util/solver.defaults.json"
)"

#define _JSON_SCHEMAS_SOLVER_ R"(
#include "../scripts/util/solver.defaults.json"
)"

const nlohmann::json json_predefined_defaults::solvers
		= nlohmann::json::parse(_JSON_DEFAULTS_SOLVER_);

const nlohmann::json json_predefined_schemas::solvers
		= nlohmann::json::parse(_JSON_SCHEMAS_SOLVER_);

#undef json_defaults_solver
#undef json_schemas_solver

};



