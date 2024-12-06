
// #define __UG4_PYBIND11_EXPERIMENTAL__
// #ifdef  __UG4_PYBIND11_EXPERIMENTAL__

// System
#include <string>

// UG4 core modules
#include "bridge/bridge.h"


// this plugin
#include "ug_pybind.h"
#include "python_user_data.h"
#include "python_views.h"

// Expose registry to python.
PYBIND11_MODULE(pyugcore, m)
{
	typedef ug::pybind::RegistryAdapter TRegistry;

	m.doc() = "My UG4 module";
	TRegistry reg(m);

	std::string grp("UG4");
	ug::pybind::RegisterStandardBridges(reg, grp);
	ug::pybind::RegisterPythonUserData(reg, grp);
	ug::pybind::RegisterPythonViews(reg, grp);
}
// #endif
