
// #define __UG4_PYBIND11_EXPERIMENTAL__
// #ifdef  __UG4_PYBIND11_EXPERIMENTAL__

// System
#include <string>

// UG4 core modules
#include "bridge/bridge.h"
#include "bindings/pybind/ug_pybind.h"

// this plugin
#include "ug_pybind.h"

// Expose registry to python.
PYBIND11_MODULE(ug4py, m)
{
	typedef ug::pybind::RegistryAdapter TRegistry;

	m.doc() = "My UG4 module";
	TRegistry reg(m);

	std::string grp("UG4");
	ug::pybind::RegisterStandardBridges(reg, grp);
}
// #endif
