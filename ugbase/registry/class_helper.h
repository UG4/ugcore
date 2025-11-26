/*
 * Copyright (c) 2010-2013:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/**
 * \file class_helper.h
 *
 * \author Martin Rupp
 *
 * \date 20.10.2010
 *
 * ClassHierarchy implementation, GetClassHierarchy, FindClass
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG_BRIDGE__CLASS_HELPER__
#define __H__UG_BRIDGE__CLASS_HELPER__

#include <string>
#include "common/ug_config.h"

#include "registry.h"
#include "class.h"
#include "global_function.h"

namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

// ClassHierarchy
//--------------
/**
 * \brief Class Hierarchy Helper Class for UG Registry
 * This class stores class names and their subclasses
 * \sa GetClassHierarchy
 */
class UG_API ClassHierarchy
{
	public:
		ClassHierarchy() : bGroup(false) {}

		/**
		 * adds the class c to the class hierarchy by attaching it to its base
		 * hierarchy (base hierarchy taken from c->class_names()). automatically
		 * creates nonexisting base hierarchy.
		 * \param	c		Class to be inserted
		 */
		void insert_class(const IExportedClass &c);

		/**
		 * searches the hierarchy for the classname name.
		 * \param	name	class name to be searched
		 * \return nullptr if class name not found,
		 * 				otherwise ClassHierarchy with the class as base
		 * 				(find_class(name)->name() == name)
		 */
		ClassHierarchy *find_class(const char *name);

		bool operator < (const ClassHierarchy &other) const
		{
			return name < other.name;
		}

		void sort();

		std::string name;
		bool bGroup;
		std::vector<ClassHierarchy> subclasses;
};
/**
 * inits hierarchy with all classes of UGBridge
 */
UG_API void GetClassHierarchy(ClassHierarchy &hierarchy, const Registry &reg);

/**
 * Finds the class classname in the default ug registry and returns
 * IExportedClass pointer if found, otherwise nullptr
 */
UG_API std::string FunctionInfo(const Registry &reg, const char *functionname);

UG_API std::string FunctionInfo(const ExportedFunctionBase &thefunc, bool isConst=false,
                       const char *classname=nullptr, const char *highlightclassname=nullptr, bool bPrintHelp=false);
UG_API std::string ConstructorInfo(const ExportedConstructor &constr, const char *classname,
		const char *highlightclassname=nullptr);

UG_API std::string ClassHierarchyString(const Registry &reg, const char *classname);
UG_API std::string ClassInfo(const IExportedClass &c);
UG_API std::string ClassInfo(const Registry &reg, const char *classname);
UG_API std::string ClassUsageExact(const Registry &reg, const char *classname, bool OutParameters);
UG_API const ExportedFunction *FindFunction(const Registry &reg, const char *functionname);
UG_API std::string ParameterToString(const ParameterInfo &par, int i);
UG_API bool IsClassInParameters(const ParameterInfo &par, const char *classname);

// end group registry
/// \}

} // end namespace
} // end namespace

#endif
