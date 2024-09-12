/*
 * Copyright (c) 2024: Goethe University Frankfurt
 * Author: Arne Naegel
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


#pragma once

#include <type_traits>
#include <string>
#include <nlohmann/json.hpp>

#include "common/util/smart_pointer.h"
#include "common/log.h"

/*
 * For each class, we must supply the following:
 *
 * - is_json_constructible = true
 * - json_defaults = <valid JSON>
 * - json_schema = <valid JSON schema >
 *
 * Optional defines:
 *
 * - json_assignment::from_json
 *
 * */


#define UG_HAS_JSON_BASICS 1 // temporary

namespace ug {

//! Marks an object, with is constructible via JSON.
/*! DEPRECATED */
struct JSONConstructible {};


/*****************************************************************************
 * Traits indicate, if objects are construtible.
 *****************************************************************************/


//! By default, an object is json_contructible, if it is default constructible
template <typename T>
struct is_json_constructible {
	constexpr static const bool value = std::is_default_constructible<T>::value;
};


// Disallow some basic types.
template <>
struct is_json_constructible<void> {
	constexpr static const bool value = false;
};

template <>
struct is_json_constructible<float> {
	constexpr static const bool value = false;
};

template <>
struct is_json_constructible<double> {
	constexpr static const bool value = false;
};

template <>
struct is_json_constructible<bool> {
	constexpr static const bool value = false;
};

template <>
struct is_json_constructible<int> {
	constexpr static const bool value = false;
};

template <>
struct is_json_constructible<size_t> {
	constexpr static const bool value = false;
};

/*
template <>
struct is_json_constructible<LuaFunctionHandle> {
	constexpr static const bool value = false;
};
*/


/*****************************************************************************
 * JSON default values.
 *****************************************************************************/

//! Default value (empty).
template <typename T>
struct json_default { static const nlohmann::json value; };

template <typename T>
const nlohmann::json json_default<T>::value = {};

template <typename T>
inline const nlohmann::json json_default_value()
{ return json_default<T>::value; }




/*****************************************************************************
 * JSON schema.
 *****************************************************************************/
//! default schema (empty)
template <typename T>
struct json_schema{
	 static constexpr const char* value = R"(
	{
  		"$schema": "http://json-schema.org/draft-07/schema#",
  		"type": "object",
 		"additionalProperties": false
	}
	)";
};


template <typename T>
inline const char* json_schema_value() {return json_schema<T>::value;}


/*****************************************************************************
 * JSON assignment.
 *****************************************************************************/
//! As a default, objects are assigned using nlohmann::from_json.
template <typename T>
struct json_assignment
{
	static void from_json(const nlohmann::json& j, T& obj)
	{ nlohmann::from_json(j, obj); }
};

/*****************************************************************************
 * JSON builder.
 *****************************************************************************/
//! This class provides basic constructors.
template <typename T>
class JSONBuilder {

public:

	JSONBuilder()
	: m_json(json_default<T>::value) {}


	JSONBuilder(nlohmann::json &j)
	: m_json(j) {}

	virtual ~JSONBuilder() {};

	//! Create a new initialized instance.
	virtual SmartPtr<T> build()
	{ return (class_is_json_constructible()) ? make_sp<T> (this->build_raw()) : SPNULL; }

	std::string get_default()
	{
		return json_default_value<T>().dump();
	}

	const nlohmann::json& get_json_default()
	{ return json_default_value<T>(); }

	const nlohmann::json& get_json()
	{ return m_json; }




protected:

	nlohmann::json m_json;

	// Returns is_json_constructible<T>::value.
	static bool class_is_json_constructible()
	{
		bool check = is_json_constructible<T>::value;
		if (!check) UG_LOG("WARNING: is_json_constructible<T>::value is false.");
		return check;
	}

	// Create new instance & initialize.
	T* build_raw()
	{
		T* obj = (class_is_json_constructible()) ? new T() : NULL;
		if (obj) json_assignment<T>::from_json(m_json, *obj);
		return obj;
	}



};


/*****************************************************************************
 * JSON default definitions.
 *****************************************************************************/
struct json_predefined_defaults {
	static const nlohmann::json solvers;
};

/*****************************************************************************
 * JSON schema definitions.
 *****************************************************************************/
struct json_predefined_schemas {
	static const nlohmann::json solvers;
};



} // namespace ug
