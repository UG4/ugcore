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

#include <nlohmann/json.hpp>

#include "common/util/smart_pointer.h"


/*
 * For each class, we must supply the following:
 *
 * - is_json_constructible = true
 * - json_defaults = <valid JSON string>
 * - json_schema = <valid JSON schema >
 *
 * Optional defines:
 *
 * - json_assignment::from_json
 * - json_factory::create_new_empty, create_new, create_smart_ptr
 * */

namespace ug {

//! Marks an object, with is constructible via JSON.
/*! DEPRECATED */
struct JSONConstructible {};



//! By default, an object is json_contructible, if it is default constructible
template <typename T>
struct is_json_constructible {
	constexpr static const bool value = std::is_default_constructible<T>::value;
};

/*****************************************************************************
 * JSON default values.
 *****************************************************************************/
//! default value (empty)
template <typename T>
struct json_default{
	 static constexpr const char* value = "{}"; // empty value
//	 static const char* get_value() { return value; }
};

template <typename T>
inline const char* json_default_value() {return json_default<T>::value;}


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
	// static const char* get_value() { return value; }
};


template <typename T>
inline const char* json_schema_value() {return json_schema<T>::value;}



// Disallow some basic types
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

template <>
struct is_json_constructible<LuaFunctionHandle> {
	constexpr static const bool value = false;
};


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
 * JSON factory helpers.
 *****************************************************************************/
//! This class provides basic constructors.
template <typename T>
struct json_factory {

	//! Create a new un-initialized instance.
	static T* create_new_empty()
	{
		return (is_json_constructible<T>::value) ? new T() : NULL;
	}

	//! Create a new initialized instance.
	static T* create_new(const nlohmann::json &j)
	{
		if (! is_json_constructible<T>::value) return NULL;

		T* obj=create_new_empty();
		if (obj) json_assignment<T>::from_json(j, *obj);
		return obj;

	}

	static T* create_new_default()
	{

	}

	//! Create a new initialized instance.
	static SmartPtr<T> create_smart_ptr(const nlohmann::json &j)
	{ return (is_json_constructible<T>::value) ? make_sp<T> (json_factory<T>::create_new(j)) : SPNULL; }

};



}
