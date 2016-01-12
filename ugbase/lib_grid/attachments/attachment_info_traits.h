/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_attachment_info_traits
#define __H__UG_attachment_info_traits

#include <typeinfo>
#include "lib_grid/common_attachments.h"

namespace ug{

#define DECLARE_ATTACHMENT_INFO_TRAITS(attachmentType, typeName)\
		template <> struct attachment_info_traits<attachmentType> {\
			static const char* type_name ()	{return typeName;}};

template <class TAttachment>
struct attachment_info_traits {
	static const char* type_name ();
};

DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<bool>, "bool");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<char>, "char");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<byte>, "byte");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<int>, "int");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<uint>, "uint");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<float>, "number");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<double>, "number");

DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<vector1>, "vector1");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<vector2>, "vector2");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<vector3>, "vector3");
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<vector4>, "vector4");


}//	end of namespace

#endif	//__H__UG_attachment_info_traits
