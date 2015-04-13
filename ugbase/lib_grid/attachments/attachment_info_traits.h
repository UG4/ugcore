// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_attachment_info_traits
#define __H__UG_attachment_info_traits

#include <typeinfo>

namespace ug{

#define DECLARE_ATTACHMENT_INFO_TRAITS(attachmentType, typeName)\
		template <> struct attachment_info_traits<attachmentType> {\
			static std::string type_name ()	{return typeName;}};

template <class TAttachment>
struct attachment_info_traits {
	static std::string type_name ();
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
