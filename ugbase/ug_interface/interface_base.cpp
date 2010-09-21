// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m09 d20

#include "interface_base.h"

using namespace std;

namespace ug{
namespace interface{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	the following section contains the different parameter-implementations

///	A helper class that simplifies parameter-class implementation
/**	\tparam TVal		the value-type.
 *	\tparam strType		a string that specifies the type-name.
 *	\tparam TCloneType	instances of this type will be returned by clone().
 *						Please note that TCloneType has to be derived from IParameter.
 */
template <class TVal, class TCloneType, int ParamTypeID>
class TParameter : public IParameter
{
	public:
		virtual IParameter* clone() const
		{
			TCloneType* param = new TCloneType;
			param->set_name(get_name());
			param->set_representation(get_representation());
			param->template m_val = m_val;
			return param;
		}
		
		virtual int get_type_id()
		{
			return ParamTypeID;
		}
		
	protected:
		TVal	m_val;
};

////////////////////////////////////////////////////////////////////////
class IntParameter : public TParameter<int, IntParameter, PTID_INT>
{
	public:
		virtual const char* get_type_string()	{return "int";}
		
		virtual int to_int(bool* bOKOut = NULL){
			if(bOKOut)
				*bOKOut = true;
			return m_val;
		}
		
		virtual bool set_int(int val){
			m_val = val;
			return true;
		}
};

////////////////////////////////////////////////////////////////////////
class DoubleParameter : public TParameter<double, DoubleParameter, PTID_DOUBLE>
{
	public:
		virtual const char* get_type_string()	{return "double";}
		
		virtual double to_double(bool* bOKOut = NULL){
			if(bOKOut)
				*bOKOut = true;
			return m_val;
		}
		
		virtual bool set_double(double val){
			m_val = val;
			return true;
		}
};

////////////////////////////////////////////////////////////////////////
class StringParameter : public TParameter<string, StringParameter, PTID_STRING>
{
	public:
		virtual const char* get_type_string()	{return "string";}
		
		virtual const char* to_string(bool* bOKOut = NULL){
			if(bOKOut)
				*bOKOut = true;
			return m_val.c_str();
		}
		
		virtual bool set_string(const char* val){
			m_val = val;
			return true;
		}
};

////////////////////////////////////////////////////////////////////////
class ObjectParameter : public TParameter<IObject*, ObjectParameter, PTID_OBJECT>
{
	public:
		virtual const char* get_type_string()	{return "object";}
		
		virtual IObject* to_object(bool* bOKOut = NULL){
			if(bOKOut)
				*bOKOut = true;
			return m_val;
		}
		
		virtual bool set_object(IObject* val){
			m_val = val;
			return true;
		}
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	the following section contains method-implementations for ParameterList.
void ParameterList::
add_int(const char* paramName, int value,
		 const char* representation)
{
	IntParameter* param = new IntParameter;
	param->set_int(value);
	param->set_name(paramName);
	param->set_representation(representation);
	
	m_interfaceTypes.push_back(string(""));
	m_params.push_back(param);
}

/**	if min > max, min and max are ignored.
*	In order to use the default representation, pass an empty string ("").*/
void ParameterList::
add_double(const char* paramName, double value,
			const char* representation)
{
	DoubleParameter* param = new DoubleParameter;
	param->set_double(value);
	param->set_name(paramName);
	param->set_representation(representation);
	
	m_interfaceTypes.push_back(string(""));
	m_params.push_back(param);
}

void ParameterList::
add_string(const char* paramName, const char* value,
			const char* representation)
{
	StringParameter* param = new StringParameter;
	param->set_string(value);
	param->set_name(paramName);
	param->set_representation(representation);
	
	m_interfaceTypes.push_back(string(""));
	m_params.push_back(param);
}

void ParameterList::
add_object(const char* paramName,
			const char* interfaceType,
			IObject* value,
			const char* representation)
{
	ObjectParameter* param = new ObjectParameter;
	param->set_object(value);
	param->set_name(paramName);
	param->set_representation(representation);
	
	m_interfaceTypes.push_back(string(interfaceType));
	m_params.push_back(param);
}

}//	end of namespace interface
}//	end of namespace ug

