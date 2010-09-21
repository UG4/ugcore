// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m09 d14

#ifndef __H__UG__INTERFACE__INTERFACE_BASE__
#define __H__UG__INTERFACE__INTERFACE_BASE__

#include <vector>
#include <string>
#include <cassert>

namespace ug{
namespace interface{


class IObject;

enum ParameterTypeID{
	PTID_UNKNOWN = 0,
	PTID_INT,
	PTID_DOUBLE,
	PTID_STRING,
	PTID_OBJECT
};

////////////////////////////////////////////////////////////////////////
class IParameter
{
	public:
		virtual ~IParameter()	{}
		
		virtual IParameter* clone() const = 0;
		
		virtual const char* get_type_string()				{return "";}
		virtual int get_type_id()							{return PTID_UNKNOWN;}
		
		virtual int to_int(bool* bOKOut = NULL)				{if(bOKOut) *bOKOut = false; return 0;}
		virtual double to_double(bool* bOKOut = NULL)		{if(bOKOut) *bOKOut = false; return 0;}
		virtual const char* to_string(bool* bOKOut = NULL)	{if(bOKOut) *bOKOut = false; return "";}
		virtual IObject* to_object(bool* bOKOut = NULL)		{if(bOKOut) *bOKOut = false; return NULL;}
		
		virtual bool set_int(int val)				{return false;}
		virtual bool set_double(double val)			{return false;}
		virtual bool set_string(const char* val)	{return false;}
		virtual bool set_object(IObject* val)		{return false;}
		
		void set_name(const char* name)				{m_name = name;}
		const char* get_name() const				{return m_name.c_str();}
		
		void set_representation(const char* rep)	{m_representation = rep;}
		const char* get_representation() const		{return m_representation.c_str();}
		
	private:
	//todo: this could be a const char* most of the time.
		std::string m_name;
	//todo: this could be a const char* most of the time.
		std::string m_representation;
};

////////////////////////////////////////////////////////////////////////
class ParameterList
{
	public:
		ParameterList()	{}
		ParameterList(const ParameterList& pl)
		{
			for(size_t i = 0; i < pl.size(); ++i){
				m_params.push_back(pl.m_params[i]->clone());
				m_interfaceTypes.push_back(pl.m_interfaceTypes[i]);
			}
		}
		
		ParameterList& operator =(const ParameterList& pl)
		{
			clear();
			for(size_t i = 0; i < pl.size(); ++i){
				m_params.push_back(pl.m_params[i]->clone());
				m_interfaceTypes.push_back(pl.m_interfaceTypes[i]);
			}
			return *this;
		}
		
		void clear()
		{
			for(size_t i = 0; i < m_params.size(); ++i)
				delete m_params[i];
			m_params.clear();
			m_interfaceTypes.clear();
		}
		
		
	/**	if min > max, min and max are ignored.
	 *	In order to use the default representation, pass an empty string ("").*/
		void add_int(const char* paramName, int value = 0,
					 const char* representation = "");

	/**	if min > max, min and max are ignored.
	 *	In order to use the default representation, pass an empty string ("").*/
		void add_double(const char* paramName, double value = 0,
						const char* representation = "");

		void add_string(const char* paramName, const char* value = "",
						const char* representation = "");

		void add_object(const char* paramName,
						const char* interfaceType = "",
						IObject* value = NULL,
						const char* representation = "");
		
		//void add_int_list
		//void add_double_list
		//void add_string_list
		//void add_object_list
		
		inline size_t size() const						{return m_params.size();}
				
	///	returns the parameter at the given index.
	/**	Please note that although the method is a const method, the returned
	 *	parameter is not const. This is important to allow a user to manipulate
	 *	the values of parameters of a const ParameterList. (constness regards
	 *	the number and type of parameters only.)*/
		inline IParameter* param(size_t index) const	{return m_params.at(index);}
		
	///	returns the interface type of the parameter at the given index.
		inline const char* get_interface_type(size_t index)	const {return m_interfaceTypes.at(index).c_str();}
		
	////////////////////////////////
	//	CONVENIANCE FUNCTIONS
	///	a conveniance function. Forwards to param(index).to_int().
		inline int to_int(size_t index)				{return param(index)->to_int();}
		
	///	a conveniance function. Forwards to param(index).to_double().
		inline double to_double(size_t index)		{return param(index)->to_double();}

	///	a conveniance function. Forwards to param(index).to_string().
		inline const char* to_string(size_t index)	{return param(index)->to_string();}
		
	///	a conveniance function. Forwards to param(index).to_object().
		inline IObject* to_object(size_t index)		{return param(index)->to_object();}
		
	///	a conveniance function. Forwards to param(index).set_int(val)
		inline void set_int(size_t index, int val)				{param(index)->set_int(val);}

	///	a conveniance function. Forwards to param(index).set_double(val)
		inline void set_double(size_t index, double val)		{param(index)->set_double(val);}
		
	///	a conveniance function. Forwards to param(index).set_int(val)
		inline void set_string(size_t index, const char* val)	{param(index)->set_string(val);}
		
	///	a conveniance function. Forwards to param(index).set_object(val)
		inline void set_object(size_t index, IObject* val)		{param(index)->set_object(val);}
		
	protected:
		std::vector<IParameter*>	m_params;
	//	for each object parameter this vector holds the associated supported type.
	//	For other parameters it holds "". Indices correspond to m_params.
		std::vector<std::string>	m_interfaceTypes;
};

////////////////////////////////////////////////////////////////////////
class MethodDesc
{
	public:
		MethodDesc(const char* name = "",
				   const char* tooltip = "",
				   const char* help = "") :
				m_name(name),
				m_tooltip(tooltip),
				m_help(help),
				m_index(0)	{}

		MethodDesc(const char* name, int index,
				   const char* tooltip = "",
				   const char* help = "") :
				m_name(name),
				m_tooltip(tooltip),
				m_help(help),
				m_index(index)	{}
				
		
	////////////////////////////////
	//	non-const methods
		std::string& name()		{return m_name;}
		std::string& tooltip()	{return m_tooltip;}
		std::string& help()		{return m_help;}
		
		ParameterList& params_in()	{return m_paramsIn;}
		ParameterList& params_out()	{return m_paramsOut;}

	////////////////////////////////
	//	const methods
		size_t index() const					{return m_index;}
		
		const std::string& name() const			{return m_name;}
		const std::string& tooltip() const		{return m_tooltip;}
		const std::string& help() const			{return m_help;}
		
		const ParameterList& params_in() const	{return m_paramsIn;}
		const ParameterList& params_out() const	{return m_paramsOut;}
		
	private:
		ParameterList	m_paramsIn;
		ParameterList	m_paramsOut;
	
		std::string		m_name;
		std::string		m_tooltip;
		std::string		m_help;
		int				m_index;
};

////////////////////////////////////////////////////////////////////////
class IGlobalFunction
{
	public:
		IGlobalFunction()	{}
		IGlobalFunction(const char* name) : m_method(name)	{}
		
		virtual ~IGlobalFunction()	{}
		
	///	executes the global function. Forwards to an overloaded protected version.
		inline void execute()	{execute(m_method.params_in(), m_method.params_out());}
	
	///	default implementation returns an empty string.
		virtual const char* group()			{return "";}
	
	///	returns the method.
		inline const MethodDesc& get_method() const
		{
			return m_method;
		}

	protected:
	///	executes the method for the given MethodDesc.
	/**	this method is pure virtual and has to be overloaded by derived classes.
	 *	Please note that only methods created by the tool itself may be passed
	 *	as parameter.
	 *
	 *	\param in	Passed for conveniance. Should be he same as m_method.params_in().
	 *	\param out	Passed for conveniance. Should be he same as m_method.params_out().
	 */
		virtual void execute(ParameterList& in, ParameterList& out) = 0;
							 
		void set_method(const MethodDesc& method);
	
	protected:
		MethodDesc	m_method;
};

////////////////////////////////////////////////////////////////////////
class IObject
{
	public:
		virtual ~IObject()	{}

	///	clones the object. Has to be implemented by derived classes.
		virtual IObject* clone() = 0;

	///	executes the method at the specified index. Forwards to the protected overload.
		inline void execute(size_t index)
		{
			assert(index < m_methods.size() && "Bad index.");
			MethodDesc& md = m_methods[index];
			execute(index, md.params_in(), md.params_out());
		}
		
	///	returns the objects type_name.
	/**	this method is pure virtual and has to be overloaded by derived classes.*/
		virtual const char* type_name() = 0;

	///	default implementation returns an empty string.
		virtual const char* group()					{return "";}
				
	///	default implementation returns "".
		virtual const char* supported_types()		{return "";}
		
	///	returns the number of methods.
	/**	methods are typically created during construction.*/
		inline size_t num_methods()	const			{return m_methods.size();}
		
	///	returns the method with the given index.
	/**	methods are typically created during construction.
	 *	Make sure that 0 <= index < num_methods().*/
		const MethodDesc& get_method(size_t index) const	{return m_methods.at(index);}
		
	protected:
	///	executes the method for the given MethodDesc.
	/**	this method can be overloaded by derived classes.
	 *	Please note that only methods created by the object itself are passed
	 *	as parameter.
	 *
	 *	\param in	Passed for conveniance. Same as m_methods[methodIndex].params_in().
	 *	\param out	Passed for conveniance. Same as m_methods[methodIndex].params_out().
	 */
		virtual void execute(size_t methodIndex,
							 ParameterList& in, ParameterList& out)	{};
		
		MethodDesc& add_method(const char* name)
		{
			m_methods.push_back(MethodDesc(name));
			return m_methods.back();
		}
		
	private:
		std::vector<MethodDesc>	m_methods;
};

}//	end of namespace interface
}//	end of namespace ug

#endif
