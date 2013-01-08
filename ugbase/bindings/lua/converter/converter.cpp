#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"
#include "lua2c.h"
#include "pclass.h"

using namespace std;
namespace ug{
namespace bridge {
	

int convert(const char *functionName)
{
	/*LUA2C f;
	f.create(functionName);
	if(f.is_valid())
	{
		double dOut[1];
		double dIn[1] = {42};
		f.f(dOut, dIn);
		cout << dOut[0] << "\n";
	}*/
	
	//bridge::Registry &reg = GetUGRegistry();
	//reg.add_function(string(functionName)+"_C", f->f, "");
	//reg.registry_changed();

    pclass parser;
    if(parser.parse_luaFunction(functionName) == false) return 0;
    
    parser.createC(cout);

	return 0;	
}


bool RegisterConverter(Registry &reg, const char* parentGroup)
{
	stringstream grpSS; grpSS << parentGroup << "/Converter";
	std::string grp = grpSS.str();

	//try
	{
		reg.add_function("convert", &convert, grp.c_str());		
	}
	//UG_REGISTRY_CATCH_THROW(grp);

	return true;
}






}}