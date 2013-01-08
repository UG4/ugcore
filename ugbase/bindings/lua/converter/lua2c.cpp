#include "pclass.h"
#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"
#include <dlfcn.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "lua2c.h"
using namespace std;

namespace ug{
namespace bridge {
    
DebugID DID_LUA2C("LUA2C");

bool LUA2C::create(const char *functionName)
{
	f=NULL;
	pclass parser;	
	if(parser.parse_luaFunction(functionName) == false) return 0;
	//parser.reduce();

	char cwd[1024];
	getcwd(cwd, 1024);

	string p = PathProvider::get_path(ROOT_PATH) + "/bin/lua2c_tmp/";
    
    if(!DirectoryExists(p.c_str()))
        CreateDirectory(p.c_str());
    
	chdir(p.c_str());

	// /Users/mrupp/NetBeansProjects/test/

	char b[255];

	fstream out("lua2c_output.c", fstream::out);		

	out << "#include <math.h>\n";
	out << "#define true 1\n";
	out << "#define false 0\n";
	parser.createC(out);
    
    iIn = parser.num_in();
	iOut = parser.num_out();
    out.close();
    
    IF_DEBUG(DID_LUA2C, 1)
    {    
        UG_LOG("--[LUA2C]-----------------------------------------------\n");
        UG_LOG("created C function from LUA function " << functionName << ":\n");
        UG_LOG(GetFileLines("lua2c_output.c", 1, -1, true) << "\n");        
		UG_LOG("--[LUA2C]-----------------------------------------------\n");
    }
	

	string c1s=string("gcc -O3 -c ") + p + "lua2c_output.c -o " + p + "lua2c_output.o";
	UG_DLOG(DID_LUA2C, 2, "compiling line: " << c1s);
	const char *c=c1s.c_str();
	//UG_LOG(c << "\n");
	system(c);


	name = functionName;
	string pDyntemplate = string(functionName) + ".dylib_XXXXXX";
	strcpy(b, pDyntemplate.c_str());
	pDyn = p+mktemp(b);
	string c2s=string("gcc -dynamiclib ") + p+"lua2c_output.o -o " + pDyn.c_str();
	UG_DLOG(DID_LUA2C, 2, "linking line: " << c2s);
	const char *c2 = c2s.c_str();
	//UG_LOG(c2 << "\n");
	system(c2);	

	//UG_LOG("pDyn = " << pDyn << "\n");

	libHandle = dlopen(pDyn.c_str(), RTLD_LAZY);

	if(libHandle!=NULL)
        f = (LUA2C_Function) dlsym(libHandle, functionName);
	
	chdir(cwd);
	//parser.createLUA(cout);	
	return f != NULL;
}

LUA2C::~LUA2C()
{
	if(libHandle)
	{
       	UG_DLOG(DID_LUA2C, 2, "removing " << name << "\n");
		dlclose(libHandle);	
		string s = string("rm ") + pDyn;
		UG_DLOG(DID_LUA2C, 2, s << "\n");
		system(s.c_str());
	}
}


}
}