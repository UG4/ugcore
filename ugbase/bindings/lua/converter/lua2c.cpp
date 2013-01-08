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
    parser.createC(cout);
	iIn = parser.num_in();
	iOut = parser.num_out();
	
	out.close();
	//system("cat lua2c_output.c");

	string c1s=string("gcc -O3 -c ") + p + "lua2c_output.c -o " + p + "lua2c_output.o";
	const char *c=c1s.c_str();
	//UG_LOG(c << "\n");
	system(c);


	name = functionName;
	string pDyntemplate = string(functionName) + ".dylib_XXXXXX";
	strcpy(b, pDyntemplate.c_str());
	pDyn = p+mktemp(b);
	string c2s=string("gcc -dynamiclib ") + p+"lua2c_output.o -o " + pDyn.c_str();
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
       	//UG_LOG("removing " << name << "\n");
		dlclose(libHandle);	
		string s = string("rm ") + pDyn;
		UG_LOG(s << "\n");
		system(s.c_str());
	}
}


}
}