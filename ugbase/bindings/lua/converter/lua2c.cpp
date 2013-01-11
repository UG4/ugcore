#include "pclass.h"
#include "bridge/bridge.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_stack_check.h"
#include "bindings/lua/info_commands.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "lua2c.h"
using namespace std;

namespace ug{
namespace bridge {
    
DebugID DID_LUA2C("LUA2C");

string MakeTmpFile(string filename, string extension, bool &bSuccess)
{
	bSuccess = true;
	const char *name = (filename+extension).c_str();
	if(!FileExists(name)) return name;
	for(int i=0;i<999999; i++)
	{
		stringstream ss;
		ss << filename << i << extension;
		name = ss.str().c_str();
		if(!FileExists(name)) return name;
	}
	bSuccess = false;
	return "";
}

bool LUA2C::create(const char *functionName)
{
	m_f=NULL;
	pclass parser;	
	if(parser.parse_luaFunction(functionName) == false) return 0;
	//parser.reduce();

	string p = PathProvider::get_path(ROOT_PATH) + "/bin/lua2c_tmp/";
    
    if(!DirectoryExists(p.c_str()))
        CreateDirectory(p.c_str());
    
	fstream out((p + "lua2c_output.c").c_str(), fstream::out);		

	out << "#include <math.h>\n";
	out << "#define true 1\n";
	out << "#define false 0\n";
	parser.createC(out);
    
    m_iIn = parser.num_in();
	m_iOut = parser.num_out();
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
	system(c1s.c_str());

	bool bTmpFileSuccess=false;
	m_name = functionName;
	
	m_pDyn = MakeTmpFile(p+string(functionName), ".dylib", bTmpFileSuccess);
	
	string c2s=string("gcc -dynamiclib ") + p+"lua2c_output.o -o " + m_pDyn.c_str();
	UG_DLOG(DID_LUA2C, 2, "linking line: " << c2s);
	system(c2s.c_str());	

	m_libHandle = OpenLibrary(m_pDyn.c_str());

	if(m_libHandle!=NULL)
        m_f = (LUA2C_Function) GetLibraryProcedure(m_libHandle, functionName);
	
	return m_f != NULL;
}

LUA2C::~LUA2C()
{
    UG_DLOG(DID_LUA2C, 2, "removing " << m_name << "\n");		
	if(m_libHandle)
	   	CloseLibrary(m_libHandle);
        
    if(m_pDyn.size() > 0)
    {
		string s = string("rm ") + m_pDyn;
		UG_DLOG(DID_LUA2C, 2, s << "\n");
		system(s.c_str());
	}
}


}
}
