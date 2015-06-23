echo "%JAVA_HOME%"

echo "%PATH%"

set CWD="%~dp0"

@set PATH="C:\Windows\SysWOW64;C:\ProgramData\Oracle\Java\javapath;C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;C:\Windows\System32\WindowsPowerShell\v1.0\;C:\Program Files (x86)\CMake 2.8\bin;C:\Program Files (x86)\Subversion\bin;C:\Program Files\mingw-w64\x86_64-4.8.4-posix-seh-rt_v3-rev0\mingw64\bin;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files\Microsoft SQL Server\120\Tools\Binn\;"C:\Program Files (x86)\Java\jre1.8.0_45\bin"" 
12:14:37 "C:\Windows\SysWOW64;C:\ProgramData\Oracle\Java\javapath;C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;C:\Windows\System32\WindowsPowerShell\v1.0\;C:\Program Files (x86)\CMake 2.8\bin;C:\Program Files (x86)\Subversion\bin;C:\Program Files\mingw-w64\x86_64-4.8.4-posix-seh-rt_v3-rev0\mingw64\bin;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files\Microsoft SQL Server\120\Tools\Binn\;"

@SET VSINSTALLDIR=%ProgramFiles(x86)%\Microsoft Visual Studio 12.0\
@SET VCINSTALLDIR=%ProgramFiles(x86)%\Microsoft Visual Studio 12.0\VC\

cd "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC"
vcvarsall.bat x86_amd64

cd "%CWD%"

cd "%WORKSPACE%"
mkdir build
cd build

cmake ..\trunk -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release -DTARGET=vrl -DLAPACK=OFF -DBLAS=OFF -DINTERNAL_BOOST=ON DEMBEDDED_PLUGINS=ON -DCOMPILE_INFO=OFF -DSmallStrainMechanics=OFF -DReceptorKinetic=OFF -Dd3f=ON -DLuaShell=ON -DProMesh=ON -G "NMake Makefiles"
nmake
