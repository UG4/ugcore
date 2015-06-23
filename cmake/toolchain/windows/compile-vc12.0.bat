@set CWD="%~dp0"

set path=%path:"=%

@SET VSINSTALLDIR=%Program Files(x86)%\Microsoft Visual Studio 12.0\
@SET VCINSTALLDIR=%Program Files(x86)%\Microsoft Visual Studio 12.0\VC\

CALL "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat" x86_amd64

cd "%CWD%"

cd "%WORKSPACE%"
mkdir build
cd build

cmake ..\trunk -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release -DTARGET=vrl -DLAPACK=OFF -DBLAS=OFF -DINTERNAL_BOOST=ON -DEMBEDDED_PLUGINS=ON -DCOMPILE_INFO=OFF -DSmallStrainMechanics=OFF -DReceptorKinetic=OFF -Dd3f=ON -DLuaShell=ON -DProMesh=ON -DCMAKE_CXX_FLAGS="/MD /O2 /Ob2 /D NDEBUG /bigobj" -G"NMake Makefiles"

nmake
