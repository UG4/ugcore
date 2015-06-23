set CWD="%~dp0"

cd "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC"
vcvarsall.bat x86_amd64

cd "%CWD%"

cd "%WORKSPACE%"
mkdir build
cd build

cmake ..\trunk -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release -DTARGET=vrl -DLAPACK=OFF -DBLAS=OFF -DINTERNAL_BOOST=ON DEMBEDDED_PLUGINS=ON -DCOMPILE_INFO=OFF -DSmallStrainMechanics=OFF -DReceptorKinetic=OFF -Dd3f=ON -DLuaShell=ON -DProMesh=ON -G "NMake Makefiles"
nmake

