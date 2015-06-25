REM ##################################################
REM # Batch script that compiles UG4 with MSVC 12.0  #
REM # Author: Michael Hoffer <info@michaelhoffer.de> #
REM ##################################################

REM ---------------------------------------------------------------------------
REM - SETUP MSVC COMPILER                                                     -
REM ---------------------------------------------------------------------------

REM set current directory
@set CWD="%cd%"

REM Fixes path issues. Unfortunately, MSVC scripts don't work with quoted paths.
set path=%path:"=%

REM Ensures that we compile with MSVC 12.0 for amd64 target
CALL "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat" x86_amd64

REM cd "%CWD%"

REM ---------------------------------------------------------------------------
REM - BUILD UG4                                                               -
REM ---------------------------------------------------------------------------

cd "%WORKSPACE%"
rmdir /s /q build
mkdir build
cd build

cmake ..\trunk -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release -DTARGET=vrl -DLAPACK=OFF -DBLAS=OFF -DINTERNAL_BOOST=ON -DEMBEDDED_PLUGINS=ON -DCOMPILE_INFO=OFF -DSmallStrainMechanics=OFF -DReceptorKinetic=OFF -Dd3f=ON -DLuaShell=ON -DProMesh=ON -DCMAKE_CXX_FLAGS="/MD /O2 /Ob2 /D NDEBUG /bigobj" -G"NMake Makefiles"

nmake


