REM ########################################################
REM # Batch script that compiles UG4 with MINGW or TDM-GCC #
REM # Author: Michael Hoffer <info@michaelhoffer.de>       #
REM ########################################################

REM ---------------------------------------------------------------------------
REM - SETUP MINGW COMPILER                                                    -
REM ---------------------------------------------------------------------------

REM set current directory
@set CWD="%cd%"

REM Set MINGW to the MinGW compiler that shall be used
set MINGW="C:\TDM-GCC-64"
set PATH=%MINGW%\bin;%path%
set path=%path:"=%

REM ---------------------------------------------------------------------------
REM - BUILD UG4                                                               -
REM ---------------------------------------------------------------------------

cd "%WORKSPACE%"
rmdir /s /q build
mkdir build
cd build

cmake ..\trunk -DDEBUG=OFF -DCMAKE_BUILD_TYPE=Release -DTARGET=vrl -DLAPACK=OFF -DBLAS=OFF -DINTERNAL_BOOST=ON -DEMBEDDED_PLUGINS=ON -DCOMPILE_INFO=OFF -DSmallStrainMechanics=OFF -DReceptorKinetic=OFF -Dd3f=ON -DLuaShell=ON -DProMesh=ON -G"MinGW Makefiles"

mingw32-make