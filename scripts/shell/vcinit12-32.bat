REM Fixes path issues. Unfortunately, MSVC scripts don't work with quoted paths.
set path=%path:"=%

call "%VS120COMNTOOLS%..\..\VC\vcvarsall.bat"
