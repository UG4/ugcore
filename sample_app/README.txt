How to compile sample_app:
- a recent version of cmake is required (min v 2.6), see www.cmake.org
- start in the root folder of the app (.../ug4/sample_app)
- create a new folder and call it build (mkdir build)
- change path to the build folder (cd build)
- execute cmake on the sample_app path (cmake ../)
- execute make (make)
- you now should find the executable in .../ug4/sample_app/bin
- if you want to recompile the executable later on, you wont have to
  execute cmake again. make will automatically recognise changes
  to the source files and to CMakeLists.txt.
