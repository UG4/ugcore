# Install script for directory: /home/lreinhardt/eclipse/cpp/ug4/ugbase

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/compile_info/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/common/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/lib_grid/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/registry/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/bindings/lua/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/lib_algebra/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/lib_disc/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/pcl/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/bridge/cmake_install.cmake")
  INCLUDE("/home/lreinhardt/eclipse/cpp/ug4/cmake/ugbase/ug_shell/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

