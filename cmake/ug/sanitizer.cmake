# Copyright (c) 2025:  Goethe University Frankfurt
# Author: Arne Naegel
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# included from ug_includes.cmake


# Support for sanitizers (Address, Memory, Undefined, Leak)
# Usage:
# 1) Enable in CMake configure step:
#    -USE_SANITIZER="ASAN;UBSAN;MSAN;LSAN"
#    (MSan cannot be combined with ASan/UBSan/LSan)
# 2) Enable on targets:
#    target_enable_sanitizers(my_lib ASAN UBSAN)
#    target_enable_sanitizers(my_app MSAN)  

# Or, if you prefer Address + Undefined instead of MSan:
# target_enable_sanitizers(my_lib ASAN UBSAN)
# target_enable_sanitizers(my_app ASAN UBSAN)


# Only include once.
include_guard (GLOBAL)

if (NOT DEFINED USE_SANITIZER OR USE_SANITIZER STREQUAL "" OR USE_SANITIZER STREQUAL "OFF")
    message (STATUS "Info: No sanitizers active (DEFAULT).")
    set (USE_SANITIZER OFF)
    # Dummy function
    # (to avoid if(USE_SANITIZER) ... else ... endif () around every call
    function (target_enable_sanitizers target)
    endfunction ()
   
else ()
    message (STATUS "Info: Using sanitizers: ${USE_SANITIZER}")
    
function (_apply_to_target _tgt _cxx_flags _link_flags)
  if (NOT TARGET "${_tgt}")
    message (FATAL_ERROR "target_enable_sanitizers: target '${_tgt}' does not exist")
  endif ()

  # Apply only to Debug builds (change to $<CONFIG:Debug,RelWithDebInfo> if desired)
  target_compile_options ("${_tgt}" PRIVATE
    $<$<CONFIG:Debug>:${_cxx_flags}>
  )
  # CMake 3.13+: target_link_options; for older CMake, fall back to LINK_FLAGS
  if (COMMAND target_link_options)
    target_link_options ("${_tgt}" PRIVATE $<$<CONFIG:Debug>:${_link_flags}> )
  else ()
    # Legacy fallback
    get_target_property (_old "${_tgt}" LINK_FLAGS)
    if (NOT _old)
      set (_old "")
    endif ()

    set_target_properties ("${_tgt}" PROPERTIES LINK_FLAGS "${_old} ${_link_flags}" )
  endif ()
endfunction ()

function (target_enable_sanitizers target)
  set (options)
  set (oneValueArgs)
  set (multiValueArgs ASAN MSAN UBSAN LSAN)
  cmake_parse_arguments (SAN "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Normalize requested sanitizers into a list (keys present -> enable)
  set (requested "")
  if (SAN_ASAN OR "ASAN" IN_LIST ARGN)  
    list (APPEND requested "ASAN")
  endif ()

  if (SAN_MSAN OR "MSAN" IN_LIST ARGN)  
    list (APPEND requested "MSAN")
  endif ()

  if (SAN_UBSAN OR "UBSAN" IN_LIST ARGN) 
    list (APPEND requested "UBSAN")
  endif ()

  if (SAN_LSAN OR "LSAN" IN_LIST ARGN) 
    list (APPEND requested "LSAN")
  endif ()

  if (requested STREQUAL "")
    message (FATAL_ERROR "target_enable_sanitizers: specify at least one of ASAN, MSAN, UBSAN, LSAN")
  endif ()

  # Detect compiler/toolchain
  set (is_clang FALSE)
  set (is_gcc   FALSE)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (is_clang TRUE)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (is_gcc TRUE)
  endif ()

  # Build up flags
  set (cxx_flags "")
  set (link_flags "")

  foreach (san IN LISTS requested)
    if (san STREQUAL "MSAN")
      if (NOT is_clang)
        message (FATAL_ERROR "MemorySanitizer (MSan) requires Clang.")
      endif ()

      if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
        message (FATAL_ERROR "MemorySanitizer is officially supported on Linux; current: ${CMAKE_SYSTEM_NAME}")
      endif ()

      # MSan specifics
      list (APPEND cxx_flags  "-fsanitize=memory" "-fsanitize-memory-track-origins=2" "-fno-omit-frame-pointer" "-fno-optimize-sibling-calls")
      list (APPEND link_flags "-fsanitize=memory")
    elseif (san STREQUAL "ASAN")
      # AddressSanitizer
      list (APPEND cxx_flags  "-fsanitize=address" "-fno-omit-frame-pointer")
      list (APPEND link_flags "-fsanitize=address")
    elseif (san STREQUAL "UBSAN")
      # UndefinedBehaviorSanitizer (trap optional: -fno-sanitize-recover=all)
      list (APPEND cxx_flags  "-fsanitize=undefined")
      list (APPEND link_flags "-fsanitize=undefined")
      # GCC sometimes needs this to keep type info for UBSan reports nicer:
      if (is_gcc)
        list (APPEND cxx_flags "-fno-omit-frame-pointer")
      endif ()

    elseif (san STREQUAL "LSAN")
      # LeakSanitizer (included in ASan on Clang by default, explicit on GCC)
      list (APPEND cxx_flags  "-fsanitize=leak")
      list (APPEND link_flags "-fsanitize=leak")
    endif ()

  endforeach ()

  # Guard against incompatible combos (MSan vs others)
  if ("MSAN" IN_LIST requested AND
      ("ASAN" IN_LIST requested OR "UBSAN" IN_LIST requested OR "LSAN" IN_LIST requested))
    message (FATAL_ERROR "MSan cannot be combined with ASan/UBSan/LSan. Choose MSAN alone.")
  endif ()

  # Apply to target
  message (STATUS "Sanitizers for target '${target}': ${requested}")
  message (STATUS "Sanitizer CXX-compile flags: ${cxx_flags}")
  message (STATUS "Sanitizer CXX-linker flags: ${link_flags}")
  # _apply_to_target(${target} "${cxx_flags}" "${link_flags}")


  add_cxx_flags ("${cxx_flags}") # ug4 function
  target_link_options (${target} PRIVATE ${link_flags})
  # Nice-to-have: disable in Release unless user forces it
  
endfunction ()


endif ()
