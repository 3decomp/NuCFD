## CMakeLists.txt
#
## Description
#
# Main CMakeLists file for the NuCFD project.
#
## LICENSE
#
# SPDX-License-Identifier: BSD-3-Clause
#

## Basic configuration
cmake_minimum_required(VERSION 3.25.2)
project(NuCFD)
enable_language(Fortran)

option(BUILD_SHARED_LIBS "Enable building as a shared library" OFF)
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE true)
endif()

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

## Set compiler flags using interface libraries
set(gcc_like_fc "$<COMPILE_LANG_AND_ID:Fortran,GNU>")
set(icc_like_fc "$<COMPILE_LANG_AND_ID:Fortran,Intel>")
set(cray_like_fc "$<COMPILE_LANG_AND_ID:Fortran,Cray>")
add_library(nucfd_compiler_flags INTERFACE)
target_compile_options(nucfd_compiler_flags INTERFACE
  "$<${gcc_like_fc}:-cpp;-std=f2018;-Wall;-Wextra;-Wpedantic;-Warray-bounds;-Wimplicit-interface;-Wimplicit-procedure;-fimplicit-none>"
  "$<${icc_like_fc}:-fpp;-std18;-warn;all>"
  "$<${cray_like_fc}:-cpp;-hfp0;-en;-emf;-eI>") # Anything greater than fp0 fails tests!
add_library(nucfd_debug_flags INTERFACE)
target_compile_options(nucfd_debug_flags INTERFACE
  "$<${gcc_like_fc}:-g;-Og;-ffpe-trap=invalid,zero;-fcheck=bounds;-fbacktrace>"
  "$<${icc_like_fc}:-g;-O0;-check;all;-fpe0;-traceback;-init:snan,arrays>"
  "$<${cray_like_fc}:-G2;-O0;-R;bc;-ei>")
add_library(nucfd_strict_flags INTERFACE)
target_compile_options(nucfd_strict_flags INTERFACE
  "$<${gcc_like_fc}:-Werror;-Wno-error=integer-division>"
  "$<${cray_like_fc}:-eN>")
if (CMAKE_BUILD_TYPE MATCHES "Debug")
  target_link_libraries(nucfd_compiler_flags nucfd_debug_flags)
endif()

## Build nucfd library
add_subdirectory(src)

## Testing
enable_testing()
add_subdirectory(tests)

## Documentation
add_custom_target(doc ford
  ${CMAKE_SOURCE_DIR}/NuCFD.md)
